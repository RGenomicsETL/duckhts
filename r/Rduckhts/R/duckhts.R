#' DuckDB HTS File Reader Extension for R
#'
#' The Rduckhts package provides an interface to the DuckDB HTS (High Throughput Sequencing)
#' file reader extension from within R. It enables reading common bioinformatics file formats
#' such as VCF/BCF, SAM/BAM/CRAM, FASTA, FASTQ, GFF, GTF, and tabix-indexed files
#' directly from R using SQL queries via DuckDB.
#'
#' @docType package
#' @name Rduckhts
#' @author DuckHTS Contributors
#' @references \url{https://github.com/RGenomicsETL/duckhts}
#' @keywords internal
#' @importFrom DBI dbExecute dbExistsTable dbRemoveTable
"_PACKAGE"

build_param_str <- function(params) {
  if (length(params) == 0) {
    return("")
  }
  param_str <- paste(paste0(names(params), " := ", params), collapse = ", ")
  if (nchar(param_str) > 0) {
    return(paste0(", ", param_str))
  }
  ""
}

#' Load DuckHTS Extension
#'
#' Loads the DuckHTS extension into a DuckDB connection. This must be called
#' before using any of the HTS reader functions.
#'
#' @param con A DuckDB connection object
#' @param extension_path Optional path to the duckhts extension file. If NULL,
#'   will try to use the bundled extension.
#'
#' @return TRUE if the extension was loaded successfully
#'
#' @examples
#' \dontrun{
#' library(DBI)
#' library(duckdb)
#'
#' con <- dbConnect(duckdb::duckdb())
#' rduckhts_load(con)
#' }
#'
#' @export
rduckhts_load <- function(con, extension_path = NULL) {
  if (is.null(extension_path)) {
    # Try to use the bundled extension
    extension_path <- system.file(
      "extdata",
      "duckhts.duckdb_extension",
      package = "Rduckhts",
      mustWork = FALSE
    )
  }

  if (!file.exists(extension_path)) {
    stop(
      "DuckHTS extension not found at: ",
      extension_path,
      "\nThis suggests the package was not built correctly during installation.",
      "\nTry recompiling the package with R CMD INSTALL --preclean Rduckhts"
    )
  }

  # Enable unsigned extensions if needed
  DBI::dbExecute(con, "SET enable_progress_bar = false")

  # Load the extension
  result <- DBI::dbExecute(con, sprintf("LOAD '%s'", extension_path))
  return(result == 0)
}

#' Create VCF/BCF Table
#'
#' Creates a DuckDB table from a VCF or BCF file using the DuckHTS extension.
#' This follows the RBCFTools pattern of creating a table that can be queried.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param table_name Name for the created table
#' @param path Path to the VCF/BCF file
#' @param region Optional genomic region (e.g., "chr1:1000-2000")
#' @param tidy_format Logical. If TRUE, FORMAT columns are returned in tidy format
#' @param overwrite Logical. If TRUE, overwrites existing table
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' con <- dbConnect(duckdb::duckdb())
#' rduckhts_load(con)
#' rduckhts_bcf(con, "variants", "file.vcf.gz")
#' dbGetQuery(con, "SELECT * FROM variants WHERE QUAL > 100 LIMIT 10")
#' }
#'
#' @export
rduckhts_bcf <- function(
  con,
  table_name,
  path,
  region = NULL,
  tidy_format = FALSE,
  overwrite = FALSE
) {
  if (!missing(table_name) && !is.null(table_name)) {
    if (DBI::dbExistsTable(con, table_name) && !overwrite) {
      stop(
        "Table '",
        table_name,
        "' already exists. Use overwrite = TRUE to replace it."
      )
    }
    if (DBI::dbExistsTable(con, table_name)) {
      DBI::dbRemoveTable(con, table_name)
    }
  }

  # Build the CREATE TABLE query
  params <- list()
  if (!is.null(region)) {
    params$region <- sprintf("'%s'", region)
  }
  if (tidy_format) {
    params$tidy_format <- "true"
  }

  param_str <- build_param_str(params)

  if (!is.null(table_name)) {
    create_query <- sprintf(
      "CREATE TABLE %s AS SELECT * FROM read_bcf('%s'%s)",
      table_name,
      path,
      param_str
    )
  } else {
    create_query <- sprintf(
      "CREATE VIEW bcf_data AS SELECT * FROM read_bcf('%s'%s)",
      path,
      param_str
    )
  }

  DBI::dbExecute(con, create_query)
  invisible(TRUE)
}

#' Create SAM/BAM/CRAM Table
#'
#' Creates a DuckDB table from SAM, BAM, or CRAM files using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param table_name Name for the created table
#' @param path Path to the SAM/BAM/CRAM file
#' @param region Optional genomic region (e.g., "chr1:1000-2000")
#' @param reference Optional reference file path for CRAM files
#' @param overwrite Logical. If TRUE, overwrites existing table
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' con <- dbConnect(duckdb::duckdb())
#' rduckhts_load(con)
#' rduckhts_bam(con, "reads", "file.bam")
#' dbGetQuery(con, "SELECT COUNT(*) FROM reads WHERE FLAG & 4 = 0")
#' }
#'
#' @export
rduckhts_bam <- function(
  con,
  table_name,
  path,
  region = NULL,
  reference = NULL,
  overwrite = FALSE
) {
  if (!missing(table_name) && !is.null(table_name)) {
    if (DBI::dbExistsTable(con, table_name) && !overwrite) {
      stop(
        "Table '",
        table_name,
        "' already exists. Use overwrite = TRUE to replace it."
      )
    }
    if (DBI::dbExistsTable(con, table_name)) {
      DBI::dbRemoveTable(con, table_name)
    }
  }

  params <- list()
  if (!is.null(region)) {
    params$region <- sprintf("'%s'", region)
  }
  if (!is.null(reference)) {
    params$reference <- sprintf("'%s'", reference)
  }

  param_str <- build_param_str(params)

  if (!is.null(table_name)) {
    create_query <- sprintf(
      "CREATE TABLE %s AS SELECT * FROM read_bam('%s'%s)",
      table_name,
      path,
      param_str
    )
  } else {
    create_query <- sprintf(
      "CREATE VIEW bam_data AS SELECT * FROM read_bam('%s'%s)",
      path,
      param_str
    )
  }

  DBI::dbExecute(con, create_query)
  invisible(TRUE)
}

#' Create FASTA Table
#'
#' Creates a DuckDB table from FASTA files using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param table_name Name for the created table
#' @param path Path to the FASTA file
#' @param overwrite Logical. If TRUE, overwrites existing table
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' con <- dbConnect(duckdb::duckdb())
#' rduckhts_load(con)
#' rduckhts_fasta(con, "sequences", "file.fa")
#' dbGetQuery(con, "SELECT NAME, length(SEQUENCE) as len FROM sequences")
#' }
#'
#' @export
rduckhts_fasta <- function(con, table_name, path, overwrite = FALSE) {
  if (!missing(table_name) && !is.null(table_name)) {
    if (DBI::dbExistsTable(con, table_name) && !overwrite) {
      stop(
        "Table '",
        table_name,
        "' already exists. Use overwrite = TRUE to replace it."
      )
    }
    if (DBI::dbExistsTable(con, table_name)) {
      DBI::dbRemoveTable(con, table_name)
    }
  }

  if (!is.null(table_name)) {
    create_query <- sprintf(
      "CREATE TABLE %s AS SELECT * FROM read_fasta('%s')",
      table_name,
      path
    )
  } else {
    create_query <- sprintf(
      "CREATE VIEW fasta_data AS SELECT * FROM read_fasta('%s')",
      path
    )
  }

  DBI::dbExecute(con, create_query)
  invisible(TRUE)
}

#' Create FASTQ Table
#'
#' Creates a DuckDB table from FASTQ files using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param table_name Name for the created table
#' @param path Path to the FASTQ file
#' @param mate_path Optional path to mate file for paired reads
#' @param interleaved Logical indicating if file is interleaved paired reads
#' @param overwrite Logical. If TRUE, overwrites existing table
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' con <- dbConnect(duckdb::duckdb())
#' rduckhts_load(con)
#' rduckhts_fastq(con, "reads", "file.fq")
#' dbGetQuery(con, "SELECT COUNT(*) FROM reads")
#' }
#'
#' @export
rduckhts_fastq <- function(
  con,
  table_name,
  path,
  mate_path = NULL,
  interleaved = FALSE,
  overwrite = FALSE
) {
  if (!missing(table_name) && !is.null(table_name)) {
    if (DBI::dbExistsTable(con, table_name) && !overwrite) {
      stop(
        "Table '",
        table_name,
        "' already exists. Use overwrite = TRUE to replace it."
      )
    }
    if (DBI::dbExistsTable(con, table_name)) {
      DBI::dbRemoveTable(con, table_name)
    }
  }

  params <- list()
  if (!is.null(mate_path)) {
    params$mate_path <- sprintf("'%s'", mate_path)
  }
  if (interleaved) {
    params$interleaved <- "true"
  }

  param_str <- build_param_str(params)

  if (!is.null(table_name)) {
    create_query <- sprintf(
      "CREATE TABLE %s AS SELECT * FROM read_fastq('%s'%s)",
      table_name,
      path,
      param_str
    )
  } else {
    create_query <- sprintf(
      "CREATE VIEW fastq_data AS SELECT * FROM read_fastq('%s'%s)",
      path,
      param_str
    )
  }

  DBI::dbExecute(con, create_query)
  invisible(TRUE)
}

#' Create GFF3 Table
#'
#' Creates a DuckDB table from GFF3 files using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param table_name Name for the created table
#' @param path Path to the GFF3 file
#' @param region Optional genomic region (e.g., "chr1:1000-2000")
#' @param attributes_map Logical. If TRUE, returns attributes as a MAP column
#' @param overwrite Logical. If TRUE, overwrites existing table
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' con <- dbConnect(duckdb::duckdb())
#' rduckhts_load(con)
#' rduckhts_gff(con, "annotations", "file.gff3.gz", attributes_map = TRUE)
#' dbGetQuery(con, "SELECT * FROM annotations WHERE feature = 'gene'")
#' }
#'
#' @export
rduckhts_gff <- function(
  con,
  table_name,
  path,
  region = NULL,
  attributes_map = FALSE,
  overwrite = FALSE
) {
  if (!missing(table_name) && !is.null(table_name)) {
    if (DBI::dbExistsTable(con, table_name) && !overwrite) {
      stop(
        "Table '",
        table_name,
        "' already exists. Use overwrite = TRUE to replace it."
      )
    }
    if (DBI::dbExistsTable(con, table_name)) {
      DBI::dbRemoveTable(con, table_name)
    }
  }

  params <- list()
  if (!is.null(region)) {
    params$region <- sprintf("'%s'", region)
  }
  if (attributes_map) {
    params$attributes_map <- "true"
  }

  param_str <- build_param_str(params)

  if (!is.null(table_name)) {
    create_query <- sprintf(
      "CREATE TABLE %s AS SELECT * FROM read_gff('%s'%s)",
      table_name,
      path,
      param_str
    )
  } else {
    create_query <- sprintf(
      "CREATE VIEW gff_data AS SELECT * FROM read_gff('%s'%s)",
      path,
      param_str
    )
  }

  DBI::dbExecute(con, create_query)
  invisible(TRUE)
}

#' Create GTF Table
#'
#' Creates a DuckDB table from GTF files using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param table_name Name for the created table
#' @param path Path to the GTF file
#' @param region Optional genomic region (e.g., "chr1:1000-2000")
#' @param attributes_map Logical. If TRUE, returns attributes as a MAP column
#' @param overwrite Logical. If TRUE, overwrites existing table
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' con <- dbConnect(duckdb::duckdb())
#' rduckhts_load(con)
#' rduckhts_gtf(con, "annotations", "file.gtf", attributes_map = TRUE)
#' dbGetQuery(con, "SELECT * FROM annotations WHERE feature = 'exon'")
#' }
#'
#' @export
rduckhts_gtf <- function(
  con,
  table_name,
  path,
  region = NULL,
  attributes_map = FALSE,
  overwrite = FALSE
) {
  if (!missing(table_name) && !is.null(table_name)) {
    if (DBI::dbExistsTable(con, table_name) && !overwrite) {
      stop(
        "Table '",
        table_name,
        "' already exists. Use overwrite = TRUE to replace it."
      )
    }
    if (DBI::dbExistsTable(con, table_name)) {
      DBI::dbRemoveTable(con, table_name)
    }
  }

  params <- list()
  if (!is.null(region)) {
    params$region <- sprintf("'%s'", region)
  }
  if (attributes_map) {
    params$attributes_map <- "true"
  }

  param_str <- build_param_str(params)

  if (!is.null(table_name)) {
    create_query <- sprintf(
      "CREATE TABLE %s AS SELECT * FROM read_gtf('%s'%s)",
      table_name,
      path,
      param_str
    )
  } else {
    create_query <- sprintf(
      "CREATE VIEW gtf_data AS SELECT * FROM read_gtf('%s'%s)",
      path,
      param_str
    )
  }

  DBI::dbExecute(con, create_query)
  invisible(TRUE)
}

#' Create Tabix-Indexed File Table
#'
#' Creates a DuckDB table from any tabix-indexed file using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param table_name Name for the created table
#' @param path Path to the tabix-indexed file
#' @param region Optional genomic region (e.g., "chr1:1000-2000")
#' @param overwrite Logical. If TRUE, overwrites existing table
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' con <- dbConnect(duckdb::duckdb())
#' rduckhts_load(con)
#' rduckhts_tabix(con, "bed_data", "file.bed.gz", region = "chr1:1000-2000")
#' dbGetQuery(con, "SELECT * FROM bed_data")
#' }
#'
#' @export
rduckhts_tabix <- function(
  con,
  table_name,
  path,
  region = NULL,
  overwrite = FALSE
) {
  if (!missing(table_name) && !is.null(table_name)) {
    if (DBI::dbExistsTable(con, table_name) && !overwrite) {
      stop(
        "Table '",
        table_name,
        "' already exists. Use overwrite = TRUE to replace it."
      )
    }
    if (DBI::dbExistsTable(con, table_name)) {
      DBI::dbRemoveTable(con, table_name)
    }
  }

  if (!is.null(region)) {
    param_str <- sprintf(", region := '%s'", region)
  } else {
    param_str <- ""
  }

  if (!is.null(table_name)) {
    create_query <- sprintf(
      "CREATE TABLE %s AS SELECT * FROM read_tabix('%s'%s)",
      table_name,
      path,
      param_str
    )
  } else {
    create_query <- sprintf(
      "CREATE VIEW tabix_data AS SELECT * FROM read_tabix('%s'%s)",
      path,
      param_str
    )
  }

  DBI::dbExecute(con, create_query)
  invisible(TRUE)
}
