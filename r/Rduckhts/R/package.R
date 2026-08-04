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
#' @useDynLib Rduckhts, .registration = TRUE
#' @importFrom DBI dbExecute dbExistsTable dbRemoveTable
#' @importFrom duckdb duckdb
"_PACKAGE"
