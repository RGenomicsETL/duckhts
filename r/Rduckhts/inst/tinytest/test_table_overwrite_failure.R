library(DBI)

(function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))

  missing_path <- tempfile(fileext = ".bed")
  expect_false(file.exists(missing_path))

  readers <- list(
    bam = rduckhts_bam,
    bcf = rduckhts_bcf,
    bigwig = rduckhts_bigwig,
    fasta = rduckhts_fasta,
    fastq = rduckhts_fastq,
    genbank = rduckhts_genbank,
    pileup = rduckhts_pileup,
    tabix = rduckhts_tabix,
    bed = rduckhts_bed,
    bam_multi = rduckhts_bam_multi,
    bcf_multi = rduckhts_bcf_multi,
    fastq_multi = rduckhts_fastq_multi,
    fasta_multi = rduckhts_fasta_multi,
    bed_multi = rduckhts_bed_multi,
    tabix_multi = rduckhts_tabix_multi,
    gff_multi = rduckhts_gff_multi,
    gtf_multi = rduckhts_gtf_multi
  )

  for (name in names(readers)) {
    dbExecute(
      con,
      "CREATE OR REPLACE TABLE preserved AS SELECT 1 AS marker UNION ALL SELECT 2"
    )
    expected <- dbGetQuery(con, "SELECT marker FROM preserved ORDER BY marker")
    args <- list(con = con, table_name = "preserved", overwrite = TRUE)
    if (name == "pileup") {
      args$region <- "chr1:1-10"
    }
    if (grepl("_multi$", name)) {
      args$files <- missing_path
    } else {
      args$path <- missing_path
    }
    expect_error(suppressWarnings(do.call(readers[[name]], args)), info = name)
    expect_true(dbExistsTable(con, "preserved"), info = name)
    if (dbExistsTable(con, "preserved")) {
      expect_identical(
        dbGetQuery(con, "SELECT marker FROM preserved ORDER BY marker"),
        expected,
        info = name
      )
    }
  }

  expect_error(
    rduckhts_bed(con, "preserved", missing_path),
    "Table 'preserved' already exists"
  )

  missing_glob <- paste0(tempfile("duckhts_missing_"), "*")
  multi_readers <- list(
    bam = rduckhts_bam_multi,
    bcf = rduckhts_bcf_multi,
    fastq = rduckhts_fastq_multi,
    fasta = rduckhts_fasta_multi,
    bed = rduckhts_bed_multi,
    tabix = rduckhts_tabix_multi,
    gff = rduckhts_gff_multi,
    gtf = rduckhts_gtf_multi
  )
  for (name in names(multi_readers)) {
    warnings <- character()
    result <- withCallingHandlers(
      tryCatch(
        do.call(multi_readers[[name]], list(
          con = con, table_name = "preserved", files = missing_glob
        )),
        error = identity
      ),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    )
    expect_inherits(result, "error", info = name)
    if (inherits(result, "error")) {
      expect_identical(
        conditionMessage(result),
        "Table 'preserved' already exists. Use overwrite = TRUE to replace it.",
        info = name
      )
    }
    expect_identical(warnings, character(), info = name)
  }
})()
