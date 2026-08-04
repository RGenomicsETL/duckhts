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

sql_quote_string <- function(x) {
  sprintf("'%s'", gsub("'", "''", x, fixed = TRUE))
}

.validate_nonnegative_integer_param <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      value < 0 || value > .Machine$integer.max || value != floor(value)) {
    stop(
      sprintf(
        "%s must be a single whole number between 0 and %d",
        name, .Machine$integer.max
      ),
      call. = FALSE
    )
  }
  as.integer(value)
}

.validate_scan_mode_param <- function(value, name = "scan_mode") {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(sprintf("%s must be 'auto' or 'sequential'", name), call. = FALSE)
  }
  value <- tolower(value)
  if (!value %in% c("auto", "sequential")) {
    stop(sprintf("%s must be 'auto' or 'sequential'", name), call. = FALSE)
  }
  value
}

sql_varchar_list_literal <- function(x, name = "value") {
  if (is.null(x) || length(x) == 0L) {
    return("[]::VARCHAR[]")
  }
  if (!is.character(x) || anyNA(x) || any(!nzchar(x))) {
    stop(name, " must be a character vector without NA or empty values", call. = FALSE)
  }
  sprintf("[%s]", paste(vapply(x, sql_quote_string, character(1)), collapse = ", "))
}

sql_map_literal <- function(x, name = "column_map", allow_empty = FALSE) {
  if (is.null(x) || length(x) == 0) {
    if (allow_empty) {
      return("map([]::VARCHAR[], []::VARCHAR[])")
    }
    stop(name, " must be non-empty", call. = FALSE)
  }
  nm <- names(x)
  if (is.null(nm) || any(!nzchar(nm))) {
    stop(name, " must be a named character vector", call. = FALSE)
  }
  keys <- paste(vapply(nm, sql_quote_string, character(1)), collapse = ", ")
  vals <- paste(vapply(as.character(x), sql_quote_string, character(1)), collapse = ", ")
  sprintf("map([%s], [%s])", keys, vals)
}

duckdb_path_exists <- function(con, path) {
  child_pattern <- paste0(sub("/+$", "", path), "/**")
  tryCatch(
    isTRUE(DBI::dbGetQuery(
      con,
      sprintf(
        "SELECT EXISTS(SELECT 1 FROM glob([%s, %s]) LIMIT 1) AS exists",
        sql_quote_string(path),
        sql_quote_string(child_pattern)
      )
    )$exists[[1]]),
    error = function(e) NA
  )
}

read_munge_column_map_file <- function(path) {
  tbl <- utils::read.delim(
    path,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  if (ncol(tbl) < 2) {
    stop("column_map_file must be a two-column TSV with source and canonical names", call. = FALSE)
  }
  out <- structure(tbl[[1]], names = toupper(tbl[[2]]))
  out[nzchar(names(out))]
}

resolve_munge_column_map <- function(raw_map, available_columns) {
  if (is.null(raw_map) || length(raw_map) == 0 || is.null(available_columns) || length(available_columns) == 0) {
    return(character())
  }
  available_columns <- as.character(available_columns)
  available_upper <- toupper(available_columns)
  canonical_names <- unique(names(raw_map))
  selected <- character(length(canonical_names))
  names(selected) <- canonical_names
  for (i in seq_along(canonical_names)) {
    canonical <- canonical_names[[i]]
    candidates <- as.character(raw_map[names(raw_map) == canonical])
    for (candidate in candidates) {
      idx <- match(toupper(candidate), available_upper)
      if (!is.na(idx)) {
        selected[[i]] <- available_columns[[idx]]
        break
      }
    }
  }
  selected[nzchar(selected)]
}

read_munge_preset_map <- function(con, preset) {
  preset_sql <- sql_quote_string(preset)
  out <- DBI::dbGetQuery(con, sprintf("SELECT duckdb_munge_preset_map(%s) AS m", preset_sql))
  if (nrow(out) != 1 || is.null(out$m[[1]])) {
    stop("duckdb_munge: unknown preset", call. = FALSE)
  }
  m <- out$m[[1]]
  if (!is.data.frame(m) || !all(c("key", "value") %in% names(m))) {
    stop("duckdb_munge: failed to read preset map", call. = FALSE)
  }
  structure(as.character(m$value), names = toupper(as.character(m$key)))
}
