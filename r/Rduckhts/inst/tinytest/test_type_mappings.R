# Test file for effective DuckDB to R type mappings
# This creates actual DuckDB tables with different types and checks R types returned

library(tinytest)
library(DBI)
library(duckdb)

test_effective_type_mappings <- function() {
  message("Testing effective DuckDB to R type mappings...")

  # Create a test connection
  drv <- duckdb::duckdb()
  con <- dbConnect(drv)

  # Test basic scalar types
  DBI::dbExecute(
    con,
    "CREATE TABLE test_basic (
    bigint_col BIGINT,
    double_col DOUBLE,
    varchar_col VARCHAR,
    boolean_col BOOLEAN,
    date_col DATE,
    timestamp_col TIMESTAMP,
    time_col TIME,
    blob_col BLOB,
    integer_col INTEGER,
    float_col FLOAT,
    smallint_col SMALLINT,
    usmallint_col USMALLINT,
    ubigint_col UBIGINT
  )"
  )

  # Insert test data
  DBI::dbExecute(
    con,
    "INSERT INTO test_basic VALUES (
    42, 3.14, 'test', TRUE, '2023-01-01', '2023-01-01 12:00:00', 
    '12:00:00', 'blob_data', 100, 2.71, 32767, 255, 18446744073709551615
  )"
  )

  # Test ARRAY types
  DBI::dbExecute(
    con,
    "CREATE TABLE test_arrays (
    varchar_array VARCHAR[],
    bigint_array BIGINT[],
    integer_array INTEGER[],
    double_array DOUBLE[],
    float_array FLOAT[],
    boolean_array BOOLEAN[]
  )"
  )

  # Insert array data
  DBI::dbExecute(
    con,
    "INSERT INTO test_arrays VALUES (
    ['a', 'b', 'c'], [1, 2, 3], [10, 20, 30], [1.1, 2.2, 3.3], 
    [0.1, 0.2, 0.3], [true, false, true]
  )"
  )

  # Test MAP types
  DBI::dbExecute(
    con,
    "CREATE TABLE test_maps (
    varchar_map MAP(VARCHAR, VARCHAR),
    bigint_map MAP(VARCHAR, BIGINT),
    double_map MAP(VARCHAR, DOUBLE)
  )"
  )

  # Insert MAP data
  DBI::dbExecute(
    con,
    "INSERT INTO test_maps VALUES (
    map(['key1', 'key2'], ['val1', 'val2']),
    map(['num1', 'num2'], [100, 200]),
    map(['dec1', 'dec2'], [1.1, 2.2])
  )"
  )

  # Test STRUCT type if available
  struct_available <- FALSE
  tryCatch(
    {
      DBI::dbExecute(
        con,
        "CREATE TABLE test_structs (
      struct_col STRUCT(name VARCHAR, age INTEGER, active BOOLEAN)
    )"
      )

      DBI::dbExecute(
        con,
        "INSERT INTO test_structs VALUES (
      struct('test', 25, true)
    )"
      )
      struct_available <- TRUE
      message("STRUCT type available in this DuckDB version")
    },
    error = function(e) {
      message("STRUCT type not available in this DuckDB version")
    }
  )

  # Test UNION type if available
  union_available <- FALSE
  tryCatch(
    {
      DBI::dbExecute(
        con,
        "CREATE TABLE test_unions (
      union_col UNION(VARCHAR, INTEGER, DOUBLE)
    )"
      )

      DBI::dbExecute(
        con,
        "INSERT INTO test_unions VALUES ('test'), (42), (3.14)"
      )
      union_available <- TRUE
      message("UNION type available in this DuckDB version")
    },
    error = function(e) {
      message("UNION type not available in this DuckDB version")
    }
  )

  # Function to get R type
  get_r_type <- function(data, col_name) {
    if (length(data[[col_name]]) == 0) {
      return("empty")
    }
    value <- data[[col_name]][1]
    if (is.list(value)) {
      # For complex types, check nested structure
      if (length(value) > 0 && is.list(value[[1]])) {
        return("data.frame")
      }
      return("list")
    }
    return(typeof(value))
  }

  # Test basic types
  message("Testing basic types...")
  basic_data <- DBI::dbGetQuery(con, "SELECT * FROM test_basic LIMIT 1")
  basic_results <- data.frame(
    duckdb_type = c(
      "BIGINT",
      "DOUBLE",
      "VARCHAR",
      "BOOLEAN",
      "DATE",
      "TIMESTAMP",
      "TIME",
      "BLOB",
      "INTEGER",
      "FLOAT",
      "SMALLINT",
      "USMALLINT",
      "UBIGINT"
    ),
    column_name = c(
      "bigint_col",
      "double_col",
      "varchar_col",
      "boolean_col",
      "date_col",
      "timestamp_col",
      "time_col",
      "blob_col",
      "integer_col",
      "float_col",
      "smallint_col",
      "usmallint_col",
      "ubigint_col"
    ),
    r_type = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(basic_results)) {
    col_name <- basic_results$column_name[i]
    basic_results$r_type[i] <- get_r_type(basic_data, col_name)
  }

  print(basic_results)

  # Test ARRAY types
  message("Testing array types...")
  array_data <- DBI::dbGetQuery(con, "SELECT * FROM test_arrays LIMIT 1")
  array_results <- data.frame(
    duckdb_type = c(
      "VARCHAR[]",
      "BIGINT[]",
      "INTEGER[]",
      "DOUBLE[]",
      "FLOAT[]",
      "BOOLEAN[]"
    ),
    column_name = c(
      "varchar_array",
      "bigint_array",
      "integer_array",
      "double_array",
      "float_array",
      "boolean_array"
    ),
    r_type = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(array_results)) {
    col_name <- array_results$column_name[i]
    array_results$r_type[i] <- get_r_type(array_data, col_name)
    # Also check what typeof returns on array elements
    if (
      is.list(array_data[[col_name]][[1]]) &&
        length(array_data[[col_name]][[1]]) > 0
    ) {
      element_type <- typeof(array_data[[col_name]][[1]][[1]])
      array_results$r_type[i] <- paste0(
        array_results$r_type[i],
        " (elements: ",
        element_type,
        ")"
      )
    }
  }

  print(array_results)

  # Test MAP types
  message("Testing map types...")
  map_data <- DBI::dbGetQuery(con, "SELECT * FROM test_maps LIMIT 1")
  map_results <- data.frame(
    duckdb_type = c(
      "MAP(VARCHAR, VARCHAR)",
      "MAP(VARCHAR, BIGINT)",
      "MAP(VARCHAR, DOUBLE)"
    ),
    column_name = c("varchar_map", "bigint_map", "double_map"),
    r_type = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(map_results)) {
    col_name <- map_results$column_name[i]
    map_results$r_type[i] <- get_r_type(map_data, col_name)

    # Check MAP structure - DuckDB MAPs come as flat lists with keys and values alternating
    if (
      is.list(map_data[[col_name]][[1]]) &&
        length(map_data[[col_name]][[1]]) > 0
    ) {
      map_list <- map_data[[col_name]][[1]]
      message(
        "MAP structure for ",
        col_name,
        ": ",
        length(map_list),
        " elements"
      )
      # MAP should be: key1, value1, key2, value2, ...
      if (length(map_list) >= 2) {
        key_type <- typeof(map_list[[1]])
        value_type <- typeof(map_list[[2]])
        map_results$r_type[i] <- paste0(
          map_results$r_type[i],
          " (keys: ",
          key_type,
          ", values: ",
          value_type,
          ")"
        )
      }
    }
  }

  print(map_results)

  # Test STRUCT types if available
  if (struct_available) {
    message("Testing struct types...")
    struct_data <- DBI::dbGetQuery(con, "SELECT * FROM test_structs LIMIT 1")
    struct_r_type <- get_r_type(struct_data, "struct_col")

    # Check STRUCT structure
    if (is.list(struct_data$struct_col[[1]])) {
      struct_elements <- struct_data$struct_col[[1]]
      message(
        "STRUCT elements: ",
        paste(names(struct_elements), collapse = ", ")
      )
      for (i in seq_along(struct_elements)) {
        message(
          "  ",
          names(struct_elements)[i],
          ": ",
          typeof(struct_elements[[i]])
        )
      }
    }

    struct_results <- data.frame(
      duckdb_type = "STRUCT(name VARCHAR, age INTEGER, active BOOLEAN)",
      column_name = "struct_col",
      r_type = struct_r_type,
      stringsAsFactors = FALSE
    )
    print(struct_results)
  }

  # Test UNION types if available
  if (union_available) {
    message("Testing union types...")
    union_data <- DBI::dbGetQuery(con, "SELECT * FROM test_unions")
    union_r_types <- sapply(1:nrow(union_data), function(i) {
      get_r_type(union_data[i, ], "union_col")
    })

    message("UNION types by row: ", paste(union_r_types, collapse = ", "))

    union_results <- data.frame(
      duckdb_type = "UNION(VARCHAR, INTEGER, DOUBLE)",
      column_name = "union_col",
      r_type = paste(unique(union_r_types), collapse = ", "),
      stringsAsFactors = FALSE
    )
    print(union_results)
  }

  # Compare with our theoretical mappings
  message("Comparing with theoretical mappings...")
  if (requireNamespace("Rduckhts", quietly = TRUE)) {
    mappings <- Rduckhts::duckdb_type_mappings()

    # Test basic types
    message("Basic type mappings:")
    for (i in 1:nrow(basic_results)) {
      db_type <- basic_results$duckdb_type[i]
      actual_r <- basic_results$r_type[i]
      expected_r <- mappings$duckdb_to_r[db_type]

      if (is.na(expected_r)) {
        message(db_type, " -> No mapping in our function, actual: ", actual_r)
      } else if (actual_r != expected_r) {
        message(
          "MISMATCH: ",
          db_type,
          " -> Expected: ",
          expected_r,
          ", Actual: ",
          actual_r
        )
      } else {
        message("MATCH: ", db_type, " -> ", actual_r)
      }
    }

    # Test array types
    message("Array type mappings:")
    for (i in 1:nrow(array_results)) {
      db_type <- array_results$duckdb_type[i]
      actual_r <- gsub(" \\(elements:.*\\)", "", array_results$r_type[i]) # Remove element type info
      expected_r <- mappings$duckdb_to_r[db_type]

      if (is.na(expected_r)) {
        message(db_type, " -> No mapping in our function, actual: ", actual_r)
      } else if (actual_r != expected_r) {
        message(
          "MISMATCH: ",
          db_type,
          " -> Expected: ",
          expected_r,
          ", Actual: ",
          actual_r
        )
      } else {
        message("MATCH: ", db_type, " -> ", actual_r)
      }
    }

    # Test MAP types
    message("MAP type mappings:")
    for (i in 1:nrow(map_results)) {
      db_type <- "MAP" # Our mapping uses generic "MAP"
      actual_r <- gsub(" \\(keys:.*\\)", "", map_results$r_type[i]) # Remove type info
      expected_r <- mappings$duckdb_to_r[db_type]

      if (is.na(expected_r)) {
        message("MAP type -> No mapping in our function, actual: ", actual_r)
      } else if (actual_r != expected_r) {
        message(
          "MISMATCH: MAP -> Expected: ",
          expected_r,
          ", Actual: ",
          actual_r
        )
      } else {
        message("MATCH: MAP -> ", actual_r)
      }
    }
  }

  # Return results for further testing
  dbDisconnect(con, shutdown = TRUE)

  list(
    basic = basic_results,
    arrays = array_results,
    maps = map_results,
    struct = if (struct_available) {
      struct_data <- DBI::dbGetQuery(con, "SELECT * FROM test_structs LIMIT 1")
      data.frame(
        duckdb_type = "STRUCT(name VARCHAR, age INTEGER, active BOOLEAN)",
        column_name = "struct_col",
        r_type = get_r_type(struct_data, "struct_col"),
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    },
    union = if (union_available) {
      union_data <- DBI::dbGetQuery(con, "SELECT * FROM test_unions")
      union_r_types <- sapply(1:nrow(union_data), function(i) {
        get_r_type(union_data[i, ], "union_col")
      })
      data.frame(
        duckdb_type = "UNION(VARCHAR, INTEGER, DOUBLE)",
        column_name = "union_col",
        r_type = paste(unique(union_r_types), collapse = ", "),
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  )
}

# Add proper assertions for the mappings
test_type_mapping_function <- function() {
  message("Testing type mapping function assertions...")

  # Test basic mappings
  mappings <- Rduckhts::duckdb_type_mappings()

  # Check that essential mappings exist
  expect_true("BIGINT" %in% names(mappings$duckdb_to_r))
  expect_true("DOUBLE" %in% names(mappings$duckdb_to_r))
  expect_true("VARCHAR" %in% names(mappings$duckdb_to_r))
  expect_true("BOOLEAN" %in% names(mappings$duckdb_to_r))

  # Check basic mappings are correct
  expect_equal(as.character(mappings$duckdb_to_r["BIGINT"]), "double")
  expect_equal(as.character(mappings$duckdb_to_r["DOUBLE"]), "double")
  expect_equal(paste0(mappings$duckdb_to_r["VARCHAR"]), "character")
  expect_equal(as.character(mappings$duckdb_to_r["BOOLEAN"]), "logical")

  # Check complex type mappings
  expect_true("MAP" %in% names(mappings$duckdb_to_r))
  expect_equal(as.character(mappings$duckdb_to_r["MAP"]), "data.frame")

  # Check array type mappings
  expect_true("VARCHAR[]" %in% names(mappings$duckdb_to_r))
  expect_equal(paste0(mappings$duckdb_to_r["VARCHAR[]"]), "list")
  expect_equal(paste0(mappings$duckdb_to_r["BIGINT[]"]), "list")

  message("Type mapping function assertions passed!")
}

# Run effective mapping test
results <- test_effective_type_mappings()

# Run function assertions
test_type_mapping_function()

message(
  "Type mapping test completed! Check output above for effective type mappings."
)
