loaded_dlls <- getLoadedDLLs()
expect_true("Rduckhts" %in% names(loaded_dlls))
marker_path <- loaded_dlls[["Rduckhts"]][["path"]]
expect_true(file.exists(marker_path))
expect_match(basename(marker_path), "^Rduckhts\\.(so|dll|dylib)$")
