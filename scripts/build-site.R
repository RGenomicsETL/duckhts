#!/usr/bin/env Rscript

if (!requireNamespace("pkgdown", quietly = TRUE) ||
    !requireNamespace("litedown", quietly = TRUE)) {
  stop("Building the site requires pkgdown and litedown", call. = FALSE)
}

site_url <- "https://rgenomicsetl.github.io/duckhts/"
package_dir <- "r/Rduckhts"
package_name <- read.dcf(file.path(package_dir, "DESCRIPTION"), fields = "Package")[[1L]]

unlink("_site", recursive = TRUE)
dir.create("_site")
file.create(file.path("_site", ".nojekyll"))
site_root <- normalizePath("_site", mustWork = TRUE)

# pkgdown requires an empty destination (or an existing pkgdown site).
destination <- file.path(site_root, package_name)
dir.create(destination)
message("== pkgdown ", package_name, " ==")
pkgdown::build_site(
  pkg = package_dir,
  new_process = FALSE,
  install = FALSE,
  preview = FALSE,
  override = list(destination = destination)
)

# Keep legacy pkgdown URLs working while the project README owns /index.html.
package_pages <- list.files(destination, pattern = "\\.html$", recursive = TRUE)
legacy_pages <- setdiff(package_pages, "index.html")
for (page in legacy_pages) {
  old_path <- file.path(site_root, page)
  if (file.exists(old_path)) {
    stop("Legacy pkgdown URL conflicts with project page: ", page, call. = FALSE)
  }
  dir.create(dirname(old_path), recursive = TRUE, showWarnings = FALSE)
  target <- paste0(site_url, package_name, "/", page)
  writeLines(c(
    "<!doctype html>",
    "<html lang=\"en\"><head><meta charset=\"utf-8\">",
    sprintf("<meta http-equiv=\"refresh\" content=\"0; url=%s\">", target),
    sprintf("<link rel=\"canonical\" href=\"%s\">", target),
    "</head><body>",
    sprintf("<a href=\"%s\">Rduckhts documentation</a>", target),
    "</body></html>"
  ), old_path)
}

landing_css <- normalizePath("scripts/site/landing.css", winslash = "/", mustWork = TRUE)
landing_header <- normalizePath("scripts/site/landing-header.html", winslash = "/", mustWork = TRUE)
docs_header <- normalizePath("scripts/site/docs-header.html", winslash = "/", mustWork = TRUE)
source_docs <- sort(list.files("docs", pattern = "\\.md$", full.names = TRUE))
if (length(source_docs) > 0L) {
  header <- readLines(landing_header, warn = FALSE)
  header <- sub(
    '    <a href="https://github',
    '    <a href="docs/">Documentation</a>\n    <a href="https://github',
    header,
    fixed = TRUE
  )
  landing_header <- tempfile(fileext = ".html")
  writeLines(header, landing_header)
}

metadata <- function(title, header) {
  c(
    "---",
    paste0("title: ", title),
    "output:",
    "  html:",
    "    options:",
    "      toc: true",
    "    meta:",
    paste0(
      "      css: [\"@default@1.14.69\", \"@article@1.14.69\", ",
      "\"@site@1.14.69\", \"", landing_css, "\"]"
    ),
    paste0("      include_before: \"", header, "\""),
    "---"
  )
}

readme <- readLines("README.md", warn = FALSE, encoding = "UTF-8")
title_heading <- which(readme == "# DuckHTS")
if (length(title_heading) > 0L) {
  readme <- readme[-title_heading[[1L]]]
}
landing_path <- file.path(site_root, "index.html")
litedown::mark(
  text = c(metadata("DuckHTS", landing_header), readme),
  output = landing_path
)
if (length(source_docs) > 0L) {
  unlink(landing_header)
}

# Relative README source links work on GitHub, but not at the Pages root.
# Pin them to the source revision that produced this deployment.
revision <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
if (length(revision) != 1L || !grepl("^[0-9a-f]{40}$", revision)) {
  stop("Cannot identify the source revision for README links", call. = FALSE)
}
landing <- xml2::read_html(landing_path)
links <- xml2::xml_find_all(landing, "//a[@href]")
for (link in links) {
  href <- xml2::xml_attr(link, "href")
  if (grepl("^(#|/|\\./|Rduckhts/|docs/|[a-z]+:)", href)) {
    next
  }
  path <- sub("#.*$", "", href)
  if (!file.exists(path)) {
    stop("README link has no source file: ", href, call. = FALSE)
  }
  xml2::xml_set_attr(
    link, "href",
    paste0("https://github.com/RGenomicsETL/duckhts/blob/", revision, "/", href)
  )
}
xml2::write_html(landing, landing_path)

if (length(source_docs) > 0L) {
  docs_destination <- file.path(site_root, "docs")
  dir.create(docs_destination)
  index <- c("# Documentation", "", vapply(source_docs, function(source) {
    sprintf("- [%s](%s.html)", tools::file_path_sans_ext(basename(source)),
            tools::file_path_sans_ext(basename(source)))
  }, character(1L)))
  litedown::mark(
    text = c(metadata("Documentation", docs_header), index),
    output = file.path(docs_destination, "index.html")
  )
  for (source in source_docs) {
    name <- tools::file_path_sans_ext(basename(source))
    markdown <- readLines(source, warn = FALSE, encoding = "UTF-8")
    litedown::mark(
      text = c(metadata(name, docs_header), markdown),
      output = file.path(docs_destination, paste0(name, ".html"))
    )
    file.copy(source, file.path(docs_destination, basename(source)), overwrite = TRUE)
  }
}
message("Legacy pkgdown redirects: ", length(legacy_pages))
