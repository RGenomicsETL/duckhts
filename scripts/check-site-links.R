#!/usr/bin/env Rscript

site_root <- normalizePath("_site", mustWork = TRUE)
site_url <- "https://rgenomicsetl.github.io/duckhts/"
pages <- list.files(site_root, pattern = "\\.html$", recursive = TRUE)
if (!all(c("index.html", "Rduckhts/index.html") %in% pages)) {
  stop("The project and package index pages are required", call. = FALSE)
}

checked <- 0L
for (page in pages) {
  document <- xml2::read_html(file.path(site_root, page))
  links <- if (page == "index.html") {
    xml2::xml_find_all(document, "//a[@href] | //img[@src]")
  } else {
    xml2::xml_find_all(document, "//nav//a[@href]")
  }
  for (link in links) {
    href <- xml2::xml_attr(link, if (xml2::xml_name(link) == "img") "src" else "href")
    if (startsWith(href, site_url)) {
      href <- substring(href, nchar(site_url) + 1L)
      path <- sub("[?#].*$", "", href)
      target <- file.path(site_root, path)
    } else if (grepl("^(//|[a-z]+:)", href)) {
      next
    } else {
      path <- sub("[?#].*$", "", href)
      target <- file.path(site_root, dirname(page), path)
    }
    target <- utils::URLdecode(target)
    if (dir.exists(target)) {
      target <- file.path(target, "index.html")
    }
    if (!file.exists(target)) {
      stop(page, ": broken link ", href, call. = FALSE)
    }
    fragment <- if (grepl("#", href, fixed = TRUE)) {
      utils::URLdecode(sub("^.*#", "", href))
    } else {
      ""
    }
    if (nzchar(fragment) && grepl("\\.html$", target)) {
      ids <- xml2::xml_attr(xml2::xml_find_all(xml2::read_html(target), "//*[@id]"), "id")
      if (!fragment %in% ids) {
        stop(page, ": missing anchor ", href, call. = FALSE)
      }
    }
    checked <- checked + 1L
  }
}

redirects <- setdiff(pages[!startsWith(pages, "Rduckhts/") &
                           !startsWith(pages, "docs/")], "index.html")
package_pages <- list.files(file.path(site_root, "Rduckhts"),
                            pattern = "\\.html$", recursive = TRUE)
if (!identical(sort(redirects), sort(setdiff(package_pages, "index.html")))) {
  stop("The legacy redirects do not cover every pkgdown page", call. = FALSE)
}
for (page in redirects) {
  document <- xml2::read_html(file.path(site_root, page))
  expected <- paste0(site_url, "Rduckhts/", page)
  canonical <- xml2::xml_attr(
    xml2::xml_find_first(document, "//link[@rel='canonical']"), "href"
  )
  refresh <- xml2::xml_attr(
    xml2::xml_find_first(document, "//meta[@http-equiv='refresh']"), "content"
  )
  if (!identical(canonical, expected) || !identical(refresh, paste0("0; url=", expected)) ||
      !file.exists(file.path(site_root, "Rduckhts", page))) {
    stop("Invalid legacy redirect: ", page, call. = FALSE)
  }
}
cat("Checked", checked, "local index/nav links and", length(redirects),
    "legacy redirects across", length(pages), "HTML pages\n")
