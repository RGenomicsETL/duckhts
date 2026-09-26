# Setup HTSlib Environment

Sets the \`HTS_PATH\` environment variable to point to the bundled
htslib plugins directory. This enables remote file access via libcurl
plugins (e.g., s3://, gs://, http://) when plugins are available.

## Usage

``` r
setup_hts_env(plugins_dir = NULL)
```

## Arguments

- plugins_dir:

  Optional path to the htslib plugins directory. When NULL, uses the
  bundled plugins directory if available.

## Value

Invisibly returns the previous value of \`HTS_PATH\` (or \`NA\` if
unset).

## Details

Call this before the process opens its first HTS file. htslib discovers
dynamic plugins on first file access and does not rescan \`HTS_PATH\`
later.

## Examples

``` r
if (FALSE) { # \dontrun{
setup_hts_env()

plugins_path <- tempfile("hts_plugins_")
dir.create(plugins_path)
setup_hts_env(plugins_dir = plugins_path)
unlink(plugins_path, recursive = TRUE)
} # }
```
