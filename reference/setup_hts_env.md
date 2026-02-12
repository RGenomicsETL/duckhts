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

Call this before querying remote URLs to allow htslib to locate its
plugins.

## Examples

``` r
# Use bundled plugins if present
setup_hts_env()

# Or set an explicit plugins directory
# plugins_path <- "/path/to/htslib/plugins"
# setup_hts_env(plugins_dir = plugins_path)
```
