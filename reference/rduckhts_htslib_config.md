# Get the Installed htslib Linking Contract

Resolves headers, an exact shared or static library, linker flags,
enabled features, and build identity from the installed Rduckhts
package. With validation enabled, the receipt and public headers are
also compared with the htslib version reported by the loaded DuckHTS
extension.

## Usage

``` r
rduckhts_htslib_config(link = NULL, validate = TRUE)
```

## Arguments

- link:

  Either \`"shared"\` or \`"static"\`. When omitted, use the link mode
  selected when this Rduckhts package was configured.

- validate:

  Whether to validate installed files, header identity, and the loaded
  htslib runtime version.

## Value

An object of class \`rduckhts_htslib_config\`. Its \`cppflags\` and
\`ldflags\` elements can be consumed by a downstream package configure
script; the remaining fields form the versioned build receipt.

## Details

The shared contract is currently available on native Unix builds. MinGW
and browser-wasm builds expose the static contract unless a shared
htslib was explicitly built. Static consumers must review
\`static_license_note\`.

## Examples

``` r
if (FALSE) { # \dontrun{
config <- rduckhts_htslib_config()
config$cppflags
config$ldflags
config$features
} # }
```
