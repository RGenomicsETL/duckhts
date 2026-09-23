# duckhts for duckdb-wasm

The [DuckHTS](https://github.com/RGenomicsETL/duckhts) DuckDB extension, packaged for
[`@duckdb/duckdb-wasm`](https://www.npmjs.com/package/@duckdb/duckdb-wasm), so a web
application can serve the extension from its own origin instead of fetching it from the
DuckDB community repository at runtime.

The `.wasm` files are the signed community-repository builds, copied byte for byte and
checked against the sha256 values in `artifacts.json`. Loading them does not require
`allow_unsigned_extensions`.

## Use

```sh
npm install duckhts @duckdb/duckdb-wasm@1.33.1-dev57.0
```

Serve `node_modules/duckhts/dist/` with the rest of your static files (for example by
copying it to `vendor/duckhts/`), then:

```js
import { loadDuckhts } from "duckhts";

const conn = await db.connect();
await loadDuckhts(conn, { baseUrl: "/vendor/duckhts/" });

const peaks = await conn.query(
  `SELECT chrom, start, "end" FROM read_bed('${location.origin}/data/peaks.bed')`,
);
```

`loadDuckhts` runs `PRAGMA platform` and loads
`<baseUrl>/<platform>/duckhts.duckdb_extension.wasm` for `wasm_mvp`, `wasm_eh` or
`wasm_threads`. The SQL functions are documented in the
[function catalog](https://github.com/RGenomicsETL/duckhts/blob/main/r/Rduckhts/inst/function_catalog/functions.md).

## Runtime support

The peer range is `1.33.1-dev57.0 || >=1.33.1`: the one tested prerelease, or any stable
release from 1.33.1. npm compares prerelease tags as text (`dev6` sorts after `dev57`),
so a `>=` floor on a prerelease would admit older builds.

| duckdb-wasm | DuckDB | Status |
|---|---|---|
| `1.33.1-dev57.0` | v1.5.4 | Tested: loads with unsigned extensions disallowed; `read_bed` / `read_gff` over same-origin HTTP |
| `1.32.0` | v1.4.3 | Fails to load: https://github.com/RGenomicsETL/duckhts/issues/247 |

## Browser file access

DuckHTS reads files through htslib using worker-synchronous XHR. Same-origin HTTP URLs
work, and cross-origin URLs need permissive CORS.

### Local files

DuckHTS builds from https://github.com/RGenomicsETL/duckhts/pull/248 onward read
`blob:` object URLs, so a page can pass a dropped or picked `File` straight to a reader.
The URL is a capability, not a path: unguessable, read-only, scoped to the page that
created it, and revocable. Any `Blob` the page holds works the same way: a file from
`<input type="file">` or a drop event, an OPFS file (`handle.getFile()`), a Blob kept in
IndexedDB, or bytes the page fetched with its own credentials. SQL stays the same:
`read_bed(url)` does not care where the bytes came from, just as native builds read
paths, `data:` URLs and `/dev/fd/N`.

**This package does not expose it yet.** The signed binaries pinned in `artifacts.json`
(DuckHTS 1.5.2) predate the handler, so a `blob:` URL fails with them. The package will
export `localFileUrl(file)` → `{ url, revoke }` in the release that pins blob-capable
builds. Until then, a page using a locally built extension (`allowUnsignedExtensions:
true`) can call `URL.createObjectURL(file)` itself:

```js
const url = URL.createObjectURL(fileInput.files[0]);
try {
  const rows = await conn.query(`SELECT chrom, start, "end" FROM read_bed('${url}')`);
} finally {
  URL.revokeObjectURL(url);
}
```

Keep the URL live until every query or result stream using it finishes, including
queries against views that retain it. Object URLs have no sibling filenames: automatic
index discovery falls back to streaming, and a region query needs a second URL for the
index passed as `index_path := '<index blob URL>'`. If a transport ignores Range, the
backend caches the full response in JavaScript, which costs browser memory
proportional to the file size.

Files registered with duckdb-wasm (`registerFileBuffer`, `registerFileHandle`,
`registerFileText`) remain **invisible** to DuckHTS readers. The file-system integration
is a separate option in https://github.com/RGenomicsETL/duckhts/issues/246.

## Development

```sh
npm ci
npm test               # loader, local-file helper and staging tests, offline
npm run test:browser   # stages the binaries, then runs Playwright against both runtimes
```

Render this README with `Rscript -e 'knitr::knit("js/README.Rmd", "js/README.md")'`
from the repository root. Local-extension blob tests run through `make wasm-playwright-test`.

`npm run stage` is the only step that uses the network. `npm pack` runs it through
`prepack`.

To move to a new DuckHTS release, update `artifacts.json` (the version, the community
DuckDB path and the three sha256 values) and the `version` in `package.json`.

## Licence

GPL-2.0-or-later: the GNU General Public License, version 2 or (at your option) any
later version. See [`LICENSE`](LICENSE). The wasm binaries also contain HTSlib,
htscodecs, libBigWig, cgranges, VariantKey, zlib, bzip2 and liblzma under their own
licences, reproduced in [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).
