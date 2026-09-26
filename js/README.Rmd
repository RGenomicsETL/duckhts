# duckhts for duckdb-wasm

The [DuckHTS](https://github.com/RGenomicsETL/duckhts) DuckDB extension, packaged for
[`@duckdb/duckdb-wasm`](https://www.npmjs.com/package/@duckdb/duckdb-wasm), so a web
application can serve the extension from its own origin instead of fetching it from the
DuckDB community repository at runtime.

Each npm release contains one pinned set of `.wasm` files, verified byte for byte against
its channel's manifest before packaging. `SIGNED` reports which channel supplied them.

## Channels

| npm dist-tag | Install | Binary source | DuckHTS version |
|---|---|---|---|
| `latest` | `npm install duckhts@latest` | signed community repository (`artifacts.json`) | released version |
| `dev` | `npm install duckhts@dev` | unsigned GitHub Actions run (`artifacts-dev.json`) | `1.5.2.9002` as npm `1.5.2-9002` |

The `dev` loader refuses to load unless DuckDB has `allow_unsigned_extensions=true`.
For duckdb-wasm, opt in explicitly when opening the database:

```js
await db.open({ allowUnsignedExtensions: true });
```

The loader checks `current_setting('allow_unsigned_extensions')`; it never enables it.
Unsigned extensions are not authenticated by DuckDB's community signature. Only opt in
when you trust the pinned commit and artifact checksums; enabling the setting applies
to the database, not just DuckHTS. The staging script verifies each download's sha256
against `artifacts-dev.json` before packaging. GitHub Actions artifacts expire; once
published to npm, the package tarball retains the verified binaries permanently.

## Use

```sh
npm install duckhts@latest @duckdb/duckdb-wasm@1.33.1-dev57.0
```

Serve `node_modules/duckhts/dist/` with the rest of your static files (for example by
copying it to `vendor/duckhts/`), then:

```js
import { loadDuckhts, SIGNED } from "duckhts";

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
| `1.33.1-dev57.0` | v1.5.4 | Both channels: HTTP readers; `dev` also reads `blob:` files after unsigned opt-in |
| `1.32.0` | v1.4.3 | `dev` loads; signed DuckHTS 1.5.2 fails: https://github.com/RGenomicsETL/duckhts/issues/247 |

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

The `dev` binaries include the handler. `localFileUrl(file)` creates an object URL and
returns `{ url, revoke }`. Keep it live while queries use it, then revoke it:

```js
import { localFileUrl, SIGNED } from "duckhts";

if (!SIGNED) {
  const local = localFileUrl(fileInput.files[0]);
  try {
    const rows = await conn.query(`SELECT chrom, start, "end" FROM read_bed('${local.url}')`);
  } finally {
    local.revoke();
  }
}
```

On `latest` (DuckHTS 1.5.2), `localFileUrl` throws because the signed binaries predate
the handler. `SIGNED` lets applications select a supported input path.

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
npm test               # loader, local-file helper and staging tests, offline (signed)
DUCKHTS_NPM_CHANNEL=dev npm test
DUCKHTS_NPM_CHANNEL=signed npm run test:browser
DUCKHTS_NPM_CHANNEL=dev npm run test:browser
DUCKHTS_NPM_CHANNEL=dev npm run check:version
```

Render this README with `Rscript -e 'knitr::knit("js/README.Rmd", "js/README.md")'`
from the repository root. Local-extension blob tests run through `make wasm-playwright-test`.

Set `DUCKHTS_NPM_CHANNEL=signed|dev` for staging, browser tests and packaging.
`npm run stage` is the only network step; `npm pack` and `npm publish` run it through
`prepack`. For offline dev staging from an existing download tree, use
`node scripts/stage.mjs dev js/artifacts-dev.json js/dist /path/to/gh-download`
from the repository root. `npm run check:version` requires a matching version in
`description.yml`, the selected manifest, `package.json` and `package-lock.json`.
It accepts `-- <channel> <pr|publish>`; without a mode it uses strict publish checks.
Pull requests check identity for the checkout's channel and skip the other channel's
identity check while testing both runtimes. The publish workflow dispatch takes
`channel=dev|signed` and `publish=true`; dispatch always requires the selected
channel to match the checkout's version form.

## Licence

GPL-2.0-or-later: the GNU General Public License, version 2 or (at your option) any
later version. See [`LICENSE`](LICENSE). The wasm binaries also contain HTSlib,
htscodecs, libBigWig, cgranges, VariantKey, zlib, bzip2 and liblzma under their own
licences, reproduced in [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).
