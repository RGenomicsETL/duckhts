// Object URLs for browser Files and Blobs, for DuckHTS readers built with the
// Emscripten blob: handler (https://github.com/RGenomicsETL/duckhts/pull/248).
//
// Not exported from the package entry yet: the signed binaries pinned in
// artifacts.json predate the handler, so this would fail with the package's own
// dist/. Export it from index.js in the release that pins blob-capable builds;
// test/exports.test.mjs enforces this.
//
// The caller owns the URL. Keep it live until every query or stream using it
// has finished, including queries against views that retain the URL.
export function localFileUrl(file) {
  const url = URL.createObjectURL(file);
  return { url, revoke: () => URL.revokeObjectURL(url) };
}
