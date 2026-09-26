// Object URLs for browser Files and Blobs, for DuckHTS readers built with the
// Emscripten blob: handler (https://github.com/RGenomicsETL/duckhts/pull/248).
//
// The caller owns the URL. Keep it live until every query or stream using it
// has finished, including queries against views that retain the URL.
export function localFileUrl(file) {
  const url = URL.createObjectURL(file);
  return { url, revoke: () => URL.revokeObjectURL(url) };
}
