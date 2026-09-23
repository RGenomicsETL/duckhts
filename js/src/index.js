// Load the DuckHTS extension into a duckdb-wasm connection from a URL the
// application controls, so no request leaves the application's own origin.
//
// `baseUrl` is the directory holding <platform>/duckhts.duckdb_extension.wasm,
// i.e. wherever the application serves this package's dist/ directory.

export const EXTENSION_FILE = "duckhts.duckdb_extension.wasm";

export const PLATFORMS = Object.freeze(["wasm_mvp", "wasm_eh", "wasm_threads"]);

export function extensionUrl(baseUrl, platform) {
  if (!PLATFORMS.includes(platform)) {
    throw new Error(
      `duckhts: no wasm build for DuckDB platform '${platform}'; ` +
        `this package ships ${PLATFORMS.join(", ")}`,
    );
  }
  const base = new URL(baseUrl, globalThis.location?.href);
  if (!base.pathname.endsWith("/")) base.pathname += "/";
  return new URL(`${platform}/${EXTENSION_FILE}`, base).href;
}

function sqlString(value) {
  return `'${value.replaceAll("'", "''")}'`;
}

export async function loadDuckhts(conn, { baseUrl }) {
  if (baseUrl === undefined) {
    throw new Error("duckhts: loadDuckhts(conn, { baseUrl }) requires baseUrl");
  }
  const result = await conn.query("PRAGMA platform");
  const platform = String(result.toArray()[0].platform);
  await conn.query(`LOAD ${sqlString(extensionUrl(baseUrl, platform))}`);
  return { platform };
}
