// Load the DuckHTS extension into a duckdb-wasm connection from a URL the
// application controls, so no request leaves the application's own origin.
//
// `baseUrl` is the directory holding <platform>/duckhts.duckdb_extension.wasm,
// i.e. wherever the application serves this package's dist/ directory.

import { BLOB_CAPABLE, SIGNED } from "./channel.js";
import { localFileUrl as createLocalFileUrl } from "./local-file.js";

export { SIGNED };

export function localFileUrl(file) {
  if (!BLOB_CAPABLE) throw new Error("duckhts: localFileUrl requires a blob-capable build; this signed build does not support blob: URLs");
  return createLocalFileUrl(file);
}

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
  if (!SIGNED) {
    const setting = await conn.query("SELECT current_setting('allow_unsigned_extensions') AS enabled");
    if (setting.toArray()[0]?.enabled !== true) {
      throw new Error("duckhts: dev channel requires allow_unsigned_extensions=true before loadDuckhts; configure DuckDB to allow unsigned extensions explicitly");
    }
  }
  const result = await conn.query("PRAGMA platform");
  const platform = String(result.toArray()[0].platform);
  await conn.query(`LOAD ${sqlString(extensionUrl(baseUrl, platform))}`);
  return { platform };
}
