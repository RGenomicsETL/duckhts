export type DuckhtsPlatform = "wasm_mvp" | "wasm_eh" | "wasm_threads";

export declare const SIGNED: boolean;

/** Create a revocable object URL for a blob-capable build; unavailable on the signed channel. */
export declare function localFileUrl(file: Blob): { url: string; revoke(): void };

export declare const EXTENSION_FILE: "duckhts.duckdb_extension.wasm";

export declare const PLATFORMS: readonly DuckhtsPlatform[];

/** URL of the DuckHTS binary for `platform` under `baseUrl`. */
export declare function extensionUrl(baseUrl: string | URL, platform: string): string;

/** The subset of an `AsyncDuckDBConnection` the loader uses. */
export interface QueryableConnection {
  query(sql: string): Promise<{ toArray(): Array<Record<string, unknown>> }>;
}

/**
 * Run `PRAGMA platform` on `conn` and `LOAD` the matching DuckHTS binary from
 * `baseUrl`, the directory where the application serves this package's `dist/`.
 * The dev channel requires `allow_unsigned_extensions=true` on the connection
 * before loading; this function never changes that setting.
 */
export declare function loadDuckhts(
  conn: QueryableConnection,
  options: { baseUrl: string | URL },
): Promise<{ platform: DuckhtsPlatform }>;
