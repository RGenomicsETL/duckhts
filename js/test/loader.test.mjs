import assert from "node:assert/strict";
import { test } from "node:test";

import { extensionUrl, loadDuckhts } from "../src/index.js";

function recordingConnection(platform) {
  const statements = [];
  return {
    statements,
    async query(sql) {
      statements.push(sql);
      return { toArray: () => (sql === "PRAGMA platform" ? [{ platform }] : []) };
    },
  };
}

test("extensionUrl joins base directory and platform", () => {
  assert.equal(
    extensionUrl("https://example.org/vendor/duckhts", "wasm_eh"),
    "https://example.org/vendor/duckhts/wasm_eh/duckhts.duckdb_extension.wasm",
  );
  assert.equal(
    extensionUrl("https://example.org/vendor/duckhts/", "wasm_mvp"),
    "https://example.org/vendor/duckhts/wasm_mvp/duckhts.duckdb_extension.wasm",
  );
});

test("extensionUrl rejects platforms without a wasm build", () => {
  assert.throws(() => extensionUrl("https://example.org/", "linux_amd64"), /linux_amd64/);
});

test("extensionUrl needs an absolute base outside a browser", () => {
  assert.throws(() => extensionUrl("vendor/duckhts", "wasm_eh"), TypeError);
});

test("loadDuckhts loads the binary for the connection's platform", async () => {
  const conn = recordingConnection("wasm_threads");
  const loaded = await loadDuckhts(conn, { baseUrl: "https://example.org/d" });
  assert.deepEqual(loaded, { platform: "wasm_threads" });
  assert.deepEqual(conn.statements, [
    "PRAGMA platform",
    "LOAD 'https://example.org/d/wasm_threads/duckhts.duckdb_extension.wasm'",
  ]);
});

test("loadDuckhts quotes the URL as a SQL string", async () => {
  const conn = recordingConnection("wasm_eh");
  await loadDuckhts(conn, { baseUrl: "https://example.org/it's/" });
  assert.equal(
    conn.statements[1],
    "LOAD 'https://example.org/it''s/wasm_eh/duckhts.duckdb_extension.wasm'",
  );
});

test("loadDuckhts refuses a native platform before issuing LOAD", async () => {
  const conn = recordingConnection("linux_amd64");
  await assert.rejects(loadDuckhts(conn, { baseUrl: "https://example.org/" }), /linux_amd64/);
  assert.deepEqual(conn.statements, ["PRAGMA platform"]);
});

test("loadDuckhts requires baseUrl", async () => {
  await assert.rejects(loadDuckhts(recordingConnection("wasm_eh"), {}), /baseUrl/);
});
