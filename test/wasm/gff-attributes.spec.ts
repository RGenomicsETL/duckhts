import { readFile } from "node:fs/promises";
import path from "node:path";
import { test, expect } from "@playwright/test";

const root = path.resolve(__dirname, "../..");

test("GFF named attributes match raw MAP values in a browser File", async ({ page }) => {
  const fixture = await readFile(path.join(root, "test/data/gff_named_attributes.gff3"), "utf8");
  await page.goto("/scripts/duckdb-wasm-local-test.html");
  const result = await page.evaluate(async (fixture) => {
    const duckdb = await import("/duckdb-browser.mjs");
    const worker = new Worker(new URL("/duckdb-browser-eh.worker.js", location.href));
    const db = new duckdb.AsyncDuckDB(new duckdb.VoidLogger(), worker);
    const url = URL.createObjectURL(new File([fixture], "attributes.gff3"));
    let conn;
    try {
      await db.instantiate(new URL("/duckdb-eh.wasm", location.href).href);
      await db.open({ allowUnsignedExtensions: true });
      conn = await db.connect();
      const extension = new URL("/duckdb-wasm/duckhts.duckdb_extension.wasm", location.href).href;
      await db.registerFileURL(extension, extension, duckdb.DuckDBDataProtocol.HTTP, false);
      await conn.query(`LOAD '${extension}'`);
      const rows = await conn.query(`SELECT ID, Parent, encoded, missing,
        ID IS NOT DISTINCT FROM attributes_map['ID'] AND
        Parent IS NOT DISTINCT FROM attributes_map['Parent'] AND
        encoded IS NOT DISTINCT FROM attributes_map['encoded'] AND
        missing IS NOT DISTINCT FROM attributes_map['missing'] AS same
        FROM read_gff('${url}', attributes := ['Parent','ID','encoded','missing'],
          attributes_map := true, scan_mode := 'sequential') ORDER BY seqname, start`);
      return rows.toArray().map((row) => row.toJSON());
    } finally {
      if (conn) await conn.close();
      await db.terminate();
      worker.terminate();
      URL.revokeObjectURL(url);
    }
  }, fixture);
  expect(result.map((row) => row.ID)).toEqual(["first", "third", null, "fourth", "last"]);
  expect(result.map((row) => row.encoded)).toEqual(["a%3Bb%20c", "%25", null, "a=b", null]);
  expect(result.every((row) => row.same && row.missing === null)).toBe(true);
});
