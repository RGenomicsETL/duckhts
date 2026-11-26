import { test, expect } from "@playwright/test";
import { bedFixture, fastqRecordCount, vcfFixture } from "./fixture-oracles";

// SET duckhts_file_system = 'duckdb' opens every input through the file system of
// the connection. In duckdb-wasm that is its browser file system, so files
// registered with registerFileBuffer, and the index beside a registered VCF, are
// readable by name. With the default htslib route the same names are not files.
test("duckhts_file_system = 'duckdb' reads registered files and same-origin HTTP", async ({ page }) => {
  const { text: bed, rows: bedRows } = await bedFixture();
  const { bytes: vcf, index, rows: vcfRows, overlapping } = await vcfFixture();
  const fastqCount = await fastqRecordCount("r/Rduckhts/inst/extdata/r1.fq");

  await page.goto("/scripts/duckdb-wasm-local-test.html");
  const result = await page.evaluate(async ({ bed, vcf, index }) => {
    const duckdb = await import("/duckdb-browser.mjs");
    const workerUrl = new URL("/duckdb-browser-base.worker.js", location.href).href;
    const worker = new Worker(workerUrl);
    const db = new duckdb.AsyncDuckDB(new duckdb.VoidLogger(), worker);
    let conn;
    try {
      await db.instantiate(new URL("/duckdb-base.wasm", location.href).href);
      await db.open({ allowUnsignedExtensions: true });
      conn = await db.connect();
      const extensionPath = new URL("/duckdb-wasm/duckhts.duckdb_extension.wasm", location.href).href;
      await db.registerFileURL(extensionPath, extensionPath, duckdb.DuckDBDataProtocol.HTTP, false);
      await conn.query(`LOAD '${extensionPath}'`);
      const rows = async (sql) => JSON.parse(JSON.stringify((await conn.query(sql)).toArray(),
        (_, value) => typeof value === "bigint" ? Number(value) : value));
      const error = async (sql) => {
        try {
          await conn.query(sql);
          return null;
        } catch (e) {
          return String(e);
        }
      };
      await db.registerFileText("route.bed", bed);
      await db.registerFileBuffer("route.vcf.gz", new Uint8Array(vcf));
      await db.registerFileBuffer("route.vcf.gz.tbi", new Uint8Array(index));
      const bedSql = (name) => `SELECT chrom, start, "end", name FROM read_bed('${name}')`;
      const vcfSql = (options) => `SELECT CHROM AS chrom, POS AS pos, ID AS id, REF AS ref,
        array_to_string(ALT, ',') AS alt FROM read_bcf('route.vcf.gz', tidy_format := false ${options})`;
      const httpFastq = new URL("/extdata/r1.fq", location.href).href;
      const fastqSql = `SELECT count(*) AS n FROM read_fastq('${httpFastq}')`;

      const htslibBed = await error(bedSql("route.bed"));
      await conn.query("SET duckhts_file_system = 'duckdb'");
      const setting = await rows("SELECT current_setting('duckhts_file_system') AS value");
      const routedBed = await rows(bedSql("route.bed"));
      const routedVcf = await rows(vcfSql(""));
      // No index_path: the search tries route.vcf.gz.csi first, which duckdb-wasm
      // opens as a one-byte buffer (https://github.com/duckdb/duckdb-wasm/issues/2332);
      // it is not a CSI, so the search goes on to the registered route.vcf.gz.tbi.
      const routedRegion = await rows(vcfSql(", region := 'chr1:12-18'"));
      const routedFastq = await rows(fastqSql);
      const missing = await error(bedSql("not-registered.bed"));
      await conn.query("RESET duckhts_file_system");
      const resetBed = await error(bedSql("route.bed"));
      return { htslibBed, setting, routedBed, routedVcf, routedRegion, routedFastq, missing, resetBed };
    } finally {
      if (conn) await conn.close();
      await db.terminate();
      worker.terminate();
    }
  }, { bed, vcf: [...vcf], index: [...index] });

  expect(result.setting).toEqual([{ value: "duckdb" }]);
  expect(result.routedBed).toEqual(bedRows);
  expect(result.routedVcf).toEqual(vcfRows);
  expect(result.routedRegion).toEqual(overlapping("chr1", 12, 18));
  expect(result.routedFastq).toEqual([{ n: fastqCount }]);
  // Failure controls: the registered name is not a file to htslib, before the
  // setting and after RESET.
  expect(result.htslibBed).toContain("read_bed: failed to open file");
  expect(result.resetBed).toContain("read_bed: failed to open file");
  // An unregistered name opens as one byte of junk and fails on its content, not
  // with "not found" (https://github.com/duckdb/duckdb-wasm/issues/2332). This
  // assertion fails once duckdb-wasm reports the missing file; then expect the name.
  expect(result.missing).toContain("read_bed: BED line has fewer than 3 tab-delimited fields");
});
