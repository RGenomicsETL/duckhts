import { readFile } from "node:fs/promises";
import path from "node:path";
import { gunzipSync } from "node:zlib";
import { test, expect } from "@playwright/test";

const root = path.resolve(__dirname, "../..");

for (const ignoreRange of [false, true]) {
  test(`local File readers (${ignoreRange ? "full-body fallback" : "browser ranges"})`, async ({ page }) => {
    const bed = await readFile(path.join(root, "test/data/fixture_mixed_regions.bed"), "utf8");
    const gff = await readFile(path.join(root, "test/data/gff_attrs.gff3"), "utf8");
    const vcf = await readFile(path.join(root, "test/data/region_union.vcf.gz"));
    const index = await readFile(path.join(root, "test/data/region_union.vcf.gz.tbi"));
    // Independent oracles: parse the physical fixture records, including duplicates.
    const bedRows = bed.trim().split("\n").map((line) => {
      const [chrom, start, end, name] = line.split("\t");
      return { chrom, start: Number(start), end: Number(end), name };
    });
    const gffRows = gff.trim().split("\n").filter((line) => !line.startsWith("#")).map((line) => {
      const [seqname, , feature, start, end] = line.split("\t");
      return { seqname, feature, start: Number(start), end: Number(end) };
    });
    const records = gunzipSync(vcf).toString("utf8").trim().split("\n")
      .filter((line) => !line.startsWith("#")).map((line) => {
        const [chrom, pos, id, ref, alt, , , info] = line.split("\t");
        const end = info.split(";").find((field) => field.startsWith("END="));
        return {
          row: { chrom, pos: Number(pos), id, ref, alt },
          end: end ? Number(end.slice(4)) : Number(pos) + ref.length - 1,
        };
      });
    const expectedRegion = records.filter(({ row, end }) => row.chrom === "chr1" && row.pos <= 18 && end >= 12);

    const browserLogs: string[] = [];
    page.on("console", (message) => browserLogs.push(message.text()));
    await page.route("**/duckhts-package.js", (route) => route.fulfill({
      path: path.join(root, "js/src/local-file.js"), contentType: "text/javascript",
    }));
    await page.goto("/scripts/duckdb-wasm-local-test.html");
    const result = await page.evaluate(async ({ bed, gff, vcf, index, ignoreRange }) => {
      const duckdb = await import("/duckdb-browser.mjs");
      const { localFileUrl } = await import("/duckhts-package.js");
      // Simulate a transport ignoring Range without using DuckHTS as the oracle.
      // The normal mode uses Chromium's XMLHttpRequest unchanged.
      const wrapper = localFileUrl(new Blob([`
        ${ignoreRange ? `
          const setHeader = XMLHttpRequest.prototype.setRequestHeader;
          XMLHttpRequest.prototype.setRequestHeader = function(name, value) {
            if (name.toLowerCase() !== 'range') setHeader.call(this, name, value);
          };
        ` : ""}
        importScripts(${JSON.stringify(new URL("/duckdb-browser-eh.worker.js", location.href).href)});
      `], { type: "text/javascript" }));
      const worker = new Worker(wrapper.url);
      const db = new duckdb.AsyncDuckDB(new duckdb.VoidLogger(), worker);
      const local = {
        bed: localFileUrl(new File([bed], "dropped.bed")),
        gff: localFileUrl(new File([gff], "dropped.gff3")),
        vcf: localFileUrl(new File([new Uint8Array(vcf)], "dropped.vcf.gz")),
        index: localFileUrl(new File([new Uint8Array(index)], "dropped.vcf.gz.tbi")),
        // Exceed hFILE's read buffer to exercise repeated ranges and cached EOF.
        large: localFileUrl(new File([bed.repeat(2048)], "large.bed")),
        // Codex review of 8d0fd3d: an empty File is a valid empty input, as a
        // zero-byte file is natively, not a missing one.
        empty: localFileUrl(new File([], "empty.bed")),
      };
      let conn;
      try {
        await db.instantiate(new URL("/duckdb-eh.wasm", location.href).href);
        await db.open({ allowUnsignedExtensions: true });
        conn = await db.connect();
        const extensionPath = new URL("/duckdb-wasm/duckhts.duckdb_extension.wasm", location.href).href;
        await db.registerFileURL(extensionPath, extensionPath,
          duckdb.DuckDBDataProtocol.HTTP, false);
        await conn.query(`LOAD '${extensionPath}'`);
        console.log("blob-tests: reads begin");
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
        // Revocation reaches the DuckDB worker asynchronously: retry (bounded)
        // until the query fails, then return that error for the assertions.
        const errorAfterRevoke = async (sql) => {
          for (let attempt = 0; attempt < 200; attempt += 1) {
            const message = await error(sql);
            if (message !== null) return message;
            await new Promise((resolve) => setTimeout(resolve, 10));
          }
          return null;
        };
        const bedSql = (url) => `SELECT chrom, start, "end", name FROM read_bed('${url}')`;
        const vcfSql = (options) => `SELECT CHROM AS chrom, POS AS pos, ID AS id, REF AS ref,
          array_to_string(ALT, ',') AS alt FROM read_bcf('${local.vcf.url}', tidy_format := false ${options})`;
        const bedResult = await rows(bedSql(local.bed.url));
        const gffResult = await rows(`SELECT seqname, feature, start, "end" FROM read_gff('${local.gff.url}')`);
        // Auto mode must stream when the synthesized blob sidecar URL is absent.
        const streamed = await rows(vcfSql(""));
        const indexed = await rows(vcfSql(`, index_path := '${local.index.url}', region := 'chr1:12-18'`));
        const large = await rows(bedSql(local.large.url));
        const emptyBed = await rows(bedSql(local.empty.url));
        const emptyGff = await rows(`SELECT seqname FROM read_gff('${local.empty.url}')`);
        console.log("blob-tests: expected errors follow");
        const missingIndex = await error(vcfSql(", region := 'chr1:12-18'"));
        local.index.revoke();
        const revokedIndex = await errorAfterRevoke(vcfSql(`, index_path := '${local.index.url}', region := 'chr1:12-18'`));
        local.bed.revoke();
        const revoked = await errorAfterRevoke(bedSql(local.bed.url));
        const unknown = await error(bedSql(`blob:${location.origin}/00000000-0000-0000-0000-000000000000`));
        // Option B is not provided by a blob scheme handler.
        const bytes = new TextEncoder().encode(bed);
        await db.registerFileBuffer("registered-buffer.bed", bytes);
        await db.registerFileText("registered-text.bed", bed);
        await db.registerFileHandle("registered-handle.bed", new File([bed], "handle.bed"),
          duckdb.DuckDBDataProtocol.BROWSER_FILEREADER, true);
        const registered = [];
        for (const name of ["buffer", "text", "handle"]) {
          registered.push(await error(bedSql(`registered-${name}.bed`)));
        }
        const alive = await rows("SELECT 42 AS answer");
        return { bedResult, gffResult, streamed, indexed, large, emptyBed, emptyGff,
          missingIndex, revokedIndex, revoked, unknown, registered, alive };
      } finally {
        if (conn) await conn.close();
        await db.terminate();
        worker.terminate();
        wrapper.revoke();
        for (const file of Object.values(local)) file.revoke();
      }
    }, { bed, gff, vcf: [...vcf], index: [...index], ignoreRange });

    expect(result.bedResult).toEqual(bedRows);
    expect(result.gffResult).toEqual(gffRows);
    expect(result.streamed).toEqual(records.map(({ row }) => row));
    expect(result.indexed).toEqual(expectedRegion.map(({ row }) => row));
    expect(result.large).toEqual(Array.from({ length: 2048 }, () => bedRows).flat());
    expect(result.emptyBed).toEqual([]);
    expect(result.emptyGff).toEqual([]);
    const readMarker = browserLogs.indexOf("blob-tests: reads begin");
    const errorMarker = browserLogs.indexOf("blob-tests: expected errors follow");
    expect(readMarker).toBeGreaterThanOrEqual(0);
    expect(errorMarker).toBeGreaterThan(readMarker);
    expect(browserLogs.slice(readMarker, errorMarker).filter((line) => line.includes("[E::"))).toEqual([]);
    expect(result.missingIndex).toContain("region query requires an index file");
    expect(result.revokedIndex).toContain("region query requires an index file");
    expect(result.revoked).toContain("read_bed: failed to open file");
    expect(result.unknown).toContain("read_bed: failed to open file");
    for (const error of result.registered) expect(error).toContain("read_bed: failed to open file");
    expect(result.alive).toEqual([{ answer: 42 }]);
  });
}

// Codex review of 9aed115: with enforceHostAllowlist, blob: URLs have an empty
// hostname, so no allowHosts entry could authorise a local file. The allowlist
// governs outbound hosts; blob: URLs make no network request, so they are exempt.
// duckdb-wasm keeps its Emscripten Module private, so the policy is set on the
// worker global. The first assertion proves it reached the extension: a same-origin
// HTTP read of a host outside allowHosts is refused.
test("host allowlist blocks non-listed hosts but not local blob: files", async ({ page }) => {
  const bed = await readFile(path.join(root, "test/data/fixture_mixed_regions.bed"), "utf8");
  const bedRows = bed.trim().split("\n").map((line) => {
    const [chrom, start, end, name] = line.split("\t");
    return { chrom, start: Number(start), end: Number(end), name };
  });
  await page.route("**/duckhts-package.js", (route) => route.fulfill({
    path: path.join(root, "js/src/local-file.js"), contentType: "text/javascript",
  }));
  await page.goto("/scripts/duckdb-wasm-local-test.html");
  const result = await page.evaluate(async ({ bed }) => {
    const duckdb = await import("/duckdb-browser.mjs");
    const { localFileUrl } = await import("/duckhts-package.js");
    const policy = { enforceHostAllowlist: true, allowHosts: ["example.org"] };
    const wrapper = localFileUrl(new Blob([`
      self.duckhtsWasmHttpConfig = ${JSON.stringify(policy)};
      importScripts(${JSON.stringify(new URL("/duckdb-browser-eh.worker.js", location.href).href)});
    `], { type: "text/javascript" }));
    const worker = new Worker(wrapper.url);
    const db = new duckdb.AsyncDuckDB(new duckdb.VoidLogger(), worker);
    const local = localFileUrl(new File([bed], "dropped.bed"));
    const httpFastq = new URL("/extdata/r1.fq", location.href).href;
    let conn;
    try {
      await db.instantiate(new URL("/duckdb-eh.wasm", location.href).href);
      await db.open({ allowUnsignedExtensions: true });
      conn = await db.connect();
      const extensionPath = new URL("/duckdb-wasm/duckhts.duckdb_extension.wasm", location.href).href;
      await db.registerFileURL(extensionPath, extensionPath, duckdb.DuckDBDataProtocol.HTTP, false);
      await conn.query(`LOAD '${extensionPath}'`);
      const rows = async (sql) => JSON.parse(JSON.stringify((await conn.query(sql)).toArray(),
        (_, value) => typeof value === "bigint" ? Number(value) : value));
      let blocked = null;
      try {
        await conn.query(`SELECT count(*) FROM read_fastq('${httpFastq}')`);
      } catch (e) {
        blocked = String(e);
      }
      const blobRows = await rows(`SELECT chrom, start, "end", name FROM read_bed('${local.url}')`);
      return { blocked, blobRows };
    } finally {
      if (conn) await conn.close();
      await db.terminate();
      worker.terminate();
      wrapper.revoke();
      local.revoke();
    }
  }, { bed });

  expect(result.blocked).toContain("r1.fq");
  expect(result.blobRows).toEqual(bedRows);
});
