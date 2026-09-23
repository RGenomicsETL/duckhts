import { readFile } from "node:fs/promises";
import { expect, test } from "@playwright/test";

// The package's loader and staged binaries in a real browser, served from one
// origin with unsigned extensions disallowed.  Expected rows come from the
// fixture files themselves, not from DuckHTS output.

const BED_FIXTURE = new URL("../../../test/data/fixture_mixed_regions.bed", import.meta.url);

async function bedFixtureRows() {
  const text = await readFile(BED_FIXTURE, "utf8");
  return text
    .trim()
    .split("\n")
    .map((line) => {
      const [chrom, start, end, name] = line.split("\t");
      return { chrom, start: Number(start), end: Number(end), name };
    });
}

async function probe(page, runtime) {
  await page.goto("/probe.html");
  await expect(page).toHaveTitle("ready");
  return page.evaluate((name) => window.probe(name), runtime);
}

test("DuckDB v1.5 runtime loads the signed binary and reads over same-origin HTTP", async ({ page }) => {
  const result = await probe(page, "next");

  expect(result.duckdb).toMatch(/^v1\.5\./);
  expect(result.unsignedAllowed).toBe(false);
  expect(result.load).toEqual({ ok: true, value: { platform: "wasm_eh" } });
  expect(result.bedHttp).toEqual({ ok: true, value: await bedFixtureRows() });
  expect(result.gffHttp).toEqual({
    ok: true,
    value: [{ seqname: "chr1", feature: "gene", start: 1, end: 10 }],
  });
});

// Files registered with duckdb-wasm live in its virtual file system, which
// htslib does not open.  This records the current contract; when
// https://github.com/RGenomicsETL/duckhts/issues/246 gains support, this test
// fails and must be turned into a positive check.
test("registered duckdb-wasm files are not visible to DuckHTS readers", async ({ page }) => {
  const result = await probe(page, "next");

  expect(result.load.ok).toBe(true);
  expect(result.bedRegisteredBuffer.ok).toBe(false);
  expect(result.bedRegisteredBuffer.error).toContain("read_bed: failed to open file");
  expect(result.bedRegisteredHandle.ok).toBe(false);
  expect(result.bedRegisteredHandle.error).toContain("read_bed: failed to open file");
});

test("DuckDB v1.4 runtime loads the binary", async ({ page }) => {
  test.fail(true, "init registration order: https://github.com/RGenomicsETL/duckhts/issues/247");
  const result = await probe(page, "stable");

  expect(result.duckdb).toMatch(/^v1\.4\./);
  expect(result.load).toEqual({ ok: true, value: { platform: "wasm_eh" } });
});
