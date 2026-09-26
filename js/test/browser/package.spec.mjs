import { readFile } from "node:fs/promises";
import { expect, test } from "@playwright/test";
import { SIGNED } from "../../src/index.js";

// Test the staged binary on the selected channel, using independent fixture rows.

const BED_FIXTURE = new URL("../../../test/data/fixture_mixed_regions.bed", import.meta.url);
const signedManifest = JSON.parse(await readFile(new URL("../../artifacts.json", import.meta.url), "utf8"));

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

async function probe(page, runtime, allowUnsigned = !SIGNED) {
  await page.goto("/probe.html");
  await expect(page).toHaveTitle("ready");
  return page.evaluate(({ name, allowUnsigned }) => window.probe(name, allowUnsigned),
    { name: runtime, allowUnsigned });
}

test("DuckDB v1.5 loads the selected binary and reads over HTTP and blob where supported", async ({ page }) => {
  const result = await probe(page, "next");

  expect(result.duckdb).toMatch(/^v1\.5\./);
  expect(result.unsignedAllowed).toBe(!SIGNED);
  expect(result.load).toEqual({ ok: true, value: { platform: "wasm_eh" } });
  expect(result.bedHttp).toEqual({ ok: true, value: await bedFixtureRows() });
  expect(result.gffHttp).toEqual({
    ok: true,
    value: [{ seqname: "chr1", feature: "gene", start: 1, end: 10 }],
  });
  if (!SIGNED) expect(result.bedBlob).toEqual({ ok: true, value: await bedFixtureRows() });
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

test("dev loader rejects an unsigned-disallowed connection before LOAD", async ({ page }) => {
  test.skip(SIGNED, "signed builds do not need an unsigned gate");
  const result = await probe(page, "next", false);
  expect(result.unsignedAllowed).toBe(false);
  expect(result.load).toEqual({ ok: false,
    error: expect.stringContaining("dev channel requires allow_unsigned_extensions=true") });
});

test("DuckDB v1.4 runtime loads the binary", async ({ page }) => {
  test.fail(SIGNED && signedManifest.duckhts === "1.5.2",
    "signed DuckHTS 1.5.2 predates the DuckDB 1.4 initialization fix");
  const result = await probe(page, "stable");

  expect(result.duckdb).toMatch(/^v1\.4\./);
  expect(result.load).toEqual({ ok: true, value: { platform: "wasm_eh" } });
});
