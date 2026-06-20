import { test, expect } from "@playwright/test";

// Drives the existing manual harness (scripts/duckdb-wasm-local-test.html)
// headlessly: instantiate duckdb-wasm, LOAD the DuckHTS wasm extension, assert
// the SIMD dispatch kernels resolve on wasm (scalar fallback), and run the
// same-origin HTTP smoke queries.  The harness sets #status to Passed/Failed.
test("DuckHTS wasm extension loads and SIMD kernels resolve (scalar fallback)", async ({ page }) => {
  const pageErrors: string[] = [];
  page.on("pageerror", (e) => pageErrors.push(String(e)));

  await page.goto("/scripts/duckdb-wasm-local-test.html");
  await page.click("#run");

  // wasm instantiation + extension load can take a while; wait for a terminal
  // state rather than a fixed delay.
  await expect(page.locator("#status")).toHaveText(/Passed|Failed/, { timeout: 150_000 });

  const status = await page.locator("#status").textContent();
  const log = await page.locator("#log").textContent();
  expect(
    status,
    `harness #status=${status}\n--- harness log ---\n${log}\n--- page errors ---\n${pageErrors.join("\n")}`,
  ).toBe("Passed");
});
