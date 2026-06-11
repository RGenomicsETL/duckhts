import { defineConfig, devices } from "@playwright/test";
import path from "node:path";

// The wasm site is assembled by `SERVE=0 scripts/start_duckdb_wasm_local_test.sh`,
// which prints DUCKHTS_WASM_SITE_ROOT / DUCKHTS_WASM_PORT.  The Makefile target
// and CI job export those; fall back to the script's defaults otherwise.
const SITE_ROOT =
  process.env.DUCKHTS_WASM_SITE_ROOT ||
  path.resolve(__dirname, "../../.duckdb-wasm-local-artifacts/site");
const PORT = process.env.DUCKHTS_WASM_PORT || "8001";
const HOST = process.env.DUCKHTS_WASM_HOST || "127.0.0.1";
const BASE_URL = `http://${HOST}:${PORT}`;

export default defineConfig({
  testDir: ".",
  fullyParallel: false,
  forbidOnly: !!process.env.CI,
  retries: 0,
  workers: 1,
  reporter: process.env.CI ? [["html", { open: "never" }], ["list"]] : "list",
  timeout: 180_000,
  use: {
    baseURL: BASE_URL,
    headless: true,
    trace: "retain-on-failure",
  },
  projects: [{ name: "chromium", use: { ...devices["Desktop Chrome"] } }],
  webServer: {
    command: `python3 -m http.server ${PORT} --bind ${HOST}`,
    cwd: SITE_ROOT,
    url: `${BASE_URL}/scripts/duckdb-wasm-local-test.html`,
    reuseExistingServer: !process.env.CI,
    timeout: 60_000,
  },
});
