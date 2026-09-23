import { defineConfig, devices } from "@playwright/test";

const PORT = process.env.DUCKHTS_NPM_TEST_PORT || "8765";
const BASE_URL = `http://127.0.0.1:${PORT}`;

export default defineConfig({
  testDir: ".",
  testMatch: "*.spec.mjs",
  fullyParallel: false,
  forbidOnly: !!process.env.CI,
  retries: 0,
  workers: 1,
  reporter: process.env.CI ? [["html", { open: "never" }], ["list"]] : "list",
  timeout: 120_000,
  use: {
    baseURL: BASE_URL,
    headless: true,
    trace: "retain-on-failure",
    // --no-sandbox: required when Chromium runs as root (dev containers, some CI).
    launchOptions: { args: ["--no-sandbox"] },
  },
  projects: [{ name: "chromium", use: { ...devices["Desktop Chrome"] } }],
  webServer: {
    command: "node serve.mjs",
    env: { PORT },
    url: `${BASE_URL}/probe.html`,
    reuseExistingServer: !process.env.CI,
    timeout: 30_000,
  },
});
