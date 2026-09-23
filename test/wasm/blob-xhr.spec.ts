import { test, expect } from "@playwright/test";

// Exercise the browser transport independently of DuckHTS. The object URL is
// created on the page, as it would be for a File from a drop event.
test("worker synchronous XHR blob transport", async ({ page, browser }) => {
  await page.goto("/scripts/duckdb-wasm-local-test.html");
  const measured = await page.evaluate(async () => {
    function probeWorker() {
      // Revocation reaches other threads asynchronously, so a request issued
      // right after revokeObjectURL() on the page may still succeed. Before the
      // revoked probe, wait (bounded) until a ranged GET is refused.
      const refused = (url) => {
        const xhr = new XMLHttpRequest();
        try {
          xhr.open("GET", url, false);
          xhr.setRequestHeader("Range", "bytes=0-0");
          xhr.send();
          return false;
        } catch (e) {
          return true;
        }
      };
      self.onmessage = async ({ data: { url, awaitRevocation } }) => {
        if (awaitRevocation) {
          let attempts = 0;
          while (!refused(url) && attempts < 200) {
            attempts += 1;
            await new Promise((resolve) => setTimeout(resolve, 10));
          }
        }
        const results = [];
        for (const [method, range] of [
          ["HEAD", null],
          ["GET", "bytes=0-0"],
          ["GET", "bytes=4-7"],
          ["GET", "bytes=14-31"],
          ["GET", "bytes=16-31"],
        ]) {
          const xhr = new XMLHttpRequest();
          let error = null;
          try {
            xhr.open(method!, url, false);
            if (range) xhr.setRequestHeader("Range", range);
            xhr.responseType = "arraybuffer";
            xhr.send();
          } catch (e) {
            error = (e as Error).name;
          }
          results.push({
            method, range, status: xhr.status, error,
            length: xhr.getResponseHeader("Content-Length"),
            contentRange: xhr.getResponseHeader("Content-Range"),
            bytes: Array.from(new Uint8Array(xhr.response || new ArrayBuffer(0))),
          });
        }
        self.postMessage(results);
      };
    }
    const workerUrl = URL.createObjectURL(new Blob([`(${probeWorker.toString()})()`]));
    const worker = new Worker(workerUrl);
    const url = URL.createObjectURL(new File([Uint8Array.from({ length: 16 }, (_, i) => i)], "probe.bin"));
    const request = (awaitRevocation) => new Promise((resolve, reject) => {
      worker.onmessage = ({ data }) => resolve(data);
      worker.onerror = reject;
      worker.postMessage({ url, awaitRevocation });
    });
    try {
      const live = await request(false);
      URL.revokeObjectURL(url);
      const revoked = await request(true);
      return { live, revoked };
    } finally {
      worker.terminate();
      URL.revokeObjectURL(workerUrl);
      URL.revokeObjectURL(url);
    }
  });
  console.log(JSON.stringify({ chromium: browser.version(), ...measured }));
  expect(measured.live).toEqual([
    { method: "HEAD", range: null, status: 0, error: "NetworkError", length: null, contentRange: null, bytes: [] },
    { method: "GET", range: "bytes=0-0", status: 206, error: null, length: "1", contentRange: "bytes 0-0/16", bytes: [0] },
    { method: "GET", range: "bytes=4-7", status: 206, error: null, length: "4", contentRange: "bytes 4-7/16", bytes: [4, 5, 6, 7] },
    { method: "GET", range: "bytes=14-31", status: 206, error: null, length: "2", contentRange: "bytes 14-15/16", bytes: [14, 15] },
    { method: "GET", range: "bytes=16-31", status: 0, error: "NetworkError", length: null, contentRange: null, bytes: [] },
  ]);
  expect(measured.revoked).toEqual(measured.live.map(({ method, range }) => ({
    method, range, status: 0, error: "NetworkError", length: null, contentRange: null, bytes: [],
  })));
});
