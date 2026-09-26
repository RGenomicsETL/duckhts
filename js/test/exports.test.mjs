import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { test } from "node:test";
import { localFileUrl, SIGNED } from "../src/index.js";

test("localFileUrl creates object URLs only for builds with the blob handler", async () => {
  const manifest = JSON.parse(await readFile(new URL(SIGNED ? "../artifacts.json" :
    "../artifacts-dev.json", import.meta.url), "utf8"));
  const [major, minor, patch] = manifest.duckhts.split(".").map(Number);
  const capable = !SIGNED || major > 1 || (major === 1 && (minor > 5 || (minor === 5 && patch >= 3)));
  const file = new File(["chr1\t0\t10\n"], "dropped.bed");
  if (!capable) {
    assert.throws(() => localFileUrl(file), /signed build does not support blob:/);
    return;
  }
  const local = localFileUrl(file);
  try {
    assert.match(local.url, /^blob:/);
    assert.equal(await (await fetch(local.url)).text(), await file.text());
  } finally {
    local.revoke();
  }
});
