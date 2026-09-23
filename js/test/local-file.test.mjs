import assert from "node:assert/strict";
import { test } from "node:test";
import { localFileUrl } from "../src/local-file.js";

test("localFileUrl exposes the File bytes without changing the File", async () => {
  const file = new File(["chr1\t0\t10\n"], "dropped.bed");
  const local = localFileUrl(file);
  try {
    assert.match(local.url, /^blob:/);
    assert.equal(await (await fetch(local.url)).text(), await file.text());
    assert.equal(file.name, "dropped.bed");
  } finally {
    local.revoke();
  }
});

test("each URL has an independent lifetime, and revoke is idempotent", async () => {
  const blob = new Blob([Uint8Array.of(0, 1, 254, 255)]);
  const first = localFileUrl(blob);
  const second = localFileUrl(blob);
  try {
    assert.notEqual(first.url, second.url);
    const { revoke } = first;
    revoke();
    revoke();
    await assert.rejects(fetch(first.url), TypeError);
    assert.deepEqual(await (await fetch(second.url)).arrayBuffer(), await blob.arrayBuffer());
  } finally {
    first.revoke();
    second.revoke();
  }
  await assert.rejects(fetch(second.url), TypeError);
});

test("localFileUrl delegates invalid inputs to createObjectURL", () => {
  assert.throws(() => localFileUrl("dropped.bed"), TypeError);
  assert.throws(() => localFileUrl(undefined), TypeError);
});
