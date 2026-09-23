import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { test } from "node:test";

import * as entry from "../src/index.js";

// localFileUrl only works with DuckHTS builds that contain the blob: handler.
// Until artifacts.json pins such a release, the package must not export it
// (Codex review P1 on https://github.com/RGenomicsETL/duckhts/pull/248). When
// the pinned release includes the handler, export it and flip this test.
const BLOB_CAPABLE_FROM = "1.5.3";

test("localFileUrl is exported only once the pinned builds support blob: URLs", async () => {
  const manifest = JSON.parse(await readFile(new URL("../artifacts.json", import.meta.url), "utf8"));
  const pinned = manifest.duckhts.split(".").map(Number);
  const capable = BLOB_CAPABLE_FROM.split(".").map(Number);
  const pinsCapable = pinned[0] > capable[0] || (pinned[0] === capable[0] &&
    (pinned[1] > capable[1] || (pinned[1] === capable[1] && pinned[2] >= capable[2])));
  assert.equal("localFileUrl" in entry, pinsCapable, `artifacts.json pins DuckHTS ${manifest.duckhts}`);
});
