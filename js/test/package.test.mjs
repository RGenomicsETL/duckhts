import assert from "node:assert/strict";
import { execFileSync } from "node:child_process";
import { test } from "node:test";
import { fileURLToPath } from "node:url";

// The published tarball must carry the licence and the notices for everything
// linked into the wasm binaries (https://github.com/RGenomicsETL/duckhts/pull/248).
test("npm pack includes LICENSE and THIRD_PARTY_NOTICES.md", () => {
  const packageRoot = fileURLToPath(new URL("..", import.meta.url));
  const [pack] = JSON.parse(execFileSync("npm", ["pack", "--dry-run", "--ignore-scripts", "--json"], {
    cwd: packageRoot, encoding: "utf8",
  }));
  const files = pack.files.map((f) => f.path);
  assert.ok(files.includes("LICENSE"), files.join(", "));
  assert.ok(files.includes("THIRD_PARTY_NOTICES.md"), files.join(", "));
});
