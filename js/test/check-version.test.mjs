import assert from "node:assert/strict";
import { execFile } from "node:child_process";
import { copyFile, mkdir, mkdtemp, rm, writeFile } from "node:fs/promises";
import { tmpdir } from "node:os";
import path from "node:path";
import { test } from "node:test";
import { fileURLToPath } from "node:url";
import { promisify } from "node:util";

const script = fileURLToPath(new URL("../scripts/check-version.mjs", import.meta.url));
const run = promisify(execFile);

test("dev version check rejects source-only development version bumps", async () => {
  const root = await mkdtemp(path.join(tmpdir(), "duckhts-version-"));
  const js = path.join(root, "js");
  await mkdir(path.join(js, "scripts"), { recursive: true });
  try {
    await copyFile(script, path.join(js, "scripts", "check-version.mjs"));
    await writeFile(path.join(js, "package.json"), JSON.stringify({ version: "1.5.2-9001" }));
    await writeFile(path.join(js, "package-lock.json"), JSON.stringify({
      version: "1.5.2-9001", packages: { "": { version: "1.5.2-9001" } },
    }));
    const manifest = path.join(js, "artifacts-dev.json");
    await writeFile(manifest, JSON.stringify({ duckhts: "1.5.2.9001" }));
    const description = path.join(root, "description.yml");
    await writeFile(description, "version: 1.5.2.9001\n");
    const command = [path.join(js, "scripts", "check-version.mjs"), "dev"];
    const env = { ...process.env, DUCKHTS_NPM_CHANNEL: "dev" };
    const valid = await run(process.execPath, command, { env });
    assert.match(valid.stdout, /dev: DuckHTS 1\.5\.2\.9001 -> npm 1\.5\.2-9001/);

    await writeFile(description, "version: 1.5.2.9002\n");
    const staleManifest = await run(process.execPath, command, { env }).catch((error) => error);
    assert.equal(staleManifest.code, 1);
    assert.match(staleManifest.stderr, /dev manifest 1\.5\.2\.9001 disagrees with description\.yml 1\.5\.2\.9002/);

    await writeFile(manifest, JSON.stringify({ duckhts: "1.5.2.9002" }));
    const stalePackage = await run(process.execPath, command, { env }).catch((error) => error);
    assert.equal(stalePackage.code, 1);
    assert.match(stalePackage.stderr, /dev npm version must be 1\.5\.2-9002/);
  } finally {
    await rm(root, { recursive: true, force: true });
  }
});
