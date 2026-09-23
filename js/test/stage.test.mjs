import assert from "node:assert/strict";
import { execFile } from "node:child_process";
import { createHash } from "node:crypto";
import { mkdtemp, readFile, rm, stat, writeFile } from "node:fs/promises";
import { createServer } from "node:http";
import { tmpdir } from "node:os";
import path from "node:path";
import { after, before, test } from "node:test";
import { fileURLToPath } from "node:url";
import { promisify } from "node:util";

// Runs scripts/stage.mjs against a local server, so the checksum contract is
// tested without touching the community repository.

const stageScript = fileURLToPath(new URL("../scripts/stage.mjs", import.meta.url));
const run = promisify(execFile);
const payload = Buffer.from("\0asm stand-in for a staged extension");
const payloadSha256 = createHash("sha256").update(payload).digest("hex");

let server;
let source;
let workDir;

before(async () => {
  server = createServer((request, response) => {
    if (request.url === "/wasm_eh/duckhts.duckdb_extension.wasm") {
      response.writeHead(200).end(payload);
    } else {
      response.writeHead(404).end();
    }
  });
  await new Promise((resolve) => server.listen(0, "127.0.0.1", resolve));
  source = `http://127.0.0.1:${server.address().port}`;
  workDir = await mkdtemp(path.join(tmpdir(), "duckhts-stage-"));
});

after(async () => {
  server.close();
  await rm(workDir, { recursive: true, force: true });
});

async function stage(name, sha256) {
  const manifest = path.join(workDir, `${name}.json`);
  const output = path.join(workDir, name);
  await writeFile(manifest, JSON.stringify({ source, platforms: { wasm_eh: sha256 } }));
  const result = await run(process.execPath, [stageScript, manifest, output]).catch((error) => error);
  return { result, target: path.join(output, "wasm_eh", "duckhts.duckdb_extension.wasm") };
}

test("stage writes a download whose sha256 matches the manifest", async () => {
  const { result, target } = await stage("match", payloadSha256);
  assert.equal(result.code, undefined, result.stderr);
  assert.deepEqual(await readFile(target), payload);
});

test("stage rejects a checksum mismatch and writes nothing", async () => {
  const { result, target } = await stage("mismatch", "0".repeat(64));
  assert.equal(result.code, 1);
  assert.match(result.stderr, /sha256 [0-9a-f]{64}, expected 0{64}/);
  await assert.rejects(stat(target), { code: "ENOENT" });
});

test("stage rejects a missing artifact", async () => {
  const manifest = path.join(workDir, "missing.json");
  await writeFile(manifest, JSON.stringify({ source, platforms: { wasm_mvp: payloadSha256 } }));
  const result = await run(process.execPath, [stageScript, manifest, path.join(workDir, "missing")]).catch(
    (error) => error,
  );
  assert.equal(result.code, 1);
  assert.match(result.stderr, /HTTP 404/);
});
