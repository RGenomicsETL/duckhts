import assert from "node:assert/strict";
import { execFile } from "node:child_process";
import { createHash } from "node:crypto";
import { chmod, mkdir, mkdtemp, readFile, rm, stat, writeFile } from "node:fs/promises";
import { createServer } from "node:http";
import { tmpdir } from "node:os";
import path from "node:path";
import { after, before, test } from "node:test";
import { fileURLToPath } from "node:url";
import { promisify } from "node:util";

const stageScript = fileURLToPath(new URL("../scripts/stage.mjs", import.meta.url));
const run = promisify(execFile);
const platforms = ["wasm_mvp", "wasm_eh", "wasm_threads"];
const fileName = "duckhts.duckdb_extension.wasm";
const payload = Buffer.from("\0asm stand-in for a staged extension");
const payloadSha256 = createHash("sha256").update(payload).digest("hex");

let server;
let source;
let workDir;

before(async () => {
  server = createServer((request, response) => {
    if (platforms.some((platform) => request.url === `/${platform}/${fileName}`)) {
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

async function stage(channel, name, badPlatform, options = {}) {
  const manifestPath = path.join(workDir, `${name}.json`);
  const output = path.join(workDir, name);
  const directory = path.join(workDir, `${name}-downloads`);
  const dev = channel === "dev";
  const pins = Object.fromEntries(platforms.map((platform) => [platform, dev ? {
    artifact: `extension-${platform}`,
    sha256: platform === badPlatform ? "0".repeat(64) : payloadSha256,
  } : platform === badPlatform ? "0".repeat(64) : payloadSha256]));
  const manifest = {
    signed: !dev,
    source: dev ? { repository: "RGenomicsETL/duckhts", run: 1 } : source,
    platforms: pins,
  };
  await writeFile(manifestPath, JSON.stringify(manifest));
  if (dev && !options.gh) {
    for (const { artifact } of Object.values(pins)) {
      await mkdir(path.join(directory, artifact), { recursive: true });
      if (!options.missing) await writeFile(path.join(directory, artifact, fileName), payload);
    }
  }
  for (const platform of options.cached ?? []) {
    const target = path.join(output, platform, fileName);
    await mkdir(path.dirname(target), { recursive: true });
    await writeFile(target, platform === options.stale ? Buffer.from("stale") : payload);
  }
  const args = [stageScript, channel, manifestPath, output, ...(dev && !options.gh ? [directory] : [])];
  const result = await run(process.execPath, args, {
    env: { ...process.env, DUCKHTS_NPM_CHANNEL: channel, ...options.env },
  }).catch((error) => error);
  return { result, output };
}

for (const channel of ["signed", "dev"]) {
  test(`${channel}: stages each pinned binary byte for byte${channel === "dev" ? " from an offline multi-artifact directory" : ""}`, async () => {
    const { result, output } = await stage(channel, `${channel}-match`);
    assert.equal(result.code, undefined, result.stderr);
    for (const platform of platforms) {
      assert.deepEqual(await readFile(path.join(output, platform, fileName)), payload);
    }
  });

  test(`${channel}: rejects a checksum mismatch without writing any binary`, async () => {
    const { result, output } = await stage(channel, `${channel}-mismatch`, "wasm_threads");
    assert.equal(result.code, 1);
    assert.match(result.stderr, /sha256 [0-9a-f]{64}, expected 0{64}/);
    for (const platform of platforms) {
      await assert.rejects(stat(path.join(output, platform, fileName)), { code: "ENOENT" });
    }
  });
}

async function mockGh(name) {
  const bin = path.join(workDir, `${name}-bin`);
  const log = path.join(workDir, `${name}-gh.log`);
  const payloadFile = path.join(workDir, `${name}-payload.wasm`);
  await mkdir(bin);
  await writeFile(payloadFile, payload);
  const executable = path.join(bin, "gh");
  await writeFile(executable, `#!/bin/sh
set -eu
[ "$#" -eq 9 ]
[ "$1" = run ] && [ "$2" = download ] && [ "$3" = 1 ]
[ "$4" = -R ] && [ "$5" = RGenomicsETL/duckhts ]
[ "$6" = -n ] && [ "$8" = -D ]
printf '%s\\n' "$7" >> "$GH_LOG"
mkdir -p "$9"
cp "$GH_PAYLOAD" "$9/${fileName}"
`);
  await chmod(executable, 0o755);
  return { log, env: { PATH: `${bin}${path.delimiter}${process.env.PATH}`,
    GH_LOG: log, GH_PAYLOAD: payloadFile } };
}

test("dev: downloads only a stale artifact beside two verified cached platforms", async () => {
  const gh = await mockGh("stale");
  const { result, output } = await stage("dev", "dev-stale", null, {
    gh: true, cached: platforms, stale: "wasm_eh", env: gh.env,
  });
  assert.equal(result.code, undefined, result.stderr);
  assert.deepEqual((await readFile(gh.log, "utf8")).trim().split("\n"), ["extension-wasm_eh"]);
  assert.match(result.stdout, /wasm_mvp: verified/);
  assert.match(result.stdout, /wasm_threads: verified/);
  for (const platform of platforms) {
    assert.deepEqual(await readFile(path.join(output, platform, fileName)), payload);
  }
});

test("dev: stages multiple missing artifacts through single-artifact gh downloads", async () => {
  const gh = await mockGh("all");
  const { result, output } = await stage("dev", "dev-all", null, { gh: true, env: gh.env });
  assert.equal(result.code, undefined, result.stderr);
  assert.deepEqual((await readFile(gh.log, "utf8")).trim().split("\n"),
    platforms.map((platform) => `extension-${platform}`));
  for (const platform of platforms) {
    assert.deepEqual(await readFile(path.join(output, platform, fileName)), payload);
  }
});

test("dev: rejects a missing artifact without fetching another channel", async () => {
  const { result } = await stage("dev", "dev-missing", null, { missing: true });
  assert.equal(result.code, 1);
  assert.match(result.stderr, /ENOENT/);
});

test("stage rejects a manifest from the other channel", async () => {
  const manifestPath = path.join(workDir, "wrong-channel.json");
  await writeFile(manifestPath, JSON.stringify({ signed: false, platforms: {} }));
  const result = await run(process.execPath, [stageScript, "signed", manifestPath, path.join(workDir, "wrong")], {
    env: { ...process.env, DUCKHTS_NPM_CHANNEL: "signed" },
  }).catch((error) => error);
  assert.equal(result.code, 1);
  assert.match(result.stderr, /signed status disagrees/);
});
