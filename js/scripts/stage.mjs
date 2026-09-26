// Stage one explicitly selected manifest. Network access is confined to this
// script; neither a failed download nor a checksum mismatch can be packaged.
// Usage: node scripts/stage.mjs signed|dev [manifest.json] [output-dir] [download-dir]
// A download-dir is an offline gh-run-download tree for dev staging.
import { execFileSync } from "node:child_process";
import { createHash } from "node:crypto";
import { mkdir, mkdtemp, readFile, rename, rm, writeFile } from "node:fs/promises";
import { tmpdir } from "node:os";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { writeChannel } from "./channel.mjs";
import { checkDevArtifacts } from "./dev-artifacts.mjs";

const packageRoot = path.dirname(path.dirname(fileURLToPath(import.meta.url)));
const channel = process.argv[2] ?? process.env.DUCKHTS_NPM_CHANNEL;
if (!["signed", "dev"].includes(channel)) {
  throw new Error("Select DUCKHTS_NPM_CHANNEL=signed|dev (or pass signed|dev as the first argument)");
}
if (process.argv[2] && process.env.DUCKHTS_NPM_CHANNEL && process.argv[2] !== process.env.DUCKHTS_NPM_CHANNEL) {
  throw new Error("Channel argument and DUCKHTS_NPM_CHANNEL disagree");
}
const manifestPath = process.argv[3] ?? path.join(packageRoot,
  channel === "dev" ? "artifacts-dev.json" : "artifacts.json");
const outputDir = process.argv[4] ?? path.join(packageRoot, "dist");
const downloadDir = process.argv[5];
if (process.argv[3] && !process.argv[4]) {
  throw new Error("A custom manifest requires an explicit output directory");
}
const manifest = JSON.parse(await readFile(manifestPath, "utf8"));
const fileName = "duckhts.duckdb_extension.wasm";
const platforms = ["wasm_mvp", "wasm_eh", "wasm_threads"];

if (channel === "dev" ? manifest.signed !== false :
    manifest.signed !== true && manifest.signed !== undefined) {
  throw new Error(`${manifestPath}: signed status disagrees with ${channel} channel`);
}
if (Object.keys(manifest.platforms).sort().join() !== [...platforms].sort().join()) {
  throw new Error(`${manifestPath}: expected exactly ${platforms.join(", ")}`);
}

function sha256(bytes) {
  return createHash("sha256").update(bytes).digest("hex");
}

async function readIfPresent(file) {
  try {
    return await readFile(file);
  } catch (error) {
    if (error.code === "ENOENT") return null;
    throw error;
  }
}

const pending = [];
for (const platform of platforms) {
  const pin = manifest.platforms[platform];
  const expected = channel === "dev" ? pin?.sha256 : pin;
  if (typeof expected !== "string" || !/^[0-9a-f]{64}$/.test(expected)) {
    throw new Error(`${manifestPath}: invalid sha256 for ${platform}`);
  }
  if (channel === "dev" && (typeof pin.artifact !== "string" ||
      !/^[\w.-]+$/.test(pin.artifact) || pin.artifact === ".." || pin.artifact === ".")) {
    throw new Error(`${manifestPath}: invalid artifact for ${platform}`);
  }
  const target = path.join(outputDir, platform, fileName);
  const cached = await readIfPresent(target);
  if (cached && sha256(cached) === expected) {
    console.log(`${platform}: verified`);
  } else {
    pending.push({ platform, expected, target, pin });
  }
}

let directory = downloadDir;
let temporary = false;
try {
  if (channel === "dev" && pending.length && !directory) {
    directory = await mkdtemp(path.join(tmpdir(), "duckhts-npm-dev-"));
    temporary = true;
    const artifacts = [...new Set(pending.map(({ pin }) => pin.artifact))];
    checkDevArtifacts(manifest, artifacts);
    for (const artifact of artifacts) {
      const destination = path.join(directory, artifact);
      await mkdir(destination, { recursive: true });
      execFileSync("gh", ["run", "download", String(manifest.source.run), "-R", manifest.source.repository,
        "-n", artifact, "-D", destination], { stdio: "inherit" });
    }
  }

  const verified = [];
  for (const { platform, expected, target, pin } of pending) {
    let bytes;
    let source;
    if (channel === "dev") {
      source = path.join(directory, pin.artifact, fileName);
      bytes = await readFile(source);
    } else {
      source = `${manifest.source}/${platform}/${fileName}`;
      const response = await fetch(source);
      if (!response.ok) throw new Error(`${source}: HTTP ${response.status}`);
      bytes = Buffer.from(await response.arrayBuffer());
    }
    const actual = sha256(bytes);
    if (actual !== expected) throw new Error(`${source}: sha256 ${actual}, expected ${expected}`);
    verified.push({ platform, target, bytes });
  }

  for (const { platform, target, bytes } of verified) {
    await mkdir(path.dirname(target), { recursive: true });
    const partial = `${target}.partial`;
    await writeFile(partial, bytes);
    await rename(partial, target);
    console.log(`${platform}: staged ${bytes.length} bytes`);
  }
  if (!process.argv[4]) await writeChannel(channel);
} finally {
  if (temporary) await rm(directory, { recursive: true, force: true });
}
