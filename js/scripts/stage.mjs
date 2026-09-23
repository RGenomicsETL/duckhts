// Stage the signed DuckHTS duckdb-wasm binaries named by artifacts.json into
// dist/<platform>/duckhts.duckdb_extension.wasm.  This is the package's only
// network step; every download must match its pinned sha256 byte for byte, so
// the published files are exactly the community-repository artifacts.
//
// Usage: node scripts/stage.mjs [manifest.json] [output-dir]
// (defaults: the package's artifacts.json and dist/).
import { createHash } from "node:crypto";
import { mkdir, readFile, rename, rm, writeFile } from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";

const packageRoot = path.dirname(path.dirname(fileURLToPath(import.meta.url)));
const manifestPath = process.argv[2] ?? path.join(packageRoot, "artifacts.json");
const outputDir = process.argv[3] ?? path.join(packageRoot, "dist");
const manifest = JSON.parse(await readFile(manifestPath, "utf8"));
const fileName = "duckhts.duckdb_extension.wasm";

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

for (const [platform, expected] of Object.entries(manifest.platforms)) {
  const target = path.join(outputDir, platform, fileName);
  const cached = await readIfPresent(target);
  if (cached && sha256(cached) === expected) {
    console.log(`${platform}: verified`);
    continue;
  }

  const url = `${manifest.source}/${platform}/${fileName}`;
  const response = await fetch(url);
  if (!response.ok) {
    throw new Error(`${url}: HTTP ${response.status}`);
  }
  const bytes = Buffer.from(await response.arrayBuffer());
  const actual = sha256(bytes);
  if (actual !== expected) {
    throw new Error(`${url}: sha256 ${actual}, expected ${expected}`);
  }

  await mkdir(path.dirname(target), { recursive: true });
  const partial = `${target}.partial`;
  await writeFile(partial, bytes);
  await rm(target, { force: true });
  await rename(partial, target);
  console.log(`${platform}: staged ${bytes.length} bytes`);
}
