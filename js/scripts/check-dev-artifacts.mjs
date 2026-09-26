// CI decides whether the dev browser/pack steps can run; publish never skips them.
import { appendFile, readFile } from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { checkDevArtifacts, UnavailableDevArtifactsError } from "./dev-artifacts.mjs";

const mode = process.argv[2];
if (!["pr", "publish"].includes(mode)) {
  throw new Error("Usage: node scripts/check-dev-artifacts.mjs pr|publish [manifest.json]");
}
const packageRoot = path.dirname(path.dirname(fileURLToPath(import.meta.url)));
const manifestPath = process.argv[3] ?? path.join(packageRoot, "artifacts-dev.json");
const manifest = JSON.parse(await readFile(manifestPath, "utf8"));
if (manifest.signed !== false) {
  throw new Error(`${manifestPath}: expected an unsigned dev manifest`);
}
const names = Object.values(manifest.platforms).map((pin) => pin.artifact);
let stage = true;
try {
  checkDevArtifacts(manifest, [...new Set(names)]);
} catch (error) {
  if (mode !== "pr" || !(error instanceof UnavailableDevArtifactsError)) throw error;
  console.log(`::notice::${error.message} Skipping dev browser and pack steps on this pull request.`);
  stage = false;
}
if (process.env.GITHUB_OUTPUT) {
  await appendFile(process.env.GITHUB_OUTPUT, `stage=${stage}\n`);
}
