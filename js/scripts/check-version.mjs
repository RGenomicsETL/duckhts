// Check the npm release identity against its pinned extension version.
import { readFile } from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";

const root = path.dirname(path.dirname(fileURLToPath(import.meta.url)));
const channel = process.argv[2] ?? process.env.DUCKHTS_NPM_CHANNEL;
if (channel !== "signed" && channel !== "dev") throw new Error("Select signed or dev");

const packageJson = JSON.parse(await readFile(path.join(root, "package.json"), "utf8"));
const lock = JSON.parse(await readFile(path.join(root, "package-lock.json"), "utf8"));
const manifest = JSON.parse(await readFile(path.join(root,
  channel === "dev" ? "artifacts-dev.json" : "artifacts.json"), "utf8"));
const description = await readFile(path.join(root, "../description.yml"), "utf8");
const sourceVersion = /^version:\s*([0-9]+\.[0-9]+\.[0-9]+\.[0-9]+)\s*$/m.exec(description)?.[1];

if (channel === "dev" && (!sourceVersion || manifest.duckhts !== sourceVersion)) {
  throw new Error(`dev manifest ${manifest.duckhts} disagrees with description.yml ${sourceVersion}`);
}
const development = channel === "dev" && /^(\d+\.\d+\.\d+)\.(9\d+)$/.exec(sourceVersion);
if (channel === "dev" && !development) {
  throw new Error(`Invalid dev version in description.yml: ${sourceVersion}`);
}
const expected = channel === "dev" ? `${development[1]}-${development[2]}` : manifest.duckhts;
if (channel === "signed" && !/^\d+\.\d+\.\d+$/.test(expected)) {
  throw new Error(`Invalid signed version in artifacts.json: ${expected}`);
}
if (packageJson.version !== expected || lock.version !== expected || lock.packages[""].version !== expected) {
  throw new Error(`${channel} npm version must be ${expected} in package.json and package-lock.json`);
}
console.log(`${channel}: DuckHTS ${manifest.duckhts} -> npm ${expected}`);
