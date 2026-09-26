// Check the npm release identity against its pinned extension version.
import { readFile } from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";

const root = path.dirname(path.dirname(fileURLToPath(import.meta.url)));
const channel = process.argv[2] ?? process.env.DUCKHTS_NPM_CHANNEL;
const mode = process.argv[3] ?? "publish";
if (channel !== "signed" && channel !== "dev") throw new Error("Select signed or dev");
if (process.argv[2] && process.env.DUCKHTS_NPM_CHANNEL && process.argv[2] !== process.env.DUCKHTS_NPM_CHANNEL) {
  throw new Error("Channel argument and DUCKHTS_NPM_CHANNEL disagree");
}
if (mode !== "pr" && mode !== "publish") throw new Error("Select pr or publish mode");

const description = await readFile(path.join(root, "../description.yml"), "utf8");
const version = /^version:[ \t]*(\d+\.\d+\.\d+)(?:\.(9\d{3}))?[ \t]*$/m.exec(description);
if (!version) throw new Error("description.yml must declare X.Y.Z or X.Y.Z.9NNN");
const sourceVersion = version[0].slice("version:".length).trim();
const sourceChannel = version[2] ? "dev" : "signed";
if (channel !== sourceChannel) {
  if (mode === "pr") {
    console.log(`${channel}: skipping identity check (description.yml ${sourceVersion} is ${sourceChannel})`);
    process.exit(0);
  }
  throw new Error(`${channel} publish requires a ${channel} version in description.yml; found ${sourceVersion}`);
}

const packageJson = JSON.parse(await readFile(path.join(root, "package.json"), "utf8"));
const lock = JSON.parse(await readFile(path.join(root, "package-lock.json"), "utf8"));
const manifest = JSON.parse(await readFile(path.join(root,
  channel === "dev" ? "artifacts-dev.json" : "artifacts.json"), "utf8"));
const expected = sourceChannel === "dev" ? `${version[1]}-${version[2]}` : sourceVersion;
if (manifest.duckhts !== sourceVersion) {
  throw new Error(`${channel} manifest ${manifest.duckhts} disagrees with description.yml ${sourceVersion}`);
}
if (packageJson.version !== expected || lock.version !== expected || lock.packages[""].version !== expected) {
  throw new Error(`${channel} npm version must be ${expected} in package.json and package-lock.json`);
}
console.log(`${channel}: DuckHTS ${sourceVersion} -> npm ${expected}`);
