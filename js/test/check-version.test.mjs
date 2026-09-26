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

async function checkout({ source, npmVersion, devPin, signedPin }) {
  const root = await mkdtemp(path.join(tmpdir(), "duckhts-version-"));
  const js = path.join(root, "js");
  await mkdir(path.join(js, "scripts"), { recursive: true });
  await copyFile(script, path.join(js, "scripts", "check-version.mjs"));
  await writeFile(path.join(root, "description.yml"), `version: ${source}\n`);
  await writeFile(path.join(js, "package.json"), JSON.stringify({ version: npmVersion }));
  await writeFile(path.join(js, "package-lock.json"), JSON.stringify({
    version: npmVersion, packages: { "": { version: npmVersion } },
  }));
  await writeFile(path.join(js, "artifacts-dev.json"), JSON.stringify({ duckhts: devPin }));
  await writeFile(path.join(js, "artifacts.json"), JSON.stringify({ duckhts: signedPin }));
  return { root, js };
}

async function check(js, channel, mode) {
  const args = [path.join(js, "scripts", "check-version.mjs"), channel];
  if (mode) args.push(mode);
  return run(process.execPath, args, {
    env: { ...process.env, DUCKHTS_NPM_CHANNEL: channel },
  }).then(
    (result) => ({ ...result, code: 0 }),
    (error) => error,
  );
}

test("development checkout checks dev identity and skips signed identity only in PR mode", async () => {
  const { root, js } = await checkout({
    source: "1.5.2.9002", npmVersion: "1.5.2-9002",
    devPin: "1.5.2.9002", signedPin: "1.5.2",
  });
  try {
    const dev = await check(js, "dev", "pr");
    assert.equal(dev.code, 0, dev.stderr);
    assert.match(dev.stdout, /dev: DuckHTS 1\.5\.2\.9002 -> npm 1\.5\.2-9002/);
    assert.equal((await check(js, "dev", "publish")).code, 0);

    const signed = await check(js, "signed", "pr");
    assert.equal(signed.code, 0, signed.stderr);
    assert.match(signed.stdout, /signed: skipping identity check \(description\.yml 1\.5\.2\.9002 is dev\)/);

    const wrongPublish = await check(js, "signed", "publish");
    assert.equal(wrongPublish.code, 1);
    assert.match(wrongPublish.stderr, /signed publish requires a signed version/);

    await writeFile(path.join(js, "artifacts-dev.json"), JSON.stringify({ duckhts: "1.5.2.9001" }));
    const stalePin = await check(js, "dev", "pr");
    assert.equal(stalePin.code, 1);
    assert.match(stalePin.stderr, /dev manifest 1\.5\.2\.9001 disagrees with description\.yml 1\.5\.2\.9002/);

    await writeFile(path.join(js, "artifacts-dev.json"), JSON.stringify({ duckhts: "1.5.2.9002" }));
    await writeFile(path.join(js, "package.json"), JSON.stringify({ version: "1.5.2-9001" }));
    const stalePackage = await check(js, "dev", "publish");
    assert.equal(stalePackage.code, 1);
    assert.match(stalePackage.stderr, /dev npm version must be 1\.5\.2-9002/);
  } finally {
    await rm(root, { recursive: true, force: true });
  }
});

test("release checkout checks signed identity and skips dev identity only in PR mode", async () => {
  const { root, js } = await checkout({
    source: "1.5.3", npmVersion: "1.5.3",
    devPin: "1.5.2.9002", signedPin: "1.5.3",
  });
  try {
    const signed = await check(js, "signed", "pr");
    assert.equal(signed.code, 0, signed.stderr);
    assert.match(signed.stdout, /signed: DuckHTS 1\.5\.3 -> npm 1\.5\.3/);
    assert.equal((await check(js, "signed", "publish")).code, 0);

    const dev = await check(js, "dev", "pr");
    assert.equal(dev.code, 0, dev.stderr);
    assert.match(dev.stdout, /dev: skipping identity check \(description\.yml 1\.5\.3 is signed\)/);

    const wrongPublish = await check(js, "dev");
    assert.equal(wrongPublish.code, 1);
    assert.match(wrongPublish.stderr, /dev publish requires a dev version/);

    await writeFile(path.join(js, "artifacts.json"), JSON.stringify({ duckhts: "1.5.2" }));
    const stalePin = await check(js, "signed", "publish");
    assert.equal(stalePin.code, 1);
    assert.match(stalePin.stderr, /signed manifest 1\.5\.2 disagrees with description\.yml 1\.5\.3/);

    await writeFile(path.join(js, "artifacts.json"), JSON.stringify({ duckhts: "1.5.3" }));
    await writeFile(path.join(js, "package.json"), JSON.stringify({ version: "1.5.2" }));
    const stalePackage = await check(js, "signed", "pr");
    assert.equal(stalePackage.code, 1);
    assert.match(stalePackage.stderr, /signed npm version must be 1\.5\.3/);

    await writeFile(path.join(js, "package.json"), JSON.stringify({ version: "1.5.3" }));
    await writeFile(path.join(js, "package-lock.json"), JSON.stringify({
      version: "1.5.3", packages: { "": { version: "1.5.2" } },
    }));
    const staleLock = await check(js, "signed", "publish");
    assert.equal(staleLock.code, 1);
    assert.match(staleLock.stderr, /signed npm version must be 1\.5\.3/);

    await writeFile(path.join(root, "description.yml"), "version: 1.5.3.90020\n");
    const malformed = await check(js, "dev", "pr");
    assert.equal(malformed.code, 1);
    assert.match(malformed.stderr, /description\.yml must declare X\.Y\.Z or X\.Y\.Z\.9NNN/);
  } finally {
    await rm(root, { recursive: true, force: true });
  }
});
