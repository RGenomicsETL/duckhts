import { execFileSync } from "node:child_process";

export class UnavailableDevArtifactsError extends Error {}

// Staging calls this only for downloads; verified local dist files work offline.
export function checkDevArtifacts(manifest, names) {
  const { repository, run } = manifest.source;
  const endpoint = `repos/${repository}/actions/runs/${run}/artifacts?per_page=100`;
  const pages = JSON.parse(execFileSync("gh", ["api", endpoint, "--paginate", "--slurp"], {
    encoding: "utf8",
  }));
  if (!Array.isArray(pages) || pages.some((page) => !Array.isArray(page?.artifacts))) {
    throw new Error(`${endpoint}: unexpected GitHub artifacts response`);
  }
  const artifacts = pages.flatMap((page) => page.artifacts);
  const unavailable = names.flatMap((name) => {
    const matches = artifacts.filter((artifact) => artifact.name === name);
    if (!matches.length) return [`${name} (missing)`];
    if (matches.some((artifact) => artifact.expired === true)) return [`${name} (expired)`];
    return [];
  });
  if (unavailable.length) {
    throw new UnavailableDevArtifactsError(
      `Dev pin ${repository} Actions run ${run}: artifact(s) ${unavailable.join(", ")} unavailable. ` +
      "Refresh the dev pin with the next X.Y.Z.9NNN development bump, or build the wasm extension " +
      "from source (README.md: Browser wasm/duckdb-wasm local setup; Makefile: wasm-playwright-test) " +
      "and stage your own manifest.",
    );
  }
  if (names.some((name) => !artifacts.some((artifact) => artifact.name === name && artifact.expired === false))) {
    throw new Error(`${endpoint}: artifact expiry status is unknown`);
  }
}
