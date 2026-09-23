// Same-origin static server for the browser test.  Everything the page needs,
// the duckdb-wasm runtimes included, comes from this origin, which is the
// deployment shape the package exists for.  Range requests are honoured
// because DuckHTS reads remote files through htslib with HTTP range reads.
import { createReadStream } from "node:fs";
import { stat } from "node:fs/promises";
import { createServer } from "node:http";
import path from "node:path";
import { fileURLToPath } from "node:url";

const packageRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "../..");
const mounts = {
  "/node_modules/": path.join(packageRoot, "node_modules"),
  "/duckhts/": path.join(packageRoot, "dist"),
  "/src/": path.join(packageRoot, "src"),
  "/data/": path.join(packageRoot, "../test/data"),
  "/": path.join(packageRoot, "test/browser/site"),
};
const contentTypes = {
  ".html": "text/html",
  ".js": "text/javascript",
  ".mjs": "text/javascript",
  ".wasm": "application/wasm",
};

function resolve(urlPath) {
  for (const [prefix, dir] of Object.entries(mounts)) {
    if (!urlPath.startsWith(prefix)) continue;
    const file = path.join(dir, decodeURIComponent(urlPath.slice(prefix.length)));
    return file.startsWith(dir + path.sep) ? file : null;
  }
  return null;
}

const port = Number(process.env.PORT || 8765);

createServer(async (request, response) => {
  const file = resolve(new URL(request.url, "http://localhost").pathname);
  const info = file && (await stat(file).catch(() => null));
  if (!info?.isFile()) {
    response.writeHead(404).end();
    return;
  }
  const headers = {
    "Accept-Ranges": "bytes",
    "Content-Type": contentTypes[path.extname(file)] || "application/octet-stream",
  };
  const range = /^bytes=(\d*)-(\d*)$/.exec(request.headers.range || "");
  if (!range) {
    response.writeHead(200, { ...headers, "Content-Length": info.size });
    if (request.method === "HEAD") return response.end();
    createReadStream(file).pipe(response);
    return;
  }
  const start = range[1] === "" ? Math.max(0, info.size - Number(range[2])) : Number(range[1]);
  const end = range[1] !== "" && range[2] !== "" ? Math.min(Number(range[2]), info.size - 1) : info.size - 1;
  if (start > end) {
    response.writeHead(416, { "Content-Range": `bytes */${info.size}` }).end();
    return;
  }
  response.writeHead(206, {
    ...headers,
    "Content-Length": end - start + 1,
    "Content-Range": `bytes ${start}-${end}/${info.size}`,
  });
  if (request.method === "HEAD") return response.end();
  createReadStream(file, { start, end }).pipe(response);
}).listen(port, "127.0.0.1");
