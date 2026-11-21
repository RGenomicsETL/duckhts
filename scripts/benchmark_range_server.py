"""HTTP range server with a request log, for remote-read benchmarks.

Serves the regular files of one directory on 127.0.0.1 and streams each body
in 64 KiB pieces. It prints "ready <port>" once it listens. Every request
appends one tab-separated line to the log:

    method  path  Range header or "-"  status  body bytes written

"Body bytes written" counts what the server handed to the socket before the
response ended or the client closed the connection. A client that abandons an
open-ended range has usually received less: the kernel still held the rest.

usage: benchmark_range_server.py <directory> <log-file> <latency-ms>
"""

import os
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

PIECE_BYTES = 64 * 1024


def parse_range(header, size):
    """(first, last) of a 'bytes=first-[last]' header, or None when unsatisfiable."""
    first_text, last_text = header.split("=", 1)[1].split("-", 1)
    first = int(first_text)
    last = min(int(last_text), size - 1) if last_text else size - 1
    return (first, last) if first <= last else None


class RangeHandler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"

    def log_message(self, *arguments):
        pass

    def record(self, status, written):
        fields = (self.command, self.path, self.headers.get("Range", "-"), str(status), str(written))
        with self.server.log_lock, open(self.server.log_path, "a") as log:
            log.write("\t".join(fields) + "\n")

    def send_status(self, status, length, content_range=None):
        self.send_response(status)
        self.send_header("Accept-Ranges", "bytes")
        if content_range:
            self.send_header("Content-Range", content_range)
        self.send_header("Content-Length", str(length))
        self.end_headers()

    def serve(self):
        time.sleep(self.server.latency_seconds)
        name = self.path.split("?", 1)[0].lstrip("/")
        target = os.path.join(self.server.directory, name)
        if "/" in name or not os.path.isfile(target):
            self.send_status(404, 0)
            self.record(404, 0)
            return
        size = os.path.getsize(target)
        status, span, content_range = 200, (0, size - 1), None
        requested = self.headers.get("Range")
        if requested:
            span = parse_range(requested, size)
            if span is None:
                self.send_status(416, 0, f"bytes */{size}")
                self.record(416, 0)
                return
            status, content_range = 206, f"bytes {span[0]}-{span[1]}/{size}"
        remaining = span[1] - span[0] + 1
        self.send_status(status, remaining, content_range)
        written = 0
        if self.command == "GET":
            try:
                with open(target, "rb") as body:
                    body.seek(span[0])
                    while remaining > 0:
                        piece = body.read(min(PIECE_BYTES, remaining))
                        if not piece:
                            break
                        self.wfile.write(piece)
                        written += len(piece)
                        remaining -= len(piece)
            except (BrokenPipeError, ConnectionResetError):
                self.close_connection = True
        self.record(status, written)

    do_GET = serve
    do_HEAD = serve


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    server = ThreadingHTTPServer(("127.0.0.1", 0), RangeHandler)
    server.directory = sys.argv[1]
    server.log_path = sys.argv[2]
    server.latency_seconds = float(sys.argv[3]) / 1000.0
    server.log_lock = threading.Lock()
    print("ready", server.server_address[1], flush=True)
    server.serve_forever()


if __name__ == "__main__":
    main()
