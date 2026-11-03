"""Verify that a LOAD registration error releases the macro-registration lock.

Run under timeout 30 with the extension path and the compiled probe path.
"""

import ctypes
from pathlib import Path
import sys

import duckdb

extension = Path(sys.argv[1]).resolve(strict=True)
probe = ctypes.CDLL(str(Path(sys.argv[2]).resolve(strict=True)))
probe.extension_init_oom_check.argtypes = [ctypes.c_char_p]
probe.extension_init_oom_check.restype = ctypes.c_int
with duckdb.connect(config={"allow_unsigned_extensions": "true", "threads": "1"}) as connection:
    connection.execute("LOAD '" + str(extension).replace("'", "''") + "'")
    result = probe.extension_init_oom_check(str(extension).encode())
    if result != 0:
        raise AssertionError(f"Native initialization OOM check failed: {result}")
    assert connection.execute("SELECT 42").fetchone() == (42,)
print("Registration lock reacquired after injected allocation-error exit")
