# Test-only vendored libraries

These libraries are used only by native C property tests under
`test/duckvep/property`.
They are not linked into the DuckDB extension.

- `greatest`: <https://github.com/silentbicycle/greatest>, v1.5.0,
  commit `11a6af1`, ISC license.
- `theft`: <https://github.com/silentbicycle/theft>, v0.4.5,
  commit `62e093d`, ISC license. The vendored copy contains a local
  `theft_random_bits_bulk()` shift-by-64 guard so the test runner itself is
  clean under UBSan; the sole source change is recorded in
  `patches/theft-shift-by-64.patch`.
