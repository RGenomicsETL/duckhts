"""Shared cache location for DuckHTS Python staging and benchmark scripts."""

from __future__ import annotations

import os
from pathlib import Path


def duckhts_cache_dir() -> Path:
    configured = os.environ.get("DUCKHTS_CACHE_DIR")
    if configured:
        return Path(configured).expanduser()
    xdg = os.environ.get("XDG_CACHE_HOME")
    if xdg:
        return Path(xdg).expanduser() / "duckhts"
    return Path.home() / ".cache" / "duckhts"


def duckhts_cache_subdir(*parts: str) -> Path:
    relative = Path(*parts)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"invalid DuckHTS cache-relative path: {relative}")
    return duckhts_cache_dir() / relative
