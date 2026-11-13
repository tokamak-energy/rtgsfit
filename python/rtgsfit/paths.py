"""Path resolution shared by replay outputs and local writers."""

from __future__ import annotations

from pathlib import Path


def resolve_output_path(path: str | Path, relative_to: Path) -> Path:
    path = Path(path).expanduser()
    return path if path.is_absolute() else relative_to / path
