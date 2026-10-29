"""JSON settings validated against the packaged schema, shared with VS Code.

Settings contain ordinary JSON values. Units, descriptions and allowed values belong to
``settings.schema.json``, not to each configuration. Archived settings are self-contained JSON.
"""

from __future__ import annotations

import copy
import json
from functools import lru_cache
from importlib import resources
from pathlib import Path
from typing import Any

import numpy as np
from jsonschema import Draft202012Validator
from jsonschema import ValidationError


@lru_cache(maxsize=1)
def settings_validator() -> Draft202012Validator:
    schema = json.loads((resources.files(__name__) / "settings.schema.json").read_text())
    Draft202012Validator.check_schema(schema)
    return Draft202012Validator(schema)


def validate_settings(values: dict[str, Any], name: str) -> None:
    """Reject unknown settings and invalid values before building or running the solver."""
    try:
        settings_validator().validate(values)
    except ValidationError as error:
        path = "/".join(str(part) for part in error.absolute_path) or "<root>"
        raise ValueError(f"invalid setting '{path}' in '{name}': {error.message}") from error
    json.dumps(values, allow_nan=False)


class Settings(dict[str, Any]):
    """Nested JSON values and their source name; serialisation reflects current effective values."""

    def __init__(self, values: dict[str, Any], name: str = "<string>") -> None:
        validate_settings(values, name)
        super().__init__(copy.deepcopy(values))
        self.name = name

    @property
    def json_text(self) -> str:
        """Self-contained, validated effective settings for storage and cache keys."""
        validate_settings(self, self.name)
        return json.dumps(self, indent=4, allow_nan=False) + "\n"

    def get_path(self, path: str) -> Any:
        """Look up a setting by a slash-separated path, e.g. ``replay/time_step``."""
        node: Any = self
        for key in path.split("/"):
            if not isinstance(node, dict) or key not in node:
                raise KeyError(f"setting '{path}' not found in '{self.name}'")
            node = node[key]
        return node

    def copy_with(self, overrides: dict[str, Any]) -> Settings:
        """Return an independently validated copy with slash-separated paths replaced."""
        values = copy.deepcopy(dict(self))
        for path, value in overrides.items():
            keys = path.split("/")
            if any(not key or ":" in key for key in keys):
                raise ValueError(f"invalid setting path '{path}'; use slash-separated JSON keys")
            node = values
            for key in keys[:-1]:
                child = node.setdefault(key, {})
                if not isinstance(child, dict):
                    raise ValueError(f"setting '{path}' traverses non-object '{key}'")
                node = child
            node[keys[-1]] = value
        return Settings(values, self.name)


def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    values: dict[str, Any] = {}
    for key, value in pairs:
        if key in values:
            raise ValueError(f"duplicate JSON setting '{key}'")
        values[key] = value
    return values


def reject_constant(value: str) -> None:
    raise ValueError(f"nonfinite JSON number '{value}' is not a valid setting")


def parse_settings(json_text: str, name: str = "<string>") -> Settings:
    """Parse strict JSON; reject duplicate keys, nonfinite numbers and schema violations."""
    values = json.loads(json_text, object_pairs_hook=unique_object, parse_constant=reject_constant)
    if not isinstance(values, dict):
        raise ValueError(f"settings in '{name}' must be a JSON object")
    return Settings(values, name)


def parse_overrides(settings: Settings, assignments: list[str]) -> dict[str, Any]:
    """Parse ``path=value`` assignments, using JSON literals or unquoted string values.

    Quote arrays as one shell argument: ``--set 'writers/names=["local"]'``.
    Schema validation catches unknown keys and wrong types, including newly added settings.
    """
    overrides: dict[str, Any] = {}
    for assignment in assignments:
        path, separator, text = assignment.partition("=")
        if not separator:
            raise ValueError(f"override {assignment!r} is not of the form path=value")
        try:
            value = json.loads(text, object_pairs_hook=unique_object, parse_constant=reject_constant)
        except json.JSONDecodeError:
            value = text
        overrides[path] = value
    settings.copy_with(overrides)
    return overrides


def replay_times(settings: Settings) -> np.ndarray:
    """Replay times from ``time_start`` through ``time_end`` in steps of ``time_step``."""
    replay = settings["replay"]
    start = float(replay["time_start"])
    end = float(replay["time_end"])
    step = float(replay["time_step"])
    if step <= 0.0 or end < start:
        raise ValueError("replay needs time_step > 0 and time_end >= time_start")
    n_time = int(np.floor((end - start) / step + 1e-9)) + 1
    return start + step * np.arange(n_time)


def available_settings() -> list[str]:
    """Names of the configurations shipped with the package."""
    root = resources.files(__name__)
    return sorted(entry.name for entry in root.iterdir() if entry.is_dir() and (entry / "rtgsfit.json").is_file())


def load_settings(name_or_path: str | Path) -> Settings:
    """Load a packaged name (e.g. ``analytic``) or a path to a JSON settings file."""
    path = Path(name_or_path).expanduser()
    if path.suffix == ".json":
        return parse_settings(path.read_text(), name=str(path))
    if path.suffix == ".xml":
        raise ValueError("XML settings are no longer supported; use an rtgsfit.json configuration")
    resource = resources.files(__name__) / str(name_or_path) / "rtgsfit.json"
    if not resource.is_file():
        raise FileNotFoundError(f"no settings '{name_or_path}'; available: {', '.join(available_settings())}")
    return parse_settings(resource.read_text(), name=str(name_or_path))
