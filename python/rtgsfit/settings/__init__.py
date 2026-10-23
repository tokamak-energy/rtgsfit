"""XML settings, in the style of the IMAS Data Dictionary `code.parameters` field.

A settings file is a tree of elements under a `<parameters>` root. Every leaf element carries a
`type` attribute that says how to parse its text:

```xml
<parameters>
  <grid>
    <n_r type="int">65</n_r>
    <r_min type="float" unit="m">9.25</r_min>
  </grid>
  <fit>
    <rogowski_names type="str_array">INIVC000 BVLT</rogowski_names>
  </fit>
</parameters>
```

Supported types: `int`, `float`, `bool`, `str`, `int_array`, `float_array`, `str_array`.
Arrays are whitespace separated. Optional `unit` and `description` attributes are kept only for
documentation. Settings are read into nested dictionaries, and the raw XML text is kept so that it
can be stored in MDSplus next to the results (`INPUT:SETTINGS`).
"""

from __future__ import annotations

import copy
import xml.etree.ElementTree as ElementTree
from importlib import resources
from pathlib import Path
from typing import Any

import numpy as np

SETTINGS_ROOT_TAG = "parameters"

_SCALAR_PARSERS = {
    "int": int,
    "float": float,
    "str": lambda text: text,
    "bool": lambda text: _parse_bool(text),
}
_ARRAY_DTYPES = {"int_array": np.int64, "float_array": np.float64}


def _parse_bool(text: str) -> bool:
    lowered = text.strip().lower()
    if lowered in ("true", "1", "yes"):
        return True
    if lowered in ("false", "0", "no"):
        return False
    raise ValueError(f"cannot parse {text!r} as bool")


def _parse_leaf(element: ElementTree.Element, path: str) -> Any:
    value_type = element.get("type")
    text = (element.text or "").strip()
    if value_type is None:
        raise ValueError(f"settings leaf '{path}' has no 'type' attribute")
    if value_type in _SCALAR_PARSERS:
        return _SCALAR_PARSERS[value_type](text)
    if value_type in _ARRAY_DTYPES:
        return np.array(text.split(), dtype=_ARRAY_DTYPES[value_type])
    if value_type == "str_array":
        return text.split()
    raise ValueError(f"settings leaf '{path}' has unknown type {value_type!r}")


def _element_to_dict(element: ElementTree.Element, path: str) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for child in element:
        child_path = f"{path}/{child.tag}"
        if child.tag in result:
            raise ValueError(f"duplicate settings element '{child_path}'")
        if len(child) > 0:
            result[child.tag] = _element_to_dict(child, child_path)
        else:
            result[child.tag] = _parse_leaf(child, child_path)
    return result


class Settings(dict[str, Any]):
    """Nested dictionary of settings, plus the XML text it was read from."""

    def __init__(self, values: dict[str, Any], xml_text: str, name: str) -> None:
        super().__init__(values)
        self.xml_text = xml_text
        self.name = name

    def get_path(self, path: str) -> Any:
        """Look up a value by a '/' separated path, e.g. `settings.get_path("grid/n_r")`."""
        node: Any = self
        for key in path.split("/"):
            if not isinstance(node, dict) or key not in node:
                raise KeyError(f"setting '{path}' not found in '{self.name}'")
            node = node[key]
        return node

    def copy_with(self, overrides: dict[str, Any]) -> Settings:
        """Return a copy with values replaced, e.g. `{"replay/pulse": 13343}`.

        The stored XML text is updated too, so the values written to MDSplus match the values used.
        """
        values = copy.deepcopy(dict(self))
        root = ElementTree.fromstring(self.xml_text)
        for path, value in overrides.items():
            keys = path.split("/")
            node = values
            for key in keys[:-1]:
                node = node.setdefault(key, {})
            node[keys[-1]] = value
            element = root
            for key in keys:
                child = element.find(key)
                if child is None:
                    child = ElementTree.SubElement(element, key)
                element = child
            element.set("type", element.get("type") or _infer_type(value))
            element.text = _format_value(value)
        return Settings(values, ElementTree.tostring(root, encoding="unicode"), self.name)


def _infer_type(value: Any) -> str:
    if isinstance(value, bool):
        return "bool"
    if isinstance(value, (int, np.integer)):
        return "int"
    if isinstance(value, (float, np.floating)):
        return "float"
    if isinstance(value, str):
        return "str"
    array = np.asarray(value)
    if array.dtype.kind in "iu":
        return "int_array"
    if array.dtype.kind == "f":
        return "float_array"
    return "str_array"


def _format_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (str, int, float, np.integer, np.floating)):
        return repr(float(value)) if isinstance(value, (float, np.floating)) else str(value)
    return " ".join(repr(float(v)) if isinstance(v, (float, np.floating)) else str(v) for v in value)


def parse_settings(xml_text: str, name: str = "<string>") -> Settings:
    """Parse settings from XML text."""
    root = ElementTree.fromstring(xml_text)
    if root.tag != SETTINGS_ROOT_TAG:
        raise ValueError(f"settings root element must be <{SETTINGS_ROOT_TAG}>, not <{root.tag}>")
    return Settings(_element_to_dict(root, ""), xml_text, name)


def parse_overrides(settings: Settings, assignments: list[str]) -> dict[str, Any]:
    """`["replay/pulse=13343"]` -> `{"replay/pulse": 13343}`, parsed with the type declared in the XML.

    New settings (not in the XML) need an explicit type: `path:type=value`, e.g. `build/debug:bool=true`.
    """
    root = ElementTree.fromstring(settings.xml_text)
    overrides: dict[str, Any] = {}
    for assignment in assignments:
        target, separator, text = assignment.partition("=")
        if not separator:
            raise ValueError(f"override {assignment!r} is not of the form path=value")
        path, _, value_type = target.partition(":")
        element = root.find(path)
        if not value_type:
            if element is None or element.get("type") is None:
                raise ValueError(f"setting '{path}' is not in '{settings.name}'; give its type as {path}:<type>=...")
            value_type = str(element.get("type"))
        leaf = ElementTree.Element("leaf", type=value_type)
        leaf.text = text
        overrides[path] = _parse_leaf(leaf, path)
    return overrides


def replay_times(settings: Settings) -> np.ndarray:
    """The time slices of a replay: `replay/time_start` to `replay/time_end` in steps of `replay/time_step`."""
    replay = settings["replay"]
    start, end, step = float(replay["time_start"]), float(replay["time_end"]), float(replay["time_step"])
    if step <= 0.0 or end < start:
        raise ValueError("replay needs time_step > 0 and time_end >= time_start")
    n_time = int(np.floor((end - start) / step + 1e-9)) + 1
    return start + step * np.arange(n_time)


def available_settings() -> list[str]:
    """Names of the settings shipped with the package."""
    root = resources.files(__name__)
    return sorted(entry.name for entry in root.iterdir() if entry.is_dir() and (entry / "rtgsfit.xml").is_file())


def load_settings(name_or_path: str | Path) -> Settings:
    """Load settings by package name (e.g. `"analytic"`) or by path to an XML file."""
    path = Path(name_or_path)
    if path.suffix == ".xml" and path.is_file():
        return parse_settings(path.read_text(), name=str(path))
    resource = resources.files(__name__) / str(name_or_path) / "rtgsfit.xml"
    if not resource.is_file():
        raise FileNotFoundError(f"no settings '{name_or_path}'; available: {', '.join(available_settings())}")
    return parse_settings(resource.read_text(), name=str(name_or_path))
