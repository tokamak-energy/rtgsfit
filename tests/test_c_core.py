"""The C core's own tests (`make -C c test`): finite differences, nulls, boundary, Poisson, LAPACK."""

from __future__ import annotations

import re
import shutil
import subprocess

import pytest
from rtgsfit.library import c_directory


@pytest.mark.skipif(shutil.which("make") is None, reason="needs make")
def test_c_unit_tests() -> None:
    completed = subprocess.run(["make", "-C", str(c_directory()), "test"], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout + completed.stderr


def test_every_c_array_parameter_has_an_adjacent_checked_length() -> None:
    """Include private helpers, not just the exported API."""
    functions = re.compile(r"\b(\w+)\s*\(([^{};]*)\)\s*\{")
    arrays = re.compile(r"\b(\w+)\s*\[\s*\]")
    paths = [
        *c_directory().glob("src/*.[ch]"),
        *c_directory().glob("tests/*.[ch]"),
    ]
    for path in sorted(paths):
        text = path.read_text()
        for function in functions.finditer(text):
            depth = 1
            end = function.end()
            while depth and end < len(text):
                depth += (text[end] == "{") - (text[end] == "}")
                end += 1
            body = text[function.end() : end]
            parameters = function[2]
            for array in arrays.finditer(parameters):
                location = f"{path.name}:{function[1]}:{array[1]}"
                length = re.match(r"\s*,\s*(?:const\s+)?int32_t\s+(\w+)", parameters[array.end() :])
                assert length, f"{location}: array must be followed immediately by its length"
                assert re.search(rf"RTGSFIT_CHECK_(?:SIZE|MIN_SIZE|ARRAY)\([^;]*\b{length[1]}\b", body), (
                    f"{location}: missing length validation"
                )
