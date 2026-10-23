"""The C core's own tests (`make -C c test`): finite differences, nulls, boundary, Poisson, LAPACK."""

from __future__ import annotations

import shutil
import subprocess

import pytest
from rtgsfit.library import c_directory


@pytest.mark.skipif(shutil.which("make") is None, reason="needs make")
def test_c_unit_tests() -> None:
    completed = subprocess.run(["make", "-C", str(c_directory()), "test"], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout + completed.stderr
