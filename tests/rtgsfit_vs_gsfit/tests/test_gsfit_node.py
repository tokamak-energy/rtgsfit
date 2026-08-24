import pytest

from rtgsfit_vs_gsfit.gsfit import gsfit_node


def test_gsfit_node_builds_expected_path():
    assert gsfit_node({"run_name": "RUN_A"}, "CONSTRAINTS.TEST") == "\\GSFIT::TOP.RUN_A.CONSTRAINTS.TEST"


def test_gsfit_node_requires_run_name():
    with pytest.raises(KeyError, match="cfg must contain 'run_name'"):
        gsfit_node({}, "CONSTRAINTS.TEST")
