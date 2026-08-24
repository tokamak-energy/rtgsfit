"""Shared helpers for building GSFIT node paths."""


def gsfit_node(cfg: dict, suffix: str) -> str:
    """Return the fully qualified GSFIT node path for the current run."""
    if "run_name" not in cfg:
        raise KeyError("cfg must contain 'run_name' to build a GSFIT node path")

    return f"\\GSFIT::TOP.{cfg['run_name']}.{suffix}"
