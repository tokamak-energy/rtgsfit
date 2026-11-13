"""RT-GSFit: real-time equilibrium reconstruction for tokamaks."""

from typing import Any

from rtgsfit.settings import load_settings

__version__ = "0.2.0"

__all__ = ["RtgsfitReplay", "__version__", "load_settings"]


def __getattr__(name: str) -> Any:
    # Imported on first use, so that `import rtgsfit` does not load the replay machinery.
    if name == "RtgsfitReplay":
        from rtgsfit.ppac import RtgsfitReplay

        return RtgsfitReplay
    raise AttributeError(f"module 'rtgsfit' has no attribute {name!r}")
