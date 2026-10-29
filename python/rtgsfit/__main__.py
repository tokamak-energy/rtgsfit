"""Command line interface: `python -m rtgsfit <command> ...` (or `rtgsfit <command> ...`).

```
rtgsfit replay analytic                                   # replay, then run the writers
rtgsfit replay st40 --set replay/pulse=13343 --set writers/local/path=13343.npz
rtgsfit preshot st40                                      # build a new preshot, then run the writers
rtgsfit write-c-preshot st40 --output-dir c/preshot       # C source of the replay preshot, for the PCS
rtgsfit write-c-preshot --pulse 52014922 --run-name RUN21 # the same, from MDSplus
rtgsfit write-csv path/to/standard_utility/csv_files/RTGSFIT.csv
rtgsfit write-data-dictionary documentation/data_dictionary.md
rtgsfit list-settings
```
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from rtgsfit import data_dictionary as dd
from rtgsfit.settings import Settings
from rtgsfit.settings import available_settings
from rtgsfit.settings import load_settings
from rtgsfit.settings import parse_overrides


def _settings(args: argparse.Namespace) -> Settings:
    settings = load_settings(args.settings)
    if args.set:
        settings = settings.copy_with(parse_overrides(settings, args.set))
    return settings


def _replay(args: argparse.Namespace) -> None:
    from rtgsfit.replay import run_replay

    run_replay(_settings(args), verbose=not args.quiet)


def _preshot(args: argparse.Namespace) -> None:
    from rtgsfit.replay import build_preshot
    from rtgsfit.replay import write_preshot

    settings = _settings(args)
    preshot = build_preshot(settings)
    print(preshot.summary())
    write_preshot(preshot, settings)


def _write_c_preshot(args: argparse.Namespace) -> None:
    from rtgsfit.database_writers.c_source import write_c_source
    from rtgsfit.replay import read_preshot

    settings = _settings(args)
    if args.pulse is not None:
        settings = settings.copy_with({
            "preshot/reader": "rtgsfit_mdsplus",
            "preshot/rtgsfit_mdsplus/pulse": int(args.pulse),
            "preshot/rtgsfit_mdsplus/run_name": str(args.run_name),
        })  # fmt: skip
    preshot = read_preshot(settings)
    for path in write_c_source(preshot, args.output_dir):
        print(f"Wrote {path}")


def _write_csv(args: argparse.Namespace) -> None:
    _write_text(dd.generate_csv(), args.output)


def _write_data_dictionary(args: argparse.Namespace) -> None:
    _write_text(dd.generate_markdown(), args.output)


def _write_text(text: str, output: str | None) -> None:
    if output is None:
        sys.stdout.write(text)
    else:
        Path(output).write_text(text)
        print(f"Wrote {output}")


def _list_settings(_: argparse.Namespace) -> None:
    print("\n".join(available_settings()))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rtgsfit", description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)

    def with_settings(command: argparse.ArgumentParser, default: str | None = None) -> None:
        command.add_argument("settings", nargs="?" if default else None, default=default,
                             help="settings name (see list-settings) or path to a JSON file")  # fmt: skip
        command.add_argument("--set", action="append", metavar="PATH=VALUE",
                             help="override a setting, e.g. replay/pulse=13343; repeatable")  # fmt: skip

    replay = commands.add_parser("replay", help="replay rtgsfit() and write the results")
    with_settings(replay)
    replay.add_argument("--quiet", action="store_true", help="do not print every time slice")
    replay.set_defaults(func=_replay)

    preshot = commands.add_parser("preshot", help="build a preshot with preshot/builder and write it")
    with_settings(preshot)
    preshot.set_defaults(func=_preshot)

    c_preshot = commands.add_parser("write-c-preshot", help="write preshot_sizes.h and preshot_data.c")
    with_settings(c_preshot, default="st40")
    c_preshot.add_argument("--pulse", type=int, help="read the preshot from this RTGSFIT MDSplus pulse")
    c_preshot.add_argument("--run-name", default="BEST", help="with --pulse: the MDSplus run")
    c_preshot.add_argument("--output-dir", default=str(Path(__file__).resolve().parents[2] / "c" / "preshot"))
    c_preshot.set_defaults(func=_write_c_preshot)

    csv = commands.add_parser("write-csv", help="write the standard_utility RTGSFIT.csv")
    csv.add_argument("output", nargs="?", help="output file; standard output if omitted")
    csv.set_defaults(func=_write_csv)

    markdown = commands.add_parser("write-data-dictionary", help="write documentation/data_dictionary.md")
    markdown.add_argument("output", nargs="?", help="output file; standard output if omitted")
    markdown.set_defaults(func=_write_data_dictionary)

    listing = commands.add_parser("list-settings", help="list the settings shipped with RT-GSFit")
    listing.set_defaults(func=_list_settings)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
