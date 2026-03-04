import argparse
import runpy
import sys

from . import __version__


def _normalize_module_name(module: str) -> str:
    if module.startswith("jcvi."):
        return module
    return f"jcvi.{module}"


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="jcvi",
        description="Run jcvi modules, e.g. `jcvi formats.fasta extract`.",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    parser.add_argument(
        "module",
        nargs="?",
        help="Module path without the `jcvi.` prefix (for example: formats.fasta).",
    )
    parser.add_argument("module_args", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)

    if not args.module:
        parser.print_help()
        return 0

    module_name = _normalize_module_name(args.module)

    original_argv = sys.argv[:]
    try:
        # Match `python -m jcvi.<module> ...` behavior for argument parsing.
        sys.argv = [module_name, *args.module_args]
        runpy.run_module(module_name, run_name="__main__", alter_sys=True)
    finally:
        sys.argv = original_argv

    return 0


if __name__ == "__main__":
    main()
