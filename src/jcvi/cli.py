# cli.py
import argparse

from . import __version__


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    args = parser.parse_args()


if __name__ == "__main__":
    main()
