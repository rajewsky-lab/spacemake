import re
import sys


def find_config():
    return "config.yaml"


def cmdline():
    from snakemake import main

    sys.argv[1] = f"--config-file={sys.argv[1]}"
    sys.argv.append(f"--rerun-incomplete")
    sys.exit(main())


if __name__ == "__main__":
    cmdline()