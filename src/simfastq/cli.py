"""
CLI for simfastq the fastq sequences generator

``[options.entry_points]`` section in ``setup.cfg``::

    console_scripts =
         fibonacci = simfastq.skeleton:run

Then run ``pip install .`` (or ``pip install -e .`` for editable mode)
which will install the command ``simfastq`` inside your current environment.

References:
    - https://setuptools.pypa.io/en/latest/userguide/entry_point.html
    - https://pip.pypa.io/en/stable/reference/pip_install
"""

import argparse
import logging
import random
import sys

from rich import pretty, print

from simfastq import __version__
from simfastq.generator import random_record
from simfastq.seed import generate_seed
from simfastq.writer import write_record

__author__ = "Davide Rambaldi"
__copyright__ = "Davide Rambaldi"
__license__ = "MIT"

_logger = logging.getLogger(__name__)


# ---- Python API ----
# The functions defined in this section can be imported by users in their
# Python scripts/interactive interpreter, e.g. via
# `from simfastq.skeleton import fib`,
# when using this Python module as a library.


def simfastq(arg_seed):
    """simfastq entry point"""
    print("[bold green]SIMFASTQ[/bold green]")
    print("[green]Simulation started[/green]")
    seed = generate_seed(arg_seed)
    print(f"Using seed: {seed}")
    random.seed(seed)
    test = random_record(sequence_length=100, x_pos=1, y_pos=1, read=1)
    print("[bold green]READ:[/bold green]")
    write_record(test, sys.stdout)


# ---- CLI ----
# The functions defined in this section are wrappers around the main Python
# API allowing them to be called directly from the terminal as a CLI
# executable/script.


def parse_args(args):
    """Parse command line parameters

    Args:
      args (List[str]): command line parameters as list of strings
          (for example  ``["--help"]``).

    Returns:
      :obj:`argparse.Namespace`: command line parameters namespace
    """
    parser = argparse.ArgumentParser(description="Generate artificial fastq sequences.")

    parser.add_argument(
        "--version",
        action="version",
        version=f"simfastq {__version__}",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        dest="loglevel",
        help="set loglevel to INFO",
        action="store_const",
        const=logging.INFO,
    )

    parser.add_argument(
        "-vv",
        "--very-verbose",
        dest="loglevel",
        help="set loglevel to DEBUG",
        action="store_const",
        const=logging.DEBUG,
    )

    parser.add_argument(
        "-n",
        "--number-of-reads",
        dest="number_of_reads",
        help="number of reads per file -> example: -n 10000",
        default=1,
    )

    parser.add_argument(
        "-p",
        "--paired-end",
        dest="paired_end",
        help="generate paired end reads.",
        action="store_true",
    )

    parser.add_argument("-s", "--seed", dest="seed", help="generation random seed")

    parser.add_argument(
        "-S", "--samples", dest="samples", help="number of samples", default=1
    )

    parser.add_argument("outdir", metavar="OUTDIR", help="output files directory.")

    return parser.parse_args(args)


def setup_logging(loglevel):
    """Setup basic logging

    Args:
      loglevel (int): minimum loglevel for emitting messages
    """
    logformat = "[%(asctime)s] %(levelname)s:%(name)s:%(message)s"
    logging.basicConfig(
        level=loglevel, stream=sys.stdout, format=logformat, datefmt="%Y-%m-%d %H:%M:%S"
    )


def main(args):
    """Wrapper allowing :func:`fib` to be called with string arguments in a CLI fashion

    Instead of returning the value from :func:`fib`, it prints the result to the
    ``stdout`` in a nicely formatted message.

    Args:
      args (List[str]): command line parameters as list of strings
          (for example  ``["--verbose", "42"]``).
    """
    pretty.install()

    args = parse_args(args)
    setup_logging(args.loglevel)
    _logger.debug("Starting crazy calculations...")

    print("[bold green]simfastq[/bold green]")
    print("Args:", args)
    simfastq(arg_seed=args.seed)
    _logger.info("Script ends here")


def run():
    """Calls :func:`main` passing the CLI arguments extracted from :obj:`sys.argv`

    This function can be used as entry point to create console scripts with setuptools.
    """
    main(sys.argv[1:])


if __name__ == "__main__":
    # ^  This is a guard statement that will prevent the following code from
    #    being executed in the case someone imports this file instead of
    #    executing it as a script.
    #    https://docs.python.org/3/library/__main__.html

    # After installing your project with pip, users can also run your Python
    # modules as scripts via the ``-m`` flag, as defined in PEP 338::
    #
    #     python -m simfastq.skeleton 42
    #
    run()
