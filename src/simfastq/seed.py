import logging
import random

__author__ = "Davide Rambaldi"
__copyright__ = "Davide Rambaldi"
__license__ = "MIT"

_logger = logging.getLogger(__name__)


def generate_seed(arg_seed):
    """
    Setup random seed or use the one provided from command line.
    If is a number, cast to int
    """
    if arg_seed is None or arg_seed == "":
        _logger.info("Setup random seed")
        return random.randint(0, 2**32 - 1)
    else:
        _logger.info(f"Using seed {arg_seed}")
        try:
            return int(arg_seed)
        except ValueError:
            return arg_seed
