import random

from simfastq.seed import generate_seed

__author__ = "Davide Rambaldi"
__copyright__ = "Davide Rambaldi"
__license__ = "MIT"


def test_set_seed():
    """Test set_seed"""
    assert generate_seed("AAAA") == "AAAA"
    assert generate_seed("1234") == 1234
    random.seed(1)
    assert generate_seed(None) == 3280387012
