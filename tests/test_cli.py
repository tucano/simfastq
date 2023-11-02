import pytest

from simfastq.cli import main, simfastq

__author__ = "Davide Rambaldi"
__copyright__ = "Davide Rambaldi"
__license__ = "MIT"


@pytest.mark.skip(reason="Need to define a UI first")
def test_simfastq():
    """API Tests"""
    simfastq("PIPPO")


def test_main(capsys, tmpdir):
    """CLI Tests"""
    # capsys is a pytest fixture that allows asserts against stdout/stderr
    # https://docs.pytest.org/en/stable/capture.html
    main([str(tmpdir)])
    captured = capsys.readouterr()
    assert "Args: Namespace" in captured.out
