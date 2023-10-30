from simfastq.cli import main, simfastq

__author__ = "Davide Rambaldi"
__copyright__ = "Davide Rambaldi"
__license__ = "MIT"


def test_simfastq():
    """API Tests"""
    simfastq()
    assert False


def test_main(capsys):
    """CLI Tests"""
    # capsys is a pytest fixture that allows asserts against stdout/stderr
    # https://docs.pytest.org/en/stable/capture.html
    main([])
    captured = capsys.readouterr()
    assert "Args: Namespace" in captured.out
