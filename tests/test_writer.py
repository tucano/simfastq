import sys

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from simfastq.writer import write_record

__author__ = "Davide Rambaldi"
__copyright__ = "Davide Rambaldi"
__license__ = "MIT"


def test_write_record(capsys):
    """Test write record API"""
    record = SeqRecord(Seq("ACGT"), letter_annotations={"phred_quality": [1, 1, 1, 1]})
    write_record(record, sys.stdout)
    captured = capsys.readouterr()
    assert "ACGT" in captured.out
