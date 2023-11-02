import random

from Bio.Seq import Seq

from simfastq.generator import (
    random_phred_quality,
    random_record,
    random_sequence,
    sequence_description,
    sequence_name,
)

__author__ = "Davide Rambaldi"
__copyright__ = "Davide Rambaldi"
__license__ = "MIT"

expected_seq = Seq("GATACCCCGA")
expected_qual = {"phred_quality": [17, 72, 8, 32, 15, 63, 57, 60, 83, 48]}
expected_qual_record = {"phred_quality": [62, 3, 49, 55, 77, 0, 89, 57, 34, 92]}


expected_record_name = "SIMFASTQ:1:SIMFASTQ_FLOWCELL_ID:1:1:1:1"
expected_record_description = "SIMFASTQ:1:SIMFASTQ_FLOWCELL_ID:1:1:1:1 1:N:0:1"
expected_description = "1:N:0:1"
expected_description_index = "1:N:0:CGACGATA"


def test_random_sequence():
    """Test random sequence generator"""
    random.seed(1)
    assert random_sequence(10) == expected_seq


def test_random_phred_quality():
    """Test random phred quality"""
    random.seed(1)
    assert random_phred_quality(10) == expected_qual


def test_random_record():
    """Test generate random record"""
    random.seed(1)
    r = random_record(sequence_length=10, x_pos=1, y_pos=1, read=1)
    assert r.seq == expected_seq
    assert r.letter_annotations == expected_qual_record
    assert r.description == expected_record_description
    assert r.name == expected_record_name
    assert r.id == expected_record_name


def test_sequence_description():
    """Set seq record description generation"""
    assert sequence_description(1, 1) == expected_description
    assert sequence_description(1, "CGACGATA") == expected_description_index


def test_sequence_name():
    """Test sequence name"""
    assert sequence_name(1, 1) == expected_record_name
    assert (
        sequence_name(1, 1, UMI="ACGTA")
        == "SIMFASTQ:1:SIMFASTQ_FLOWCELL_ID:1:1:1:1:ACGTA"
    )
