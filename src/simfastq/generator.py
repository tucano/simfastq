import random

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

__author__ = "Davide Rambaldi"
__copyright__ = "Davide Rambaldi"
__license__ = "MIT"


MAX_PHRED_QUALITY = 93


def random_sequence(sequence_length):
    """
    Generate random sequences
    """
    return Seq("".join([random.choice("AGTC") for x in range(sequence_length)]))


def random_phred_quality(sequence_length):
    """
    Generate random phred quality
    """
    return {
        "phred_quality": [
            random.randrange(MAX_PHRED_QUALITY) for x in range(sequence_length)
        ]
    }


def random_record(
    sequence_length,
    x_pos,
    y_pos,
    read,
    UMI=None,
    sample_index=1,
    instrument="SIMFASTQ",
    run_number=1,
    flowcell_ID="SIMFASTQ_FLOWCELL_ID",
    lane=1,
    tile=1,
    is_filtered="N",
    control_number=0,
):
    """
    Generate random seq record
    """
    seq = random_sequence(sequence_length)
    qual = random_phred_quality(sequence_length)
    name = sequence_name(
        x_pos=x_pos,
        y_pos=y_pos,
        UMI=UMI,
        instrument=instrument,
        run_number=run_number,
        flowcell_ID=flowcell_ID,
        lane=lane,
        tile=tile,
    )

    desc = sequence_description(
        read=read,
        sample_index=sample_index,
        is_filtered=is_filtered,
        control_number=control_number,
    )

    desc = f"{name} {desc}"

    record = SeqRecord(
        seq, letter_annotations=qual, description=desc, id=name, name=name
    )
    return record


def sequence_name(
    x_pos,
    y_pos,
    UMI=None,
    instrument="SIMFASTQ",
    run_number=1,
    flowcell_ID="SIMFASTQ_FLOWCELL_ID",
    lane=1,
    tile=1,
):
    """Generate a sequence id and name"""
    seq_id = f"{instrument}:{run_number}:{flowcell_ID}:{lane}:{tile}:{x_pos}:{y_pos}"  # noqa: E231, E501
    if UMI is not None and UMI != "":
        seq_id += f":{UMI}"  # noqa: E231
    return seq_id


def sequence_description(read, sample_index=1, is_filtered="N", control_number=0):
    """
    Generate a sequence description
    """
    return f"{read}:{is_filtered}:{control_number}:{sample_index}"  # noqa: E231
