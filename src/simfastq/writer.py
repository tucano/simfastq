import logging
import sys

from Bio import SeqIO

__author__ = "Davide Rambaldi"
__copyright__ = "Davide Rambaldi"
__license__ = "MIT"

_logger = logging.getLogger(__name__)


def write_record(seq_record):
    """Write a seqRecord to fastq format"""
    SeqIO.write(seq_record, sys.stdout, "fastq")
