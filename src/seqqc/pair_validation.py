
# Import the FASTQ parser we already developed.
from itertools import zip_longest
from seqqc.fastq_parser import read_fastq


def get_spot_id(read_id):
    """
    Extract an SRA spot identifier from a FASTQ header.

    Example:
        "SRR7992936.1 1 length=222"
        becomes "SRR7992936.1"
    """

    # Keep only the part before the first space.
    return read_id.split()[0]


def check_read_pairs(r1_path, r2_path):
    """Summarize shared spot IDs and record-by-record alignment."""
    r1_ids = set()
    r2_ids = set()
    position_matches = 0
    position_mismatches = 0
    extra_r1_records = 0
    extra_r2_records = 0

    for r1_record, r2_record in zip_longest(
        read_fastq(r1_path), read_fastq(r2_path)
    ):
        r1_id = get_spot_id(r1_record[0]) if r1_record is not None else None
        r2_id = get_spot_id(r2_record[0]) if r2_record is not None else None

        if r1_id is not None:
            r1_ids.add(r1_id)
        if r2_id is not None:
            r2_ids.add(r2_id)

        if r1_id is None:
            extra_r2_records += 1
        elif r2_id is None:
            extra_r1_records += 1
        elif r1_id == r2_id:
            position_matches += 1
        else:
            position_mismatches += 1

    return {
        "r1_unique_ids": len(r1_ids),
        "r2_unique_ids": len(r2_ids),
        "complete_pairs": len(r1_ids & r2_ids),
        "r1_only": len(r1_ids - r2_ids),
        "r2_only": len(r2_ids - r1_ids),
        "position_matches": position_matches,
        "position_mismatches": position_mismatches,
        "extra_r1_records": extra_r1_records,
        "extra_r2_records": extra_r2_records,
    }