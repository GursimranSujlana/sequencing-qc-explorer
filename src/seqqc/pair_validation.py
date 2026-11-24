
# Import the FASTQ parser we already developed.
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
    """
    Compare the spot identifiers in two FASTQ files.

    Returns counts of matching and unmatched reads.
    """

    # Collect identifiers from Read 1.
    r1_ids = set()

    for read_id, sequence, quality in read_fastq(r1_path):
        r1_ids.add(get_spot_id(read_id))

    # Collect identifiers from Read 2.
    r2_ids = set()

    for read_id, sequence, quality in read_fastq(r2_path):
        r2_ids.add(get_spot_id(read_id))

    # Find identifiers present in both files.
    paired_ids = r1_ids.intersection(r2_ids)

    # Find identifiers appearing only in R1.
    r1_only = r1_ids.difference(r2_ids)

    # Find identifiers appearing only in R2.
    r2_only = r2_ids.difference(r1_ids)

    # Return summary statistics.
    return {
        "r1_unique_ids": len(r1_ids),
        "r2_unique_ids": len(r2_ids),
        "complete_pairs": len(paired_ids),
        "r1_only": len(r1_only),
        "r2_only": len(r2_only),
    }
