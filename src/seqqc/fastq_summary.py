
# Import the FASTQ parser we developed earlier.
from seqqc.fastq_parser import read_fastq

# Import our Phred quality conversion function.
from seqqc.qc_metrics import phred_scores


def summarize_fastq(file_path):
    """
    Calculate quality-control statistics for a FASTQ file.

    Parameters
    ----------
    file_path : str
        Path to the FASTQ file.

    Returns
    -------
    dict
        Dictionary containing sequencing QC statistics.
    """

    # Initialise counters before reading the file.
    total_reads = 0
    total_bases = 0
    total_gc_bases = 0
    total_valid_bases = 0
    total_quality = 0
    q20_bases = 0
    q30_bases = 0

    # Process one sequencing record at a time.
    for read_id, sequence, quality in read_fastq(file_path):

        # Count each read.
        total_reads += 1

        # Count nucleotides in the read.
        total_bases += len(sequence)

        # Convert the sequence to uppercase for consistent counting.
        sequence = sequence.upper()

        # Count G and C nucleotides.
        total_gc_bases += sequence.count("G") + sequence.count("C")

        # Count unambiguous DNA bases for GC calculation.
        total_valid_bases += sum(
            sequence.count(base) for base in "ATGC"
        )

        # Convert quality characters to Phred scores.
        scores = phred_scores(quality)

        # Add individual Phred scores to the running total.
        total_quality += sum(scores)

        # Count bases with Phred scores >= 20.
        q20_bases += sum(score >= 20 for score in scores)

        # Count bases with Phred scores >= 30.
        q30_bases += sum(score >= 30 for score in scores)

    # An empty FASTQ file cannot produce a meaningful summary.
    if total_reads == 0:
        raise ValueError("FASTQ file contains no reads")

    # Calculate the summary statistics.
    summary = {
        "read_count": total_reads,
        "total_bases": total_bases,
        "mean_read_length": total_bases / total_reads,

        # Exclude ambiguous bases from the GC denominator.
        "gc_percentage": (
            100 * total_gc_bases / total_valid_bases
            if total_valid_bases else 0.0
        ),

        # Average Phred score across all sequenced bases.
        "mean_phred_quality": (
            total_quality / total_bases
            if total_bases else 0.0
        ),

        # Percentage of bases meeting Q20.
        "q20_percentage": (
            100 * q20_bases / total_bases
            if total_bases else 0.0
        ),

        # Percentage of bases meeting Q30.
        "q30_percentage": (
            100 * q30_bases / total_bases
            if total_bases else 0.0
        ),
    }

    return summary
