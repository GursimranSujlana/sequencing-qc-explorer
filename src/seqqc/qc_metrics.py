
def phred_scores(quality_string):
    """
    Convert FASTQ Phred+33 quality characters
    into numerical Phred quality scores.

    Example:
        "I5!" -> [40, 20, 0]
    """

    # Create an empty list to store numerical scores
    scores = []

    # Process each character in the quality string
    for character in quality_string:

        # ord() returns the character's ASCII code
        ascii_value = ord(character)

        # Subtract 33 for Phred+33 encoding
        quality_score = ascii_value - 33

        # Store the calculated quality score
        scores.append(quality_score)

    # Return the numerical quality scores
    return scores


def calculate_gc_content(sequence):
    """
    Calculate GC percentage in a DNA sequence.

    Only A, T, G and C are included in the denominator.
    Ambiguous bases such as N are excluded.
    """

    # Convert lowercase nucleotides to uppercase.
    sequence = sequence.upper()

    # Count G and C nucleotides.
    gc_count = sequence.count("G") + sequence.count("C")

    # Count only unambiguous DNA bases.
    valid_bases = sum(sequence.count(base) for base in "ATGC")

    # Avoid division by zero.
    if valid_bases == 0:
        return 0.0

    # Calculate GC percentage.
    return (gc_count / valid_bases) * 100


def calculate_mean_quality(quality_string):
    """
    Calculate the arithmetic mean of Phred+33 scores.
    """

    # Convert quality characters to numerical scores.
    scores = phred_scores(quality_string)

    # Handle an empty quality string.
    if not scores:
        return 0.0

    # Calculate arithmetic mean.
    return sum(scores) / len(scores)


def calculate_q_percentage(quality_string, threshold):
    """
    Calculate the percentage of bases meeting
    a minimum Phred quality threshold.
    """

    # Convert ASCII quality characters to numerical scores.
    scores = phred_scores(quality_string)

    # Prevent division by zero.
    if not scores:
        return 0.0

    # Count bases meeting or exceeding the threshold.
    passing_bases = sum(
        score >= threshold for score in scores
    )

    # Return percentage.
    return (passing_bases / len(scores)) * 100


import numpy as np


def calculate_length_statistics(read_lengths):
    """
    Calculate read-length summary statistics.

    Parameters
    ----------
    read_lengths : list
        Lengths of sequencing reads.

    Returns
    -------
    dict
        Minimum, maximum, mean, median,
        standard deviation and read N50.
    """

    # Convert read lengths to a NumPy array.
    lengths = np.asarray(read_lengths)

    # Check whether the list is empty.
    if lengths.size == 0:
        raise ValueError("No read lengths provided")

    # Sort read lengths from longest to shortest.
    sorted_lengths = np.sort(lengths)[::-1]

    # Calculate cumulative sequenced bases.
    cumulative_bases = np.cumsum(sorted_lengths)

    # Find the position where half the bases are covered.
    halfway = np.sum(sorted_lengths) / 2

    # Find the first length reaching the halfway point.
    n50_index = np.searchsorted(cumulative_bases, halfway)
    n50 = sorted_lengths[n50_index]

    return {
        "minimum": int(np.min(lengths)),
        "maximum": int(np.max(lengths)),
        "mean": float(np.mean(lengths)),
        "median": float(np.median(lengths)),
        "standard_deviation": float(np.std(lengths)),
        "n50": int(n50),
    }
