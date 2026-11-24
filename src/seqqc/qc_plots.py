
# Import Matplotlib for creating scientific plots
import matplotlib.pyplot as plt

# Import our existing FASTQ parser
from seqqc.fastq_parser import read_fastq


def plot_read_lengths(r1_path, r2_path, output_path):
    """
    Compare read-length distributions from two FASTQ files.

    Parameters
    ----------
    r1_path : str
        Path to the Read 1 FASTQ file.
    r2_path : str
        Path to the Read 2 FASTQ file.
    output_path : str
        Location where the plot will be saved.
    """

    # Store the read lengths from R1 and R2 separately
    r1_lengths = []
    r2_lengths = []

    # Read each R1 FASTQ record and calculate its length
    for read_id, sequence, quality in read_fastq(r1_path):
        r1_lengths.append(len(sequence))

    # Read each R2 FASTQ record and calculate its length
    for read_id, sequence, quality in read_fastq(r2_path):
        r2_lengths.append(len(sequence))

    # Create a new plotting area
    plt.figure(figsize=(10, 6))

    # Import NumPy at the top of the Python file
    import numpy as np

    # Use shared 5-bp intervals for both read distributions
    maximum_length = max(max(r1_lengths), max(r2_lengths))

    bins = np.arange(0, maximum_length + 6, 5)

    # Plot R1 read-length distribution
    plt.hist(
        r1_lengths,
        bins=bins,
        alpha=0.6,
        label="Read 1",
        density=True,
    )

    # Plot R2 read-length distribution
    plt.hist(
        r2_lengths,
        bins=bins,
        alpha=0.6,
        label="Read 2",
        density=True,
    )

    # Add descriptive labels and a title
    plt.xlabel("Read length (bp)")
    plt.ylabel("Probability density")
    plt.title("Illumina MiSeq Read-Length Distribution")

    # Display the legend
    plt.legend()

    # Adjust layout to prevent text from being cut off
    plt.tight_layout()

    # Save the plot as a high-resolution image
    plt.savefig(output_path, dpi=300)

    # Close the figure to free memory
    plt.close()
