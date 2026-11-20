
def read_fastq(file_path):
    """
    Read records from an uncompressed FASTQ file.

    Each FASTQ record contains four lines:
    1. Read identifier
    2. Nucleotide sequence
    3. Separator (+)
    4. Quality string
    """

    # Open the file in reading mode.
    # 'with' automatically closes the file when finished.
    with open(file_path, "r", encoding="utf-8") as handle:

        # Continue reading until the end of the FASTQ file.
        while True:

            # Read the four lines belonging to one record.
            header = handle.readline()
            if not header:
                break  # End of file

            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()

            # Check that all four lines are present.
            if not all([sequence, separator, quality]):
                raise ValueError("Incomplete FASTQ record")

            # FASTQ headers must begin with @.
            if not header.startswith("@"):
                raise ValueError("FASTQ header must start with @")

            # The separator line must begin with +.
            if not separator.startswith("+"):
                raise ValueError("FASTQ separator must start with +")

            # Remove newline characters.
            read_id = header.rstrip("\r\n")[1:]
            sequence = sequence.rstrip("\r\n")
            quality = quality.rstrip("\r\n")

            # Every nucleotide needs one quality character.
            if len(sequence) != len(quality):
                raise ValueError(
                    "Sequence and quality lengths differ"
                )

            # Return one record at a time.
            yield read_id, sequence, quality
