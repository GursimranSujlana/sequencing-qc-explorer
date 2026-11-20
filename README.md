
# SeqQC Explorer

A Python toolkit for quality assessment and analysis of Illumina and Oxford Nanopore sequencing data.

## Project Status

Early development (v0.1.0).

## Current Features

- FASTQ record parsing
- Phred+33 quality-score conversion
- GC-content calculation
- Mean Phred quality calculation
- Q20 and Q30 percentage calculations
- Automated unit tests using pytest

## Planned Features

- Sample-level sequencing QC summaries
- Read-length and quality distributions
- Compressed FASTQ support
- Illumina and ONT comparative analysis
- QC visualisations and reports
- Validation against established QC tools

## Installation

Create a Python 3.12 environment and run:

```bash
python -m pip install -e ".[dev]"
```

## Testing

```bash
python -m pytest tests/ -v
```

## Data

The repository uses small synthetic datasets for testing. Larger public sequencing datasets will be downloaded separately using documented instructions.

## License

To be determined.
