
# Import the functions we developed
from seqqc.qc_metrics import (
    phred_scores,
    calculate_gc_content,
    calculate_mean_quality,
    calculate_q_percentage,
)


def test_phred_conversion():
    # I = Q40, 5 = Q20, ! = Q0
    assert phred_scores("I5!") == [40, 20, 0]


def test_gc_content():
    # Two of the four known bases are G or C
    assert calculate_gc_content("ATGC") == 50.0


def test_gc_with_ambiguous_bases():
    # Our default method excludes N from the denominator
    assert calculate_gc_content("ATGCNN") == 50.0


def test_mean_quality():
    # (40 + 20 + 0) / 3 = 20
    assert calculate_mean_quality("I5!") == 20.0


def test_q20_percentage():
    # Two of three bases have quality >= 20
    result = calculate_q_percentage("I5!", 20)
    assert round(result, 2) == 66.67


def test_q30_percentage():
    # One of three bases has quality >= 30
    result = calculate_q_percentage("I5!", 30)
    assert round(result, 2) == 33.33
