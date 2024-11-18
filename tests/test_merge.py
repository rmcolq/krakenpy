import pytest
from krakenpy import merge
import filecmp
import os


def test_second_unclassified():
    """Test merge when second file pair all unclassified."""
    input_prefix = "tests/data/taxid_630"
    input_assignment1 = f"{input_prefix}/PlusPF-8.kraken_assignments.tsv"
    input_assignment2 = f"{input_prefix}/Viral.kraken_assignments.tsv"
    input_report1 = f"{input_prefix}/PlusPF-8.kraken_report.txt"
    input_report2 = f"{input_prefix}/Viral.kraken_report.txt"


    expected_report = f"{input_prefix}/PlusPF-8.kraken_report.txt"
    expected_assignment = f"{input_prefix}/PlusPF-8.kraken_assignments.tsv"

    output_prefix = f"{input_prefix}/merged"
    out_assignment = f"{output_prefix}.kraken_assignments.tsv"
    out_report = f"{output_prefix}.kraken_report.txt"
    merge.merge([input_assignment1, input_assignment2], [input_report1, input_report2], output_prefix)

    #assert(filecmp.cmp(out_assignment, expected_assignment, shallow=False))
    assert(filecmp.cmp(out_report, expected_report, shallow=False))

    #os.unlink(out_assignment)
    os.unlink(out_report)

def test_second_unclassified_inverted():
    """Test merge when second file pair all unclassified."""
    input_prefix = "tests/data/taxid_630"
    input_assignment1 = f"{input_prefix}/PlusPF-8.kraken_assignments.tsv"
    input_assignment2 = f"{input_prefix}/Viral.kraken_assignments.tsv"
    input_report1 = f"{input_prefix}/PlusPF-8.kraken_report.txt"
    input_report2 = f"{input_prefix}/Viral.kraken_report.txt"


    expected_report = f"{input_prefix}/PlusPF-8.kraken_report.txt"
    expected_assignment = f"{input_prefix}/PlusPF-8.kraken_assignments.tsv"

    output_prefix = f"{input_prefix}/merged_inverted"
    out_assignment = f"{output_prefix}.kraken_assignments.tsv"
    out_report = f"{output_prefix}.kraken_report.txt"
    merge.merge([input_assignment2, input_assignment1], [input_report2, input_report1], output_prefix)

    #assert(filecmp.cmp(out_assignment, expected_assignment, shallow=False))
    assert(filecmp.cmp(out_report, expected_report, shallow=False))

    #os.unlink(out_assignment)
    os.unlink(out_report)

def test_second_more_precise():
    """Test merge when second file pair gives an additional level of specificity."""
    input_prefix = "tests/data/taxid_1003835"
    input_assignment1 = f"{input_prefix}/PlusPF-8.kraken_assignments.tsv"
    input_assignment2 = f"{input_prefix}/Viral.kraken_assignments.tsv"
    input_report1 = f"{input_prefix}/PlusPF-8.kraken_report.txt"
    input_report2 = f"{input_prefix}/Viral.kraken_report.txt"

    expected_report = f"{input_prefix}/Viral.kraken_report.txt"
    expected_assignment = f"{input_prefix}/Viral.kraken_assignments.tsv"

    output_prefix = f"{input_prefix}/merged"
    #out_assignment = f"{output_prefix}.kraken_assignments.tsv"
    out_report = f"{output_prefix}.kraken_report.txt"
    merge.merge([input_assignment1, input_assignment2], [input_report1, input_report2], output_prefix)

    #assert (filecmp.cmp(out_assignment, expected_assignment, shallow=False))
    assert (filecmp.cmp(out_report, expected_report, shallow=False))

    #os.unlink(out_assignment)
    os.unlink(out_report)

def test_second_more_precise_inverted():
    """Test merge when second file pair gives an additional level of specificity."""
    input_prefix = "tests/data/taxid_1003835"
    input_assignment1 = f"{input_prefix}/PlusPF-8.kraken_assignments.tsv"
    input_assignment2 = f"{input_prefix}/Viral.kraken_assignments.tsv"
    input_report1 = f"{input_prefix}/PlusPF-8.kraken_report.txt"
    input_report2 = f"{input_prefix}/Viral.kraken_report.txt"

    output_prefix = f"{input_prefix}/merged_inverted"
    out_assignment = f"{output_prefix}.kraken_assignments.tsv"
    out_report = f"{output_prefix}.kraken_report.txt"
    merge.merge([input_assignment2, input_assignment1], [input_report2, input_report1], output_prefix)

    expected_report = f"{input_prefix}/expected.kraken_report.txt"
    expected_assignment = f"{input_prefix}/expected.kraken_assignments.tsv"

    #assert (filecmp.cmp(out_assignment, expected_assignment, shallow=False))
    assert (filecmp.cmp(out_report, expected_report, shallow=False))

    #os.unlink(out_assignment)
    os.unlink(out_report)