import pytest
from krakenpy.taxonomy import *
import filecmp
import os


def test_taxonentry():
    """Test TaxonEntry."""
    output = TaxonEntry()
    assert (output.name == "unclassified")
    assert (output.taxon_id == "0")
    assert (output.rank == "U")

    output = TaxonEntry("1", "root", "R")
    assert (output.name == "root")
    assert (output.taxon_id == "1")
    assert (output.rank == "R")

def test_krakenassignmententry_equals():
    """Test equality."""
    entry1 = TaxonEntry("1", "root", "R")
    entry2 = TaxonEntry("1", "root", "R1")
    entry3 = TaxonEntry("2", "root", "R")
    entry4 = TaxonEntry("1", "classified", "R")
    entry5 = "1\tclassified"

    assert (entry1 == entry1)
    assert (entry2 == entry2)
    assert not (entry1 == entry2)
    assert not (entry2 == entry1)
    assert (entry3 == entry3)
    assert not (entry3 == entry1)
    assert not (entry3 == entry2)
    assert not (entry2 == entry3)
    assert not (entry1 == entry4)
    assert not (entry2 == entry4)
    assert not (entry3 == entry4)
    assert not (entry2 == entry5)
    assert not (entry5 == entry4)

def test_taxonentry_print():
    """Test print."""
    entry1 = TaxonEntry("1", "root", "R")
    entry1.print()


def test_taxonomy():
    """Test Taxonomy."""
    output = Taxonomy()
    assert (len(output.entries) == 0)
    assert (len(output.parents) == 0)
    assert (len(output.children) == 0)

    taxonomy_dir = "tests/data/taxonomy"
    output = Taxonomy(taxonomy_dir)
    print(len(output.parents), len(output.children))
    assert (len(output.entries) == 0)
    assert (len(output.parents) == 1792)
    assert (len(output.children) == 659)

    output = Taxonomy(taxonomy_dir, taxon_ids=["2759", "129875", "232094", "131567", "377627", "232100", "9606"])
    print(len(output.entries), len(output.parents), len(output.children))
    assert (len(output.entries) == 1792)
    assert (len(output.parents) == 1792)
    assert (len(output.children) == 659)

def test_taxonomy_equals():
    """Test Taxonomy."""
    entry1 = Taxonomy()
    taxonomy_dir = "tests/data/taxonomy"
    entry2 = Taxonomy(taxonomy_dir)
    entry3 = Taxonomy(taxonomy_dir, taxon_ids=["2759", "129875", "232094", "131567", "377627", "232100", "9606"])
    entry4 = "Im a string"

    assert (entry1 == entry1)
    assert (entry2 == entry2)
    assert not (entry1 == entry2)
    assert not (entry2 == entry1)
    assert (entry3 == entry3)
    assert not (entry3 == entry1)
    assert not (entry3 == entry2)
    assert not (entry2 == entry3)
    assert not (entry4 == entry2)
    assert not (entry2 == entry4)

def test_taxonomy_missing_file():
    """Test missing file."""
    taxonomy_dir = "tests/data/taxonomy_missing"
    with pytest.raises(SystemExit) as e:
        output = Taxonomy(taxonomy_dir)
    assert e.value.code == 4

def test_taxonomy_corrupt_file():
    """Test missing file."""
    taxonomy_dir = "tests/data/taxonomy_missing2"
    taxon_ids = ["2759", "129875", "232094", "131567", "377627", "232100", "9606"]
    with pytest.raises(SystemExit) as e:
        output = Taxonomy(taxonomy_dir, taxon_ids)
    assert e.value.code == 4

def test_taxonomy_missing_file2():
    """Test missing file."""
    taxonomy_dir = "tests/data/taxonomy_missing"
    output = Taxonomy()

    output.load_entries_from_nodes(taxonomy_dir, [])
    assert(len(output.entries) == 0)

    taxon_ids = ["2759", "129875", "232094", "131567", "377627", "232100", "9606"]
    with pytest.raises(SystemExit) as e:
        output.load_entries_from_nodes(taxonomy_dir, taxon_ids)
    assert e.value.code == 4

def test_taxonomy_missing_file3():
    """Test missing file."""
    taxonomy_dir = "tests/data/taxonomy_missing"
    output = Taxonomy()

    output.load_entries_from_nodes(taxonomy_dir, [])
    assert (len(output.entries) == 0)

    taxon_ids = ["2759", "129875", "232094", "131567", "377627", "232100", "9606"]
    with pytest.raises(SystemExit) as e:
        output.load_entries_from_names(taxonomy_dir, taxon_ids)
    assert e.value.code == 4

def test_taxonomy_missing_file4():
    """Test missing file."""
    taxonomy_dir = "tests/data/taxonomy_missing"
    output = Taxonomy()
    with pytest.raises(SystemExit) as e:
        output.load_parents_and_children(taxonomy_dir)
    assert e.value.code == 4

def test_taxonomy_get_taxon_id_map():
    """Test taxon_id_map."""
    taxonomy_dir = "tests/data/taxonomy"
    taxon_ids = ["2759", "129875", "232094", "131567", "377627", "232100", "9606"]
    output = Taxonomy(taxonomy_dir)
    result = output.get_taxon_id_map(taxon_ids)
    expected = {'2759': {'2759', '131567'}, '129875': {'129875'}, '232094': {'232094', '129875', '377627'}, '131567': {'131567'}, '377627': {'129875', '377627'}, '232100': {'129875', '377627', '232100'}, '9606': {'9606'}, '63221': {'9606'}, '741158': {'9606'}, '232108': {'129875', '377627'}, '232096': {'129875', '377627'}, '2157': {'131567'}, '2': {'131567'}, '2698737': {'2759', '131567'}, '2598132': {'2759', '131567'}, '2608240': {'2759', '131567'}, '1401294': {'2759', '131567'}, '2795258': {'2759', '131567'}, '33154': {'2759', '131567'}, '2763': {'2759', '131567'}, '33090': {'2759', '131567'}, '554915': {'2759', '131567'}, '554296': {'2759', '131567'}, '2611352': {'2759', '131567'}, '2608109': {'2759', '131567'}, '2686027': {'2759', '131567'}, '42452': {'2759', '131567'}, '2611341': {'2759', '131567'}, '2683617': {'2759', '131567'}, '38254': {'2759', '131567'}, '2489521': {'2759', '131567'}, '3004206': {'2759', '131567'}, '3027': {'2759', '131567'}, '61964': {'2759', '131567'}, '28282': {'129875'}, '1509400': {'129875'}, '1069441': {'129875'}, '10528': {'129875'}, '10529': {'129875'}}
    assert(result == expected)

    result = output.get_taxon_id_map(taxon_ids, include_unclassified=True)
    expected = {'0':{"2759", "129875", "232094", "131567", "377627", "232100", "9606"}, '2759': {'2759', '131567'}, '129875': {'129875'}, '232094': {'232094', '129875', '377627'},
                '131567': {'131567'}, '377627': {'129875', '377627'}, '232100': {'129875', '377627', '232100'},
                '9606': {'9606'}, '63221': {'9606'}, '741158': {'9606'}, '232108': {'129875', '377627'},
                '232096': {'129875', '377627'}, '2157': {'131567'}, '2': {'131567'}, '2698737': {'2759', '131567'},
                '2598132': {'2759', '131567'}, '2608240': {'2759', '131567'}, '1401294': {'2759', '131567'},
                '2795258': {'2759', '131567'}, '33154': {'2759', '131567'}, '2763': {'2759', '131567'},
                '33090': {'2759', '131567'}, '554915': {'2759', '131567'}, '554296': {'2759', '131567'},
                '2611352': {'2759', '131567'}, '2608109': {'2759', '131567'}, '2686027': {'2759', '131567'},
                '42452': {'2759', '131567'}, '2611341': {'2759', '131567'}, '2683617': {'2759', '131567'},
                '38254': {'2759', '131567'}, '2489521': {'2759', '131567'}, '3004206': {'2759', '131567'},
                '3027': {'2759', '131567'}, '61964': {'2759', '131567'}, '28282': {'129875'}, '1509400': {'129875'},
                '1069441': {'129875'}, '10528': {'129875'}, '10529': {'129875'}}
    assert (result == expected)
