import pytest
from krakenpy.assignment import *
from krakenpy.taxonomy import *
import filecmp
import os


def test_trim_read_id1():
    """Test trim_read_id."""
    read_id = "example/1"
    expected = "example"
    output = trim_read_id(read_id)
    assert(expected == output)

    read_id = "example/2"
    expected = "example"
    output = trim_read_id(read_id)
    assert (expected == output)

    read_id = "example/12"
    expected = "example/12"
    output = trim_read_id(read_id)
    assert (expected == output)

    read_id = "example.1"
    expected = "example.1"
    output = trim_read_id(read_id)
    assert (expected == output)

def test_krakenassignmententry():
    """Test KrakenAssignmentEntry."""
    output = KrakenAssignmentEntry()
    assert (output.classified == "U")
    assert (output.read_id == "")
    assert (output.taxon_id == "0")
    assert (output.length == 0)
    assert (output.kmer_string == "")

    line = ("C\tcadc9752-bcc4-af2c-be48-d30a9f06e364\t2748958\t6306\t0:31 2748958:2 0:45 1003835:3 0:82 1003835:2 0:35")
    output = KrakenAssignmentEntry(line)
    assert (output.classified == "C")
    assert (output.read_id == "cadc9752-bcc4-af2c-be48-d30a9f06e364")
    assert (output.taxon_id == "2748958")
    assert (output.length == 6306)
    assert (output.kmer_string == "0:31 2748958:2 0:45 1003835:3 0:82 1003835:2 0:35")

    line = ("U\tad88b02c-8dc5-c9cd-4e62-33270ccb9b2f/1\t0\t653\t0:619")
    output = KrakenAssignmentEntry(line)
    assert (output.classified == "U")
    assert (output.read_id == "ad88b02c-8dc5-c9cd-4e62-33270ccb9b2f")
    assert (output.taxon_id == "0")
    assert (output.length == 653)
    assert (output.kmer_string == "0:619")

    line = ("C\tartificial_read\tA\t653\t81077:619")
    output = KrakenAssignmentEntry(line)
    assert (output.classified == "C")
    assert (output.read_id == "artificial_read")
    assert (output.taxon_id == "81077")
    assert (output.length == 653)
    assert (output.kmer_string == "81077:619")

def test_krakenassignmententry_equals():
    """Test equality."""
    line = ("C\tartificial_read\tA\t653\t81077:619")
    entry1 = KrakenAssignmentEntry(line)
    line = ("U\tad88b02c-8dc5-c9cd-4e62-33270ccb9b2f/1\t0\t653\t0:619")
    entry2 = KrakenAssignmentEntry(line)
    assert (entry1 == entry1)
    assert (entry2 == entry2)
    assert not (entry1 == entry2)
    assert not (entry2 == entry1)

    line = ("U\tad88b02c-8dc5-c9cd-4e62-33270ccb9b2f/1\t0\t654\t0:619")
    entry3 = KrakenAssignmentEntry(line)
    assert (entry3 == entry3)
    assert not (entry3 == entry1)
    assert not (entry3 == entry2)
    assert not (entry2 == entry3)

    entry4 = "U\tad88b02c-8dc5-c9cd-4e62-33270ccb9b2f/1\t0\t654\t0:619"
    assert not (entry1 == entry4)
    assert not (entry2 == entry4)
    assert not (entry3 == entry4)

def test_krakenassignmententry_print():
    """Test print."""
    line = ("C\tartificial_read\tA\t653\t81077:619")
    entry1 = KrakenAssignmentEntry(line)
    entry1.print()

def test_krakenassignmententry_declassify():
    """Test printing."""
    line = ("C\tcadc9752-bcc4-af2c-be48-d30a9f06e364\t2748958\t6306\t0:31 2748958:2 0:45 1003835:3 0:82 1003835:2 0:35")
    output = KrakenAssignmentEntry(line)
    output.declassify()
    assert (output.classified == "U")
    assert (output.read_id == "cadc9752-bcc4-af2c-be48-d30a9f06e364")
    assert (output.taxon_id == "0")
    assert (output.length == 6306)
    assert (output.kmer_string == "0:31 2748958:2 0:45 1003835:3 0:82 1003835:2 0:35")

def test_krakenassignmententry_line():
    """Test printing."""
    line = "C\tcadc9752-bcc4-af2c-be48-d30a9f06e364\t2748958\t6306\t0:31 2748958:2 0:45 1003835:3 0:82 1003835:2 0:35"
    output = KrakenAssignmentEntry(line)
    output_line = output.get_line()
    assert (output_line == line)

    output = KrakenAssignmentEntry()
    output.add_line(line)
    output_line = output.get_line()
    assert (output_line == line)

def test_krakenassignmententry_bad_line():
    """Test printing."""
    line = "C\tcadc9752-bcc4-af2c-be48-d30a9f06e364\t2748958  6306\t0:31 2748958:2 0:45 1003835:3 0:82 1003835:2 0:35"
    with pytest.raises(SystemExit) as e:
        output = KrakenAssignmentEntry(line)
    assert e.value.code == 11

def test_krakenassignments():
    """Test KrakenAssignments."""
    file_name = "not_a_real_file"
    output = KrakenAssignments(file_name)
    assert (len(output.entries) == 0)
    assert (output.file_name == file_name)

    input_prefix = "tests/data/taxid_630"
    input_assignment = f"{input_prefix}/Viral.kraken_assignments.tsv"
    output = KrakenAssignments(input_assignment)
    assert (len(output.entries) == 0)
    assert (output.file_name == input_assignment)

    output = KrakenAssignments(input_assignment, load=True)
    assert (len(output.entries) == 2500)
    assert (output.file_name == input_assignment)
    assert ("2deb3d12-e44f-e13b-0076-74f80fe13193" in output.entries)
    expected_entry = KrakenAssignmentEntry("U\t2deb3d12-e44f-e13b-0076-74f80fe13193\t0\t1168\t0:1134")
    assert(output.entries["2deb3d12-e44f-e13b-0076-74f80fe13193"] == expected_entry)

def test_krakenassignments_load():
    """Test KrakenAssignments load function."""
    input_prefix = "tests/data/paired"
    input_assignment = f"{input_prefix}/Viral.kraken_assignments.tsv"
    output = KrakenAssignments(input_assignment, load=True)

    assert (len(output.entries) == 1000)
    assert (output.file_name == input_assignment)
    assert ("Human_adenovirus_A|129875_9628_10178_2:0:0_5:0:0_0" in output.entries)
    expected_entry = KrakenAssignmentEntry("U\tHuman_adenovirus_A|129875_29788_30302_1:0:0_2:0:0_1/1\t0\t301\t0:116 A:34 0:117")
    assert (output.entries["Human_adenovirus_A|129875_29788_30302_1:0:0_2:0:0_1"] == expected_entry)

    input_assignment = f"{input_prefix}/small.kraken_assignments.tsv"
    output = KrakenAssignments(input_assignment, load=True)

    assert (len(output.entries) == 7)
    assert (output.file_name == input_assignment)
    classifieds = ["U","C","C","C","C","U","U","U"]
    taxon_ids = ["0", "81077", "129875", "1", "1","0","0","0"]
    for i in range(7):
        assert (f"Human_adenovirus_A|129875_{i}" in output.entries)
        assert (output.entries[f"Human_adenovirus_A|129875_{i}"].classified == classifieds[i])
        assert (output.entries[f"Human_adenovirus_A|129875_{i}"].taxon_id == taxon_ids[i])

def test_krakenassignments_equals():
    """Test equality."""
    input_prefix = "tests/data/taxid_630"
    input_assignment = f"{input_prefix}/Viral.kraken_assignments.tsv"
    entry1 = KrakenAssignments(input_assignment)
    entry2 = KrakenAssignments(input_assignment, load=True)

    input_prefix = "tests/data/taxid_630"
    input_assignment = f"{input_prefix}/PlusPF-8.kraken_assignments.tsv"
    entry3 = KrakenAssignments(input_assignment, load=True)

    entry4 = KrakenAssignmentEntry()

    assert (entry1 == entry1)
    assert (entry2 == entry2)
    assert (entry3 == entry3)
    assert not (entry1 == entry2)
    assert not (entry2 == entry1)
    assert not (entry1 == entry3)
    assert not (entry3 == entry1)
    assert not (entry1 == entry4)
    assert not (entry4 == entry1)
    assert not (entry2 == entry3)
    assert not (entry3 == entry2)

def test_krakenassignments_update():
    """Test KrakenAssignments load function."""
    input_prefix = "tests/data/paired"
    input_assignment = f"{input_prefix}/small.kraken_assignments.tsv"
    output = KrakenAssignments(input_assignment, load=True)

    input_assignment = f"{input_prefix}/additional.kraken_assignments.tsv"
    new = KrakenAssignments(input_assignment, load=True)

    changes = output.update(new)
    assert (len(output.entries) == 9)

    expected_changes = {"0":{"20":2,"0":1},"81077":{"20":1}}
    assert(changes == expected_changes)

    classifieds = ["C","C","C","C","C","U","U","C","U"]
    taxon_ids = ["20", "20", "129875", "1", "1","0","0","20","0"]
    for i in range(8):
        print(i)
        assert (f"Human_adenovirus_A|129875_{i}" in output.entries)
        assert (output.entries[f"Human_adenovirus_A|129875_{i}"].classified == classifieds[i])
        assert (output.entries[f"Human_adenovirus_A|129875_{i}"].taxon_id == taxon_ids[i])

def test_krakenassignments_save():
    input_prefix = "tests/data/taxid_630"
    input_assignment = f"{input_prefix}/Viral.kraken_assignments.tsv"
    output = KrakenAssignments(input_assignment, load=True)
    out_assignment = f"{input_prefix}/test.kraken_assignments.tsv"
    output.file_name = out_assignment
    output.save()

    expected = input_assignment
    assert (filecmp.cmp(out_assignment, expected, shallow=False))
    os.unlink(out_assignment)

def test_update_and_save():
    input_prefix = "tests/data/paired"
    input_assignment = f"{input_prefix}/small.kraken_assignments.tsv"
    output = KrakenAssignments(input_assignment, load=True)
    input_assignment = f"{input_prefix}/additional.kraken_assignments.tsv"
    new = KrakenAssignments(input_assignment, load=True)
    output.update(new)

    out_assignment = f"{input_prefix}/test.kraken_assignments.tsv"
    output.file_name = out_assignment
    output.save()

    expected = f"{input_prefix}/expected_small_and_additional.kraken_assignments.tsv"
    assert (filecmp.cmp(out_assignment, expected, shallow=False))
    os.unlink(out_assignment)

def test_update_and_save2():
    input_prefix = "tests/data/taxid_630"
    input_assignment = f"{input_prefix}/Viral.kraken_assignments.tsv"
    output = KrakenAssignments(input_assignment, load=True)
    input_assignment = f"{input_prefix}/PlusPF-8.kraken_assignments.tsv"
    new = KrakenAssignments(input_assignment, load=True)
    output.update(new)

    out_assignment = f"{input_prefix}/test.kraken_assignments.tsv"
    output.file_name = out_assignment
    output.save()

    expected = f"{input_prefix}/expected_merged_inverted.kraken_assignments.tsv"
    assert (filecmp.cmp(out_assignment, expected, shallow=False))
    os.unlink(out_assignment)
def test_update_and_save3():
    input_prefix = "tests/data/taxid_630"
    input_assignment = f"{input_prefix}/PlusPF-8.kraken_assignments.tsv"
    output = KrakenAssignments(input_assignment, load=True)
    input_assignment = f"{input_prefix}/Viral.kraken_assignments.tsv"
    new = KrakenAssignments(input_assignment, load=True)
    output.update(new)

    out_assignment = f"{input_prefix}/test.kraken_assignments.tsv"
    output.file_name = out_assignment
    output.save()

    expected = f"{input_prefix}/PlusPF-8.kraken_assignments.tsv"
    assert (filecmp.cmp(out_assignment, expected, shallow=False))
    os.unlink(out_assignment)

def test_krakenassignments_get_read_map():
    """Test KrakenAssignments get_read_map function."""
    input_prefix = "tests/data/taxid_630"
    input_assignment = f"{input_prefix}/PlusPF-8.kraken_assignments.tsv"
    output = KrakenAssignments(input_assignment, load=True)
    taxon_ids = ["9605"]
    read_map = output.get_read_map(taxon_ids)
    assert (len(read_map) == 0)

    loaded_taxonomy = Taxonomy("/Users/rmcolq/Work/git/scylla/store_dir/taxonomy_dir", taxon_ids)
    read_map = output.get_read_map(taxon_ids, loaded_taxonomy.parents)
    print(read_map)
    assert (len(read_map) == 6)

    taxon_ids = ["9606"]
    read_map = output.get_read_map(taxon_ids)
    assert (len(read_map) == 6)

    input_prefix = "tests/data/paired"
    input_assignment = f"{input_prefix}/small.kraken_assignments.edited.tsv"
    output = KrakenAssignments(input_assignment, load=True)
    taxon_ids = ["129875", "1", "10528", "81077"]
    read_map = output.get_read_map(taxon_ids, loaded_taxonomy.parents)
    print(read_map)
    expected_read_map = {'Human_adenovirus_A|129875_1': '81077', 'Human_adenovirus_A|129875_2': '129875',
                         'Human_adenovirus_A|129875_3': '1', 'Human_adenovirus_A|129875_4': '1',
                         'Human_adenovirus_A|129875_5': '1', 'Human_adenovirus_A|129875_7': '1',
                         'Human_adenovirus_A|129875_9': '1', 'Human_adenovirus_A|129875_10': '129875',
                         'Human_adenovirus_A|129875_11': '129875'}
    assert (read_map == expected_read_map)
    taxon_ids = ["129875", "2", "81077"]
    read_map = output.get_read_map(taxon_ids, loaded_taxonomy.parents)
    print(read_map)
    expected_read_map = {'Human_adenovirus_A|129875_1': '81077', 'Human_adenovirus_A|129875_2': '129875',
                         'Human_adenovirus_A|129875_10': '129875', 'Human_adenovirus_A|129875_11': '129875'}
    assert (read_map == expected_read_map)