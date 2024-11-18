import pytest
from krakenpy.report import *


def test_krakenentry():
    """Test KrakenEntry."""
    output = KrakenEntry()
    assert (output.taxon_id == "0")
    assert (output.name == "unclassified")
    assert (output.rank == "U")
    assert (output.depth == 0)
    assert (output.count == 0)
    assert (output.ucount == 0)
    assert (output.domain == None)
    assert (output.parent == None)
    assert (output.children == set())
    assert (output.sibling_rank == 0)
    assert (output.hierarchy == [])

    row={"% of Seqs":4.0,"Clades":20,"Taxonomies":10,"Rank":"S","Taxonomy ID":"630","Scientific Name":"          Yersinia enterocolitica"}
    output = KrakenEntry(row)
    assert (output.taxon_id == "630")
    assert (output.name == "Yersinia enterocolitica")
    assert (output.rank == "S")
    assert (output.depth == 5)
    assert (output.count == 20)
    assert (output.ucount == 10)
    assert (output.domain == None)
    assert (output.parent == None)
    assert (output.children == set())
    assert (output.sibling_rank == 0)
    assert (output.hierarchy == [])

    row = {"% of Seqs": 4.0, "Clades": 20, "Taxonomies": 10, "Rank": "S", "Taxonomy ID": "630",
           "Scientific Name": "          Yersinia enterocolitica"}
    output = KrakenEntry(row, domain="Bacteria", hierarchy=["1","2","1224","1236","629"])
    assert (output.domain == "Bacteria")
    assert (output.hierarchy == ["1","2","1224","1236","629"])

    row = {"% of Seqs": 4.0, "Clades": 20, "Taxonomies": 100, "Rank": "S", "Taxonomy ID": "630",
           "Scientific Name": "          Yersinia enterocolitica"}
    output = KrakenEntry(row)
    assert (output.count == 100)
    assert (output.ucount == 20)

def test_krakenentry_equals():
    row = {"% of Seqs": 4.0, "Clades": 20, "Taxonomies": 10, "Rank": "S", "Taxonomy ID": "630",
           "Scientific Name": "          Yersinia enterocolitica"}
    entry1 = KrakenEntry(row)
    row = {"% of Seqs": 4.1, "Clades": 30, "Taxonomies": 10, "Rank": "S", "Taxonomy ID": "630",
           "Scientific Name": "          Yersinia enterocolitica"}
    entry2 = KrakenEntry(row)

    assert(entry1 == entry1)
    assert (entry2 == entry2)
    assert (entry1 != entry2)

def test_krakenentry_update():
    row = {"% of Seqs": 4.0, "Clades": 20, "Taxonomies": 10, "Rank": "S", "Taxonomy ID": "630",
           "Scientific Name": "            Yersinia enterocolitica"}
    entry1 = KrakenEntry(row, domain="Bacteria", hierarchy=["1","2","1224","1236","629"])

    with pytest.raises(AssertionError):
        entry2 = KrakenEntry(row, domain="Viruses", hierarchy=["1","2","1224","1236","629"])
        entry1.update(entry2)

    with pytest.raises(AssertionError):
        entry2 = KrakenEntry(row, domain="Bacteria", hierarchy=["1","2","1224","629"])
        entry1.update(entry2)

    with pytest.raises(AssertionError):
        entry2 = KrakenEntry(row, domain="Bacteria", hierarchy=["1","2","1224","629"])
        entry2.add_parent("629")
        entry1.update(entry2)

    with pytest.raises(AssertionError):
        entry2 = KrakenEntry(row, domain="Bacteria", hierarchy=["1","2","1224","629"])
        entry2.depth == 1
        entry1.update(entry2)

    row = {"% of Seqs": 4.1, "Clades": 30, "Taxonomies": 10, "Rank": "S1", "Taxonomy ID": "630",
           "Scientific Name": "            Yersinia enterocolitica2"}
    entry2 = KrakenEntry(row, domain="Bacteria", hierarchy=["1","2","1224","1236","629"])
    entry2.add_child("9999")
    entry1.update(entry2)

    assert (entry1.name == "Yersinia enterocolitica2")
    assert (entry1.rank == "S1")
    assert (entry1.children == set(["9999"]))

def test_krakenreport():
    """Test KrakenReport."""
    output = KrakenReport()
    assert (len(output.entries) == 0)
    assert (output.total == 0)
    assert (output.unclassified == 0)
    assert (output.classified == 0)
    assert (len(output.domains) == 0)
    assert (output.file_name == None)

    input_prefix = "tests/data/paired"
    input_report = f"{input_prefix}/Viral.kraken_report.txt"
    output = KrakenReport(input_report)
    assert (len(output.entries) == 13)
    assert (output.total == 1000)
    assert (output.unclassified == 775)
    assert (output.classified == 225)
    assert (len(output.domains) == 1)
    assert (output.file_name == input_report)

    input_prefix = "tests/data/taxid_630"
    input_report = f"{input_prefix}/Viral.kraken_report.txt"
    output = KrakenReport(input_report)
    assert (len(output.entries) == 1)
    assert (output.total == 2500)
    assert (output.unclassified == 2500)
    assert (output.classified == 0)
    assert (len(output.domains) == 0)
    assert (output.file_name == input_report)

def test_krakenreport_equals():
    """Test KrakenReport."""
    input_prefix = "tests/data/paired"
    input_report = f"{input_prefix}/Viral.kraken_report.txt"
    output1 = KrakenReport(input_report)

    input_prefix = "tests/data/taxid_630"
    input_report = f"{input_prefix}/Viral.kraken_report.txt"
    output2 = KrakenReport(input_report)
    assert (output1 == output1)
    assert (output2 == output2)
    assert (output1 != output2)
    assert (output2 != output1)

def test_krakenreport_add_parent_child():
    """Test KrakenReport."""
    input_prefix = "tests/data/paired"
    input_report = f"{input_prefix}/Viral.kraken_report.txt"
    output = KrakenReport(input_report)

    output.add_parent_child("10509","129875")
    assert("129875" in output.entries["10509"].children)
    assert(output.entries["129875"].parent == "10509")


def test_krakenreport_set_sibling_ranks():
    """Test KrakenReport."""
    input_prefix = "tests/data/taxid_630"
    input_report = f"{input_prefix}/Viral.kraken_report.txt"
    output = KrakenReport(input_report)
    output.set_sibling_ranks()

    ranks = {"2":1,"2759":2,"91347":1,"629":2, "0":0}
    for taxon_id,entry in output.entries.items():
        if taxon_id in ranks:
            assert(entry.sibling_rank == ranks[taxon_id])
        else:
            print(taxon_id, entry.sibling_rank)
            assert (entry.sibling_rank == 1)

