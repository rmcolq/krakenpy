from collections import defaultdict
import csv
import sys


class KrakenEntry:
    """
    A class representing a line in a kraken report.

    Attributes:
        taxon_id (str): The NCBI taxon identifier.
        name (string): The scientific name associated with this taxon.
        rank (str): A letter coding the rank of this taxon.
        depth (int): The number of indentations this entry had in the kraken file (related to hierarchy).
        count (int): The count of reads assigned to this taxon and its descendants.
        name (int): The count of reads assigned specifically to this taxon.
        domain (str): The domain this taxon is a member of.
        parent (str): The taxon id associated with the taxonomic parent.
        children (set): A set of taxon ids associated with the direct taxonomic children.
        sibling_rank (int): An integer representing the ranking among direct siblings (share the parent) based on count.
        hierarchy (list): An ordered list of taxon ids representing the parents to taxonomic root.
    """
    def __init__(self, row=None, domain=None, hierarchy=[]):
        self.taxon_id = "0"
        self.name = "unclassified"
        self.rank = "U"
        self.depth = 0
        self.count = 0  # inclusive count
        self.ucount = 0  # unique_count
        self.domain = domain
        self.parent = None
        self.children = set()
        self.sibling_rank = 0
        self.hierarchy = hierarchy
        if row is not None:
            self.add_row(row)

    def print(self):
        print(
            f"{self.taxon_id},{self.name},{self.rank},{self.depth},{self.count},{self.ucount},{self.domain},{self.parent},{self.children},{self.sibling_rank},{self.hierarchy}")

    def parse_depth(self, name):
        parse_name = name.split(" ")
        depth = 0
        for i in parse_name:
            if i != "":
                break
            depth += 1
        depth = int(depth / 2)
        return depth

    def add_row(self, row):
        self.taxon_id = row["Taxonomy ID"]
        self.name = row["Scientific Name"].strip()
        self.depth = self.parse_depth(row["Scientific Name"])
        self.rank = row["Rank"]
        self.count = int(row["Clades"])  # inclusive count
        self.ucount = int(row["Taxonomies"])  # unique_count
        if self.count < self.ucount:
            self.count, self.ucount = self.ucount, self.count
        self.hierarchy = self.hierarchy[:self.depth]

    def add_parent(self, parent):
        self.parent = parent

    def add_child(self, child):
        self.children.add(child)

    def set_sibling_rank(self, rank):
        self.sibling_rank = rank


class KrakenReport:
    """
    A class representing a kraken report.

    Attributes:
        entries (dict): A dict with keys for taxon ids and values for KrakenEntry.
        total (int): Total number of reads in report.
        unclassified (int): Number of unclassified reads.
        classified (int): Number of classified reads.
        domains (int): A dict with keys for names of domains and values for associated taxon id.
    """
    def __init__(self, file_name=None):
        self.entries = defaultdict(KrakenEntry)
        self.total = 0
        self.unclassified = 0
        self.classified = 0
        self.domains = defaultdict(int)
        if file_name:
            self.load_df(file_name)
            self.unclassified = self.entries[0].count
            self.classified = self.entries[1].count
            self.total = self.classified + self.unclassified

    def print(self):
        print(f"Report has {len(self.entries)} taxon entries corresponding to {self.classified} classified and {self.unclassified} unclassified reads.")

    def add_parent_child(self, parent_id, child_id):
        self.entries[child_id].add_parent(parent_id)
        self.entries[parent_id].add_child(child_id)

    def set_sibling_ranks(self):
        for entry_id, entry in self.entries.items():
            if entry.sibling_rank > 0 or entry.parent in [None, 1, 131567]:
                continue
            if entry.rank in ["D", "R", "R1"]:
                entry.set_sibling_rank(1)
            elif len(self.entries[entry.parent].children) == 1:
                entry.set_sibling_rank(1)
            else:
                sibling_dict = {i: self.entries[i].count for i in self.entries[entry.parent].children}
                sorted_counts = sorted(sibling_dict.values(), reverse=True)
                for i, c in sibling_dict.items():
                    rank = sorted_counts.index(c) + 1
                    self.entries[i].set_sibling_rank(rank)

    def check_sibling_ranks(self):
        for entry_id in self.entries:
            if self.entries[entry_id].sibling_rank == 0:
                print(entry_id)

    def check_report(self, file_name):
        with open(file_name, 'r') as handle:
            line = handle.readline()
            num_fields = len(line.split("\t"))
            if num_fields not in [6,8]:
                sys.stderr.write(
                    f"Kraken report file {file_name} badly formatted - must have 6 or 8 columns"
                    )
                sys.exit(2)
            if line.startswith("%"):
                return True, num_fields
            else:
                return False, num_fields

    def load_df(self, file_name):
        csvfile = open(file_name, newline='')
        df = {}
        report_has_header, num_fields = self.check_report((file_name))
        if report_has_header:
            df = csv.DictReader(csvfile, delimiter="\t")
        elif num_fields == 6:
            df = csv.DictReader(csvfile, delimiter="\t", fieldnames=["% of Seqs","Clades","Taxonomies","Rank","Taxonomy ID","Scientific Name"])
        elif num_fields == 8:
            df = csv.DictReader(csvfile, delimiter="\t", fieldnames=["% of Seqs","Clades","Taxonomies","Read Minimizers", "Taxon Minimizers", "Rank","Taxonomy ID","Scientific Name"])

        hierarchy = []
        domain = None
        for row in df:
            try:
                if row["Rank"] == "D":
                    domain = row["Scientific Name"].strip()
                    self.domains[domain] = row["Taxonomy ID"]
                entry = KrakenEntry(row=row, domain=domain, hierarchy=hierarchy)

            except:
                print(f"Found badly formatted row:\n{row}")
                print(f"Quitting load of {file_name}")
                break

            self.entries[entry.taxon_id] = entry
            hierarchy = entry.hierarchy.copy()
            if len(hierarchy) > 0:
                self.add_parent_child(hierarchy[-1], entry.taxon_id)
            if entry.taxon_id != "0":
                hierarchy.append(entry.taxon_id)
        csvfile.close()
        self.set_sibling_ranks()
        # self.check_sibling_ranks()

    def get_domains(self):
        domains = []
        for entry_id, entry in self.entries.items():
            if entry.rank == "D":
                domains.append(entry)
                entry.print()
        return domains

    def get_tips(self):
        tips = []
        for entry_id, entry in self.entries.items():
            if len(entry.children) == 0:
                tips.append(entry)
                entry.print()
        return tips

    def get_rank_entries(self, rank):
        subset = []
        for entry_id, entry in self.entries.items():
            if entry.rank == rank:
                subset.append(entry)
                entry.print()
        return subset

    def get_percentage(self, taxon_id, domain=None):

        if domain and self.entries[taxon_id].domain != domain:
            return 0.0

        denominator = self.classified
        if domain:
            denominator = self.entries[self.domains[domain]].count

        count = self.entries[taxon_id].count
        return float(count) / denominator

    def to_source_target_df(self, max_rank=None, domain=None, trace_ids=[]):
        records = []
        ignore = set()
        skip = set()
        for entry_id, entry in self.entries.items():

            # we don't want the connections to cellular organisms, root etc
            if entry.sibling_rank == 0:
                continue

            # filter by domain where required
            if domain and entry.domain != domain:
                continue

            # filter by rank when specified
            if max_rank:
                if entry.sibling_rank > max_rank:
                    ignore.add(entry_id)
                    if entry_id in trace_ids:
                        print(entry_id, 1)
                    continue
                elif entry.parent in ignore:
                    ignore.add(entry_id)
                    if entry_id in trace_ids:
                        print(entry_id, 2)
                    continue

            # filter if an intermediate rank
            if entry.rank not in ["K", "D", "D1", "D2", "P", "C", "O", "F", "G", "S", "S1", "S2"]:
                skip.add(entry_id)
                if entry_id in trace_ids:
                    print(entry_id, 3)
                continue

            index = 1
            while index < len(entry.hierarchy) and entry.hierarchy[-index] in skip:
                index += 1
            source_id = entry.hierarchy[-index]
            if entry_id in trace_ids:
                print(entry_id, 4, source_id)
            records.append({"source": self.entries[source_id].name, "target": entry.name, "value": entry.count,
                            "percentage": self.get_percentage(entry_id, domain=domain)})

        with open("source_target.csv", 'w', newline='') as csvfile:
            fieldnames = ["source", "target", "value", "percentage"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for row in records:
                writer.writerow(row)

        print(len(records), len(ignore), len(skip))

        return records

    def to_df(self, sample_id="sample_id", ranks=[]):
        if not ranks or len(ranks) == 0:
            taxon_ids = [e for e in self.entries.keys()]
        else:
            taxon_ids = [e for e in self.entries.keys() if self.entries[e].rank in ranks]
        return pd.DataFrame({sample_id: [self.entries[e].count for e in taxon_ids]}, index=taxon_ids)

    def check_host(self, host_dict):
        for host_id, max_host_count in host_dict.items():
            if host_id in self.entries and self.entries[host_id].count > max_host_count:
                sys.stderr.write(
                    f"ERROR: found {self.entries[host_id].count} reads corresponding to host {self.entries[host_id].name} with taxon_id {host_id}, max allowed is {max_host_count}\n"
                    )
                sys.exit(2)