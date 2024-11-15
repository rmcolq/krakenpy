from collections import defaultdict
import sys

from krakenpy.report import KrakenReport
from krakenpy.assignment import KrakenAssignments

def merge_all_assignments(list_assignment_files, output_file):
    kraken_assignments = KrakenAssignments(output_file)
    changes = defaultdict(lambda: defaultdict(int))

    for assignment_file in list_assignment_files:
        changes = kraken_assignments.update(assignment_file, changes)

    kraken_assignments.save()
    return changes

def check_pair(kraken_assignment_file, kraken_report_file):
    report_stem = kraken_report_file.split("/")[-1].split("kraken")[0]
    assignment_stem = kraken_assignment_file.split("/")[-1].split("kraken")[0]
    if report_stem != assignment_stem:
        print(f"Found report stem {report_stem} and assignment stem {assignment_stem} from files {kraken_report_file} and {kraken_assignment_file}")
    assert(report_stem == assignment_stem)

    kreport = KrakenReport(kraken_report_file)
    kassignments = KrakenAssignments(kraken_assignment_file, load=True)

    counts = defaultdict(int)
    for read_id,entry in kassignments.entries.items():
        counts[entry.taxon_id] += 1
    #print(counts)

    for taxon_id in kreport.entries:
        if counts[taxon_id] != kreport.entries[taxon_id].ucount:
            print(f"A: Taxon id {taxon_id} has {kreport.entries[taxon_id].ucount} counts in report and {counts[taxon_id]} counts in assignment file")

        if taxon_id in counts:
            del counts[taxon_id]
    for taxon_id in counts:
        print(
            f"B: Taxon id {taxon_id} has {kreport.entries[taxon_id].ucount} counts in report and {counts[taxon_id]} counts in assignment file")


def merge(kraken_assignment_files, kraken_report_files, out_prefix):
    print("Initialize merged KrakenAssignments and KrakenReport")
    merged_assignments = KrakenAssignments(f"{out_prefix}.kraken_assignments.tsv")
    merged_reports = KrakenReport()

    assert(len(kraken_assignment_files) == len(kraken_report_files))

    pairs = zip(kraken_assignment_files, kraken_report_files)
    for assignment_file, report_file in pairs:
        print(f"Update with pair {assignment_file} and {report_file}")
        check_pair(assignment_file, report_file)

        changes = merged_assignments.update(assignment_file)
        merged_reports.update(report_file, changes)
        merged_assignments.save()
        merged_reports.save(f"{out_prefix}.kraken_report.txt")

    print(f"Save results to {out_prefix}.kraken_assignments.tsv and {out_prefix}.kraken_report.txt")
    merged_assignments.save()
    merged_reports.save(f"{out_prefix}.kraken_report.txt")