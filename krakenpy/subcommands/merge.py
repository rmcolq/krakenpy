from krakenpy.merge import *

def run(options):

    merge(options.in_assignments,
        options.in_reports,
        options.out_prefix
        )