"""
Greg McInes
Altman Lab
gmcinnes@stanford.edu
"""

from gseq import fasta_extract, bed_extract, parse_vcf_line, chr_greater_than_chr, get_vcf_subject_ids, is_indel
import argparse

class MergeMultiallelic(object):
    def __init__(self, vcf, debug=False):
        self.debug = debug
        self.vcf = vcf
        self.debug = debug
        self.run()

    def run(self):

        vcf = open(self.vcf)
        next = vcf.readline()
        last = None

        while next:
            if next.startswith("#"):
                print(next.rstrip())
                next = vcf.readline()
                continue

            v = parse_vcf_line(next)

            if last and last.pos == v.pos and last.ref == v.ref:
                # Last becomes a combination of the two
                last.make_multiallelic(v)
            else:
                if last:
                    last.print_row()
                last = v
            next = vcf.readline()

"""
Parse the command line
"""
def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'This is a script I wrote')
    parser.add_argument("-f", "--vcf", help="VCF file to expand")
    parser.add_argument("-d", "--debug", action='store_true', default=False,
                                help="Output debugging messages.  May be very verbose.")
    options = parser.parse_args()
    return options


"""
Main
"""
if __name__ == "__main__":
    options = parse_command_line()
    MergeMultiallelic(options.vcf, options.debug)

