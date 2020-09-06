"""
Greg McInes
Altman Lab
gmcinnes@stanford.edu
"""

from gseq import fasta_extract, bed_extract, parse_vcf_line, chr_greater_than_chr, get_vcf_subject_ids, is_indel
import argparse

class VCF2gVCF(object):
    def __init__(self, vcf, bed, fasta, format=None, buffer_deletions=False, debug=False):
        self.debug = debug
        self.vcf = vcf
        self.bed = bed
        self.fasta = fasta
        self.format = format
        if self.format is None:
            self.format = "GT"
        self.buffer_deletions = buffer_deletions
        self.debug = debug
        self.run()

    def run(self):
        subjects = get_vcf_subject_ids(self.vcf)
        bed = bed_extract(self.bed)
        vcf = open(self.vcf)
        next_vcf_line = self.get_next_vcf_line(vcf)
        # Get the first line of the vcf that's not a comment
        while next_vcf_line is None:
            next_vcf_line = self.get_next_vcf_line(vcf)

        for b in range(bed.count):
            bed_chrom, bed_start, bed_end = bed.retrieve_entry(b)
            # Check that the current vcfline is in range
            chrom, start, end = bed.retrieve_entry(b)
            start = int(start)
            end = int(end)
            ref_seq = fasta_extract(self.fasta, chrom, start, end)

            next_vcf_line = self.update_vcf_row(vcf, bed_chrom, bed_start, next_vcf_line)

            skip = 0
            for i in range(start, end):
                if skip > 0:
                    skip -= 1
                    continue

                if next_vcf_line is not None and i == next_vcf_line.pos:
                    if self.buffer_deletions and len(next_vcf_line.ref) > 1:
                        skip = len(next_vcf_line.ref) - 1

                    next_vcf_line.print_row()

                    next_vcf_line = self.update_vcf_row(vcf, bed_chrom, i+1, next_vcf_line)
                else:
                    gt = self.format_genotype()
                    genotypes = "\t".join([gt] * len(subjects))

                    print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (bed_chrom, i, ".", ref_seq[i-start].upper(), ".", "30",
                                                                  "PASS", '.', self.format, genotypes))
                    pass
                next_vcf_line = self.update_vcf_row(vcf, bed_chrom, i + 1, next_vcf_line)

    def format_genotype(self):
        fields = self.format.split(":")
        if len(fields) == 1:
            gt = "0/0"
        else:
            extra = ":".join(["."] * (len(fields) - 1))
            gt = "0/0:%s" % (extra)
        return gt

    def get_next_vcf_line(self, vcf):
        next = vcf.readline()
        if next.startswith("#"):
            print(next.rstrip())
            return None
        if len(next) == 0:
            return None
        #print(next)
        next = parse_vcf_line(next)
        try:
            next.chrom = int(next.chrom)
        except:
            if not next:
                return ['chrM', 99999999999]
            return next
        return next

    def update_vcf_row(self, vcf_object, bed_chrom, bed_pos, current=None):
        if current and self.check_range(bed_chrom, current.chrom, bed_pos, current.pos):
            return current
        complete = False
        while not complete:
            next_row = self.get_next_vcf_line(vcf_object)
            if next_row is None:
                return next_row
            complete = self.check_range(bed_chrom, next_row.chrom, bed_pos, next_row.pos)
        return next_row

    def check_range(self, bed_chrom, vcf_chrom, bed_pos, vcf_pos):
        if chr_greater_than_chr(bed_chrom, vcf_chrom):
            return False
        if bed_pos > vcf_pos:
            return False
        return True


"""
Parse the command line
"""
def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'This is a script I wrote')
    parser.add_argument("-f", "--vcf", help="VCF file to expand")
    parser.add_argument("-b", "--bed", help="Bed file with coordinates")
    parser.add_argument("-F", "--fasta", help="Fasta file to pull reference from")
    parser.add_argument("--format", default=None, help="format field")
    parser.add_argument("--buffer_deletions", action='store_true', default=False, help="Remove lines following a deletion")
    parser.add_argument("-d", "--debug", action='store_true', default=False,
                                help="Output debugging messages.  May be very verbose.")
    options = parser.parse_args()
    return options


"""
Main
"""
if __name__ == "__main__":
    options = parse_command_line()
    VCF2gVCF(options.vcf, options.bed, options.fasta, options.format, options.buffer_deletions, options.debug)

