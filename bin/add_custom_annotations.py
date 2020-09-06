"""
Greg McInes
Altman Lab
gmcinnes@stanford.edu
"""

from gseq import fasta_extract, bed_extract, parse_vcf_line, chr_greater_than_chr, get_vcf_subject_ids, is_indel, \
    is_rare, is_deleterious, is_exonic, high_at, high_gc, is_tf_binding_site, is_deleterious_2, \
    is_dnase_hypersensitivity_site, is_eQTL, is_methylation_site, variant_class, get_gene, is_high_impact, get_impact
import argparse



class AddCustomAnnotations(object):
    def __init__(self, vcf, annotations=False, fasta=None, debug=False):
        self.vcf = vcf
        self.user_annotations = annotations
        self.fasta = fasta
        self.debug = debug

        self.run()

    def run(self):
        annotations = self.set_annotations()

        with open(self.vcf) as f:
            for line in f:
                if line.startswith("#"):
                    print(line.rstrip())
                    continue
                v = parse_vcf_line(line, annovar=True)
                v = self.add_annotations(v, annotations)
                v.print_row()


    def print_annotations(self, vcf_row, annotations):
        ref_annotations = []
        alt_annotations = []

        ref = vcf_row.ref
        alt = vcf_row.alt
        pos = vcf_row.pos

        for a in annotations:

            if "%s_%s" % (a, ref) in vcf_row.info:
                ref_annotations.append(vcf_row.info["%s_%s" % (a, ref)])
            else:
                ref_annotations.append(0)

            if "%s_%s" % (a, alt) in vcf_row.info:
                alt_annotations.append(vcf_row.info["%s_%s" % (a, alt)])
            else:
                alt_annotations.append(0)

        ref_out = [pos, ref] + ref_annotations
        alt_out = [pos, alt] + alt_annotations

        print(",".join(str(x) for x in ref_out))
        print(",".join(str(x) for x in alt_out))



    def add_annotations(self, vcf_row, annotations):
        if "rare" in annotations:
            if "gnomAD_genome_ALL" in vcf_row.info:
                af = vcf_row.info["gnomAD_genome_ALL"]
                vcf_row.add_info("rare_%s" % vcf_row.alt, 1 if is_rare(af) else 0)
            else:
                if self.debug: print("Allele frequency missing!")
                vcf_row.add_info("rare_%s" % vcf_row.alt, 0)

        if "deleterious" in annotations:
            vcf_row.add_info("deleterious_%s" % vcf_row.alt, 1 if is_deleterious_2(vcf_row.info, vcf_row.ref, vcf_row.alt) else 0)

        if "coding"  in annotations:
            vcf_row.add_info("coding_%s" % vcf_row.alt, 1 if is_exonic(vcf_row.info) else 0)
            vcf_row.add_info("coding_%s" % vcf_row.ref, 1 if is_exonic(vcf_row.info) else 0)

        if "indel" in annotations:
            vcf_row.add_info("indel_%s" % vcf_row.alt, 1 if is_indel(vcf_row.ref, vcf_row.alt) else 0)

        if "high_gc" in annotations:
            vcf_row.add_info("gc_%s" % vcf_row.alt, 1 if high_gc(self.fasta, vcf_row.chrom, vcf_row.pos-10, vcf_row.pos+10) else 0)
            vcf_row.add_info("gc_%s" % vcf_row.ref,
                             1 if high_gc(self.fasta, vcf_row.chrom, vcf_row.pos - 10, vcf_row.pos + 10) else 0)

        if "high_at" in annotations:
            vcf_row.add_info("at_%s" % vcf_row.alt, 1 if high_at(self.fasta, vcf_row.chrom, vcf_row.pos-10, vcf_row.pos+10) else 0)
            vcf_row.add_info("at_%s" % vcf_row.ref,
                             1 if high_at(self.fasta, vcf_row.chrom, vcf_row.pos - 10, vcf_row.pos + 10) else 0)

        if "active_site" in annotations:
            if "cyp2d6_active_site" in vcf_row.info:
                site = vcf_row.info["cyp2d6_active_site"]
                if site != ".":
                    vcf_row.add_info("active_site_%s" % vcf_row.alt, 1)
                    vcf_row.add_info("active_site_%s" % vcf_row.ref, 1)
                else:
                    vcf_row.add_info("active_site_%s" % vcf_row.alt, 0)
                    vcf_row.add_info("active_site_%s" % vcf_row.ref, 0)

        if "tfbs" in annotations:
            vcf_row.add_info("tfbs_%s" % vcf_row.ref, 1 if is_tf_binding_site(vcf_row.info) else 0)
            vcf_row.add_info("tfbs_%s" % vcf_row.alt, 1 if is_tf_binding_site(vcf_row.info) else 0)

        if "eqtl" in annotations:
            vcf_row.add_info("eqtl_%s" % vcf_row.alt, 1 if is_eQTL(vcf_row.info) else 0)

        if "eqtlpos" in annotations:
            vcf_row.add_info("eqtlpos_%s" % vcf_row.alt, 1 if is_eQTLpos(vcf_row.info) else 0)

        if "eqtlneg" in annotations:
            vcf_row.add_info("eqtlneg_%s" % vcf_row.alt, 1 if is_eQTLneg(vcf_row.info) else 0)

        if "dnase" in annotations:
            vcf_row.add_info("dnase_%s" % vcf_row.ref, 1 if is_dnase_hypersensitivity_site(vcf_row.info) else 0)
            vcf_row.add_info("dnase_%s" % vcf_row.alt, 1 if is_dnase_hypersensitivity_site(vcf_row.info) else 0)

        if "methyl" in annotations:
            vcf_row.add_info("methyl_%s" % vcf_row.ref, 1 if is_methylation_site(vcf_row.info) else 0)
            vcf_row.add_info("methyl_%s" % vcf_row.alt, 1 if is_methylation_site(vcf_row.info) else 0)

        if "gene" in annotations:
            gene = get_gene(vcf_row.info)
            vcf_row.add_info("gene", gene)

        if "impact" in annotations:
            impact = get_impact(vcf_row.info)
            vcf_row.add_info("impact", impact)

        if "high_impact" in annotations:
            vcf_row.add_info("highimpact", 1 if is_high_impact(vcf_row.info) else 0)

        return vcf_row

    def set_annotations(self):
        if self.user_annotations is not None:
            annotations = self.user_annotations.split(",")
        else:
            annotations = ['deleterious', 'rare'] # add more here as you come up with them
        return annotations
"""
Parse the command line
"""
def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'VCF to gVCF script.  Fills in reference calls for a VCF')
    parser.add_argument("--vcf", help="Input vcf")
    parser.add_argument("--annotations", default=None)
    parser.add_argument("--output", default=None)
    parser.add_argument("--fasta", default=None)
    parser.add_argument("-d", "--debug", action='store_true', default=False,
                                help="Output debugging messages.  May be very verbose.")
    options = parser.parse_args()
    return options


"""
Main
"""
if __name__ == "__main__":
    options = parse_command_line()
    AddCustomAnnotations(options.vcf, options.annotations, options.fasta, options.debug)

