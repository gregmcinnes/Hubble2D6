"""
Greg McInes
Altman Lab
gmcinnes@stanford.edu
"""

import argparse
from gseq import fasta_extract, get_vcf_subject_ids, chr_greater_than_chr, parse_vcf_line, bed_extract, is_deletion
from glib import g_open
#from config import FASTA
import pandas as pd
import numpy as np
import itertools
import sys

CHROM = 0
START = 1
RSID = 2
REF = 3
ALT = 4
INFO = 7

class gVCF2seq(object):
    def __init__(self, vcf, zero_pad=10, one_hot=False, single_strand=False, max_pad=0,
                 annotations=None, generate=False, genosum=False, prefix="gvcf2seq", no_indels=True, unphase=False,
                 debug=False):
        self.vcf = vcf
        self.zero_pad = zero_pad
        self.one_hot = one_hot
        self.single_strand = single_strand
        self.max_pad = max_pad
        self.annotations = None # set in create_annotation_df
        self.generate = generate
        self.genosum = genosum
        self.prefix = prefix
        self.no_indels = no_indels
        self.unphase = unphase
        self.debug = debug

        self.annotation_df = self.create_annotation_df(annotations)
        self.annotation_lookup = self.annotation_df_dict()
        self.annotation_header = list(self.annotation_df)
        self.null_annotation = self.get_annotation_index('0')

        self.run()

    def create_annotation_df(self, columns):

        if self.genosum:
            genosum = [np.array([0]),
                       np.array([1]),
                       np.array([2])]

            full = []
            if columns is not None:
                user_specified = columns.split(",")
                self.annotations = user_specified
                annotation_columns = ['key', 'gt'] + user_specified
                K = len(user_specified)
                all_annotation_combos = [np.reshape(np.array(i), (1, K)) for i in itertools.product([0, 1], repeat=K * 1)]
                count = 0
                for row in genosum:
                    for combo in all_annotation_combos:
                        new_row = np.concatenate([[count], np.concatenate([row, combo[0]])])
                        full.append(new_row)
                        count += 1
            else:
                annotation_columns = ['key', 'gt']
                count = 0
                for row in genosum:
                    new_row = np.concatenate([[count], row])
                    full.append(new_row)
                    count += 1

        else:
            one_hot = [np.array([1, 0, 0, 0]),
                    np.array([0, 1, 0, 0]),
                    np.array([0, 0, 1, 0]),
                    np.array([0, 0, 0, 1]),
                    np.array([0, 0, 0, 0])]

            full = []
            if columns is not None:
                user_specified = columns.split(",")
                self.annotations = user_specified
                annotation_columns = ['key', 'A', 'C', 'T', 'G'] + user_specified
                K = len(user_specified)
                all_annotation_combos = [np.reshape(np.array(i), (1, K)) for i in itertools.product([0, 1], repeat=K * 1)]
                count = 0
                for row in one_hot:
                    for combo in all_annotation_combos:
                        new_row = np.concatenate([[count], np.concatenate([row, combo[0]])])
                        full.append(new_row)
                        count += 1
            else:
                annotation_columns = ['key', 'A', 'C', 'T', 'G']
                count = 0
                for row in one_hot:
                    new_row = np.concatenate([[count], row])
                    full.append(new_row)
                    count += 1

        df = pd.DataFrame.from_records(full, columns=annotation_columns)
        df.set_index('key', inplace=True)

        return df

    def process_subject(self, subject, index, generate=False):
        seq_a = self.initiate_seq()
        seq_b = self.initiate_seq()
        indices = {}

        annotation_counts = {}
        if self.annotations is not None:
            for a in self.annotations:
                annotation_counts[a] = 0

        with open(self.vcf) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                v = parse_vcf_line(line)

                if self.genosum and v.is_multiallelic():
                    continue

                if generate:
                    # generate data based on allele frequencies

                    a = None # create function that generates these
                    b = None

                    pass
                else:
                    genotype = v.calls[index].split(":")[0]
                    if "/" in genotype:
                        if self.unphase:
                            a, b = sorted(genotype.split("/"))
                        else:
                            a, b = genotype.split("/")
                    elif "|" in genotype:
                        if self.unphase:
                            a, b = sorted(genotype.split("|"))
                        else:
                            a, b = genotype.split("|")
                    else:
                        print("Unrecognized genotype delimiter! %s" % genotype)
                        exit()

                if v.alt.startswith("<"):
                    a = '0'
                    b = '0'

                allele_a = self.get_allele(a, v)
                allele_b = self.get_allele(b, v)

                allele_a_annotations = self.fetch_allele_annotation(allele_a, v)
                allele_b_annotations = self.fetch_allele_annotation(allele_b, v)

                genotype = v.calls[index].split(":")[0]
                split_string = "/"
                if genotype[1] == "|":
                    split_string = "|"
                a, b = genotype.split(split_string)

                if self.annotations is not None:
                    for anno in self.annotations:
                        if str(a) == "1":
                            annotation_counts[anno] += int(allele_a_annotations[anno])
                        if str(b) == "1":
                            annotation_counts[anno] += int(allele_b_annotations[anno])

                if self.genosum:
                    if a == ".":
                        a = 0
                    if b == ".":
                        b = 0
                    genosum = int(a) + int(b)
                    annotations = allele_a_annotations
                    for a in allele_b_annotations:
                        annotations[a] = max(int(allele_a_annotations[a]), int(allele_b_annotations[a]))

                    seq_a = self.update_seq(seq_a, genosum, annotations)
                else:
                    buffer_a = 0
                    buffer_b = 0

                    if self.no_indels:
                        # If we're not outputing the actual sequence then just output the first nucleotide from
                        # the observed allele.  So if it's an indel and the length of the ref or alt is more than 1,
                        # just output the first base.  No need to track buffers and shit like that.
                        allele_a = allele_a[0]
                        allele_b = allele_b[0]

                    else:

                        # Check if the current position has the possibility of an insertion, if so set the buffer to the
                        # difference between the length of the allele present and the length of the largest possible
                        # insertion
                        if v.max_insertion() > len(allele_a):
                            buffer_a = v.max_insertion() - len(allele_a) #+ 1 #+ len(v.ref)

                        if v.max_insertion() > len(allele_b):
                            buffer_b = v.max_insertion() - len(allele_b) #+ 1

                        # Check if the current variant is a deletion, if so set the buffer to the difference between the
                        # lengths of the alternate and reference alleles

                        if is_deletion(v.ref, allele_a):
                            offset=max(len(v.ref), v.max_insertion())
                            buffer_a = offset - len(allele_a)

                        if is_deletion(v.ref, allele_b):
                            offset = max(len(v.ref), v.max_insertion())
                            buffer_b = offset - len(allele_b)

                        if (len(allele_a) + buffer_a) != (len(allele_b) + buffer_b):
                            print("Unequal sequence additions!", file=sys.stderr)
                            print("Ref: %s" % v.ref, file=sys.stderr)
                            print("Allele A: %s, Buffer A: %s, Sum: %s" % (allele_a, buffer_a, (len(allele_a) + buffer_a)), file=sys.stderr)
                            print("Allele B: %s, Buffer B: %s, Sum: %s" % (allele_b, buffer_b, (len(allele_b) + buffer_b)), file=sys.stderr)
                            print("Max insertion: %s" % v.max_insertion(), file=sys.stderr)
                            print("Is deletion: %s, %s" % (is_deletion(v.ref, allele_a), is_deletion(v.ref, allele_b)), file=sys.stderr)
                            print(v.alts(), file=sys.stderr)
                            v.max_insertion2()

                    seq_a = self.update_seq(seq_a, allele_a, allele_a_annotations, buffer_a)
                    seq_b = self.update_seq(seq_b, allele_b, allele_b_annotations, buffer_b)

                    indices[v.pos] = len(seq_a) - 1
        f.close()

        if len(seq_a) != len(seq_b):
            print("UNEQUAL SEQUENCE LENGTHS!", file=sys.stderr)
            print(subject, file=sys.stderr)
            return

        if self.genosum:
            print(subject + "\t" + "\t".join(str(x) for x in seq_a))
        else:
            # max pad if necessary
            if self.max_pad != 0:
                seq_a = self.pad_to_max(seq_a)
                seq_b = self.pad_to_max(seq_b)

            # Merge the sequences together
            seq_a_out = ",".join(str(x) for x in seq_a)
            seq_b_out = ",".join(str(x) for x in seq_b)
            print(subject + "\t" + seq_a_out + "\t" + seq_b_out)
        return annotation_counts, indices

    def pad_to_max(self, seq):
        if self.debug:
            print("Max padding sequences")
        if len(seq) >= self.max_pad:
            if self.debug:
                print("Seq A too long, %s.  Cutting!" % len(seq))
            seq = seq[0:self.max_pad]
        else:
            if self.debug: print("Seq A too short, %s.  Adding!" % len(seq))
            while len(seq) < self.max_pad:
                seq.append(self.null_annotation)
        return seq

    def update_seq(self, seq, allele, annotations, buffer=0):
        # a buffer is used when a sequence undergoes a deletion or when another sequence has an insertion
        if self.genosum:
            i = self.get_annotation_index(allele, annotations)
            seq.append(i)
        else:
            for a in allele:
                # get the index
                i = self.get_annotation_index(a, annotations)
                seq.append(i)
            for i in range(buffer):
                seq.append(self.null_annotation)
        return seq

    def get_allele(self, index, vcf):
        if index == '0':
            return vcf.ref
        if index == '.':
            return '0'
        return vcf.alts()[int(index)-1]

    def fetch_allele_annotation(self, allele, vcf):
        # Get annotations
        annotations = {}
        if self.annotations is None:
            return annotations
        for a in self.annotations:
            key = "%s_%s" % (a, allele)
            if key in vcf.info:
                annotations[a] = vcf.info[key]
            else:
                annotations[a] = 0
        return annotations

    def initiate_seq(self):
        seq = []
        for x in range(self.zero_pad):
            seq.append(self.null_annotation)
        return seq

    def get_annotation_index(self, nucleotide, annotations=None):
        # Get zero vector
        if nucleotide == '0':
            if self.genosum:
                query = "gt == 0"
            else:
                query = "A == 0 and G == 0 and T == 0 and C == 0"
            if self.annotations is not None:
                for a in self.annotations:
                    query += " and %s == 0" % a
            result = self.annotation_df.query(query)
            nrow = result.shape[0]
            if nrow != 1:
                print("Invalid query results")
                print(query)
                exit()
            index = result.index.tolist()[0]

            return(index)

        adf_key_list = []
        nucleotides = ["A", "C", "T", "G"]
        for n in nucleotides:
            if n == nucleotide:
                adf_key_list.append(1)
            else:
                adf_key_list.append(0)

        for a in self.annotation_header:
            if a in nucleotides:
                continue
            adf_key_list.append(annotations[a])

        adf_key = "".join(str(x) for x in adf_key_list)

        dict_index = self.annotation_lookup[adf_key]

        return(dict_index)

    def annotation_df_dict(self):
        annotation_df_list = self.annotation_df.values.tolist()
        annotation_dict = {}
        for i in range(len(annotation_df_list)):
            list_str = str("".join(str(x) for x in annotation_df_list[i]))
            annotation_dict[list_str] = i
        return(annotation_dict)

    def run(self):
        subjects = get_vcf_subject_ids(self.vcf)

        annotation_counts = {}

        for s in range(0, len(subjects)):
            subject = subjects[s]

            if self.debug:
                print("NEW SUBJECT")
                print(subject)
                print(s)

            annotation_counts[subject], index = self.process_subject(subject, s, generate=self.generate)

        self.annotation_df.to_csv("%s.annotation_embeddings.csv" % self.prefix)

        if self.annotations is not None:
            file = open("%s.annotation_counts.txt" % self.prefix, "w")
            header = "subject," + ",".join(self.annotations) + "\n"
            file.write(header)
            for s in annotation_counts:
                output = [s]
                for a in self.annotations:
                    output.append(annotation_counts[s][a])
                file.write(",".join(str(x) for x in output) + "\n")
            file.close()

        file = open("%s.index.txt" % self.prefix, "w")
        for i in index.keys():
            file.write("%s,%s\n" % (i, index[i]))

"""
Parse the command line
"""
def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'VCF to gVCF script.  Fills in reference calls for a VCF')
    parser.add_argument("--vcf", help="Input vcf")
    parser.add_argument("--zero_pad_length", type=int, default=10, help="Zero pad between bed regions")
    parser.add_argument("--one_hot", action="store_true", help="One hot encode genotypes")
    parser.add_argument("--genosum", action="store_true")
    parser.add_argument("--no_indels", action="store_true")
    parser.add_argument("--max_pad", type=int, default=0)
    parser.add_argument("--annotations", default=None)
    parser.add_argument("--unphase", action="store_true")
    parser.add_argument("--prefix", default="gvcf2seq")
    parser.add_argument("-d", "--debug", action='store_true', default=False,
                                help="Output debugging messages.  May be very verbose.")
    options = parser.parse_args()
    return options

"""
Main
"""
if __name__ == "__main__":
    options = parse_command_line()
    gVCF2seq(vcf=options.vcf,
             zero_pad=options.zero_pad_length,
             one_hot=options.one_hot,
             max_pad=options.max_pad,
             annotations=options.annotations,
             genosum=options.genosum,
             prefix=options.prefix,
             no_indels=options.no_indels,
             unphase=options.unphase,
             debug=options.debug)

