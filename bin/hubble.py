


# take in vcf?
# create embeddings
#   - use *1 as reference

#   - For each sample in the vcf
#   - Look up embeddings for all listed variants - need external file
#   - identify any variants without embeddings
#   - Get the index where that embedding belongs - need external file
#   - Replace the embedding at the location
#   - output seq
#   - pass through model
#   - have a force options that will ignore indels or something.  not needed in v1

"""
Greg McInes
Altman Lab
gmcinnes@stanford.edu
"""

from glib import hubble_welcome
hubble_welcome()

import os
import sys
import argparse
import numpy as np

from gseq import get_vcf_subject_ids, parse_vcf_line
from predict import HubblePredict

class Hubble(object):
    def __init__(self, vcf=None, seq=None, output=None, force=False, debug=False):
        self.vcf = vcf
        self.seq = seq
        self.output = output
        self.force = force
        self.debug = debug

        self.run()

    def run(self):
        if self.seq is not None:
            self.seq_run()
        elif self.vcf is not None:
            self.vcf_run()
        else:
            print("No input file selected!  Either a VCF or a seq file is required.", file=sys.stderr)
            exit(1)


    def seq_run(self):
        if self.debug:
            print("Running from seq file")
        hubble = HubblePredict(seq_file=self.seq, output=self.output, debug=self.debug)
        hubble.run()

        if self.debug:
            print("Complete!")

    def vcf_run(self):
        if self.debug:
            print("Running from VCF.  Mapping variants to pre-made embeddings.")

        # Generate the seq data
        seqs = self.build_seqs()

        # Create the data object
        data = self.format_seqs(seqs)

        # Predict
        hubble = HubblePredict(output=self.output, debug=self.debug)
        hubble.run(data=data)

    def build_seqs(self):
        # Build a dictionary of seqs for each sample using the reference
        samples = get_vcf_subject_ids(self.vcf)
        sample_seqs = self.sample_seq_init()
        embeddings = self.precomputed_embeddings()

        # Read in file variant by variant
        with open(self.vcf) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                vcf_row = parse_vcf_line(line)

                # Get the embedding index for the current variant
                ref_key = "%s_%s" % (vcf_row.pos, vcf_row.ref)
                alt_key = "%s_%s" % (vcf_row.pos, vcf_row.alt)

                ref_embedding = self.get_embedding(embeddings, ref_key)
                alt_embedding = self.get_embedding(embeddings, alt_key)

                if ref_embedding is None or alt_embedding is None:
                    print("Variant does not have a precomputed embedding!")
                    vcf_row.print_row(summary=True)
                    if self.force:
                        print("Skipping!")
                        continue
                    else:
                        print("Use vcf2seq.sh to create embeddings or --force to skip.")
                        exit(1)

                # Update with the appropriate embedding index
                #vcf_row.print_row(summary=True)
                for i in range(len(samples)):
                    s = samples[i]
                    gt = vcf_row.calls[i]

                    if gt == "0/0" or gt == "0|0":
                        position_index = ref_embedding['position_index']
                        embedding_index = ref_embedding['embedding_index']
                        #old = sample_seqs[s][position_index]
                        sample_seqs[s][position_index] = embedding_index
                        #new = sample_seqs[s][position_index]

                    elif gt == "1/1" or gt == "1|1":
                        position_index = alt_embedding['position_index']
                        embedding_index = alt_embedding['embedding_index']
                        #old = sample_seqs[s][position_index]
                        sample_seqs[s][position_index] = embedding_index
                        #new = sample_seqs[s][position_index]

                    else:
                        print("Unsupported genotype: %s" % gt)
                        print("Use homozygous reference or homozygous alt to indicate allele status of haplotypes.")
                        exit(1)

        return sample_seqs

    def format_seqs(self, seq_data):

        sample_names = []
        seqs = []

        # This value is the index for the null vector, i.e. no annotations at all
        null_vector = 2048

        for k in seq_data.keys():
            sample_names.append(k)
            seq = seq_data[k]
            null = [null_vector] * 32
            full_seq = seq.copy() + null + seq.copy()
            seqs.append(full_seq)

        # Convert the indices to embeddings
        X_ind = np.array(seqs)

        hub = HubblePredict()
        X = hub.indices2embeddings(X_ind)

        # Put the data in a dictionary and return
        data = {
            "X": X,
            "sample_names": sample_names
        }

        return data

    def get_embedding(self, embeddings, key):
        if key in embeddings:
            return embeddings[key]
        else:
            return None

    def sample_seq_init(self):
        samples = get_vcf_subject_ids(self.vcf)
        ref_seq = self.reference_seq()

        sample_dict = {}
        for s in samples:
            sample_dict[s] = ref_seq.copy()
        return sample_dict

    def reference_seq(self):
        wd = sys.path[0]
        file = os.path.join(wd, "../data/ref.seq")

        with open(file) as f:
            for line in f:
                fields = line.rstrip().split()
                seq = [int(x) for x in fields[1].split(',')]
        return seq

    def precomputed_embeddings(self):
        wd = sys.path[0]
        file = os.path.join(wd, "../data/embeddings.txt")

        embeddings = {}
        with open(file) as f:
            for line in f:
                fields = line.rstrip().split()
                key = "%s_%s" % (fields[1], fields[2])
                embeddings[key] = {"position": int(fields[1]),
                                   "allele": fields[2],
                                   "position_index": int(fields[3]),
                                   "embedding_index": int(fields[4])}

        return embeddings



"""
Parse the command line
"""
def parse_command_line():
    hubble_welcome()
    parser = argparse.ArgumentParser(
        description = 'Hubble.2D6 predicts CYP2D6 haplotype function. Provide a VCF with positions contained within '
                      'each haplotype or star allele.  Indicate allele status with a homozygous alternate allele (1/1). '
                      'Currently only hg19 is supported. See our manuscript for details on Hubble: '
                      'doi.org/10.1101/684357')

    parser.add_argument("-v", "--vcf",
                        help="VCF file input.  SNVs and common INDELs will be mapped to embeddings automatically. VCF "
                             "can be multisample for multiple haplotypes.  Haplotypes carrying a variant should "
                             "have a homozygous alternate allele in the appropriate column (1/1). Multiallelic sites not "
                             "yet supported, split multiallelic sites into different files.")
    parser.add_argument("-s", "--seq",
                        help="SEQ file input.  File containing embedding indices.  Used for annotation of uncommon INDELs.")
    parser.add_argument("-o", "--output", help="Output file. Default=hubble_output.txt")
    parser.add_argument("-f", "--force", action='store_true', default=False,
                        help="Ignore warnings for unmapped variants")
    parser.add_argument("-d", "--debug", action='store_true', default=False,
                        help="Output debugging messages.  May be very verbose.")
    options = parser.parse_args()

    return options


"""
Main
"""
if __name__ == "__main__":
    options = parse_command_line()
    Hubble(vcf=options.vcf,
           seq=options.seq,
           output=options.output,
           force=options.force,
           debug=options.debug)

