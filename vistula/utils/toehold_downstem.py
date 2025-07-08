import argparse
import subprocess
import math
from Bio import SeqIO
from Bio.Seq import Seq
import re, os
import shutil

KOZAK = "GCCACCAUGG"
LINKER = "AACCUGGCGGCAGCGCAAAAG"
S = ["G", "C"]
W = ["A", "U"]

parser = argparse.ArgumentParser(
    prog="SwitchForge",
    description="""Designs toehold switches based on the reporter and transcript
    sequences provided by user""",
    epilog="Made by HepaSwitch iGem 2025 team"
)

parser.add_argument('transcript_file')
parser.add_argument('reporter_file')

args = parser.parse_args()

transcript_file = args.transcript_file
reporter_file = args.reporter_file

# PREDICT THE SECONDARY STRUCTURE OF REPORTER

# result = subprocess.run(
#     ["RNAfold", "--noPS", reporter_file],  # --noPS pomija generowanie pliku .ps
#     check=True,
#     capture_output=True,
#     text=True
# )

# rep_sec_str = result.stdout.split("\n")[2].split(" ")[0]

# RUN DESI RNA ON DIFFERENT INPUTS


tmp_dir = "rnaup_tmp"

if not os.path.exists(tmp_dir):
    os.makedirs(tmp_dir)

def split_toehold(kmer_lengths):
    """
    For each k-mer length, splits it into two parts:
    - `part1`: 5' region that lies outside the hairpin — upstream of the structured region (not base-paired)
    - `part2`: 3' region that becomes the base-paired lower stem (ds), including the start codon (AUG)

    This division allows selecting the portion of the trigger sequence that will participate
    in forming the lower hairpin structure of the toehold switch.

    Returns a list of tuples (part1, part2) for each k-mer length.

    Note:
        This does NOT split the toehold structure itself.
        It partitions the input k-mer (target sequence) into:
            - a non-hybridizing region (ignored),
            - and a functional region (used to build the switch).
    """
    pair_lengths = []
    for n in kmer_lengths:
        part1 = math.floor(2 * n / 5)
        part2 = n - part1 
        pair_lengths.append((part1, part2))
    return pair_lengths

def split_with_fixed_middle(n):
    """
    Given a sequence length `n` (including a fixed 3-nt AUG codon in the center),
    returns the index ranges (slices) for each of the three parts:
    - upstream: before AUG
    - AUG: fixed start codon (3 nt)
    - downstream: after AUG

    This is useful when slicing a sequence directly, e.g.:
        seq[upstream[0]:upstream[1]], seq[aug[0]:aug[1]], seq[downstream[0]:downstream[1]]

    Parameters:
        n (int): Total length of the sequence, including AUG.

    Returns:
        tuple: Three tuples of indices: (upstream_range, aug_range, downstream_range)
    """
    remainder = n - 3
    n1 = remainder // 2
    
    upstream = (0, n1)
    aug = (n1, n1 + 3)
    downstream = (n1 + 3, n)

    return upstream, aug, downstream

kmer_lengths = [i for i in range(20, 36)]
pair_lengths = split_toehold(kmer_lengths) # lengths of 5' toehold part and the stem part

def update_toehold(toehold_seq, toehold_struct, seq, struct):
    """
    Append a nucleotide sequence and its corresponding secondary structure symbols
    to the running toehold sequence and structure strings.

    Parameters:
        toehold_seq (str): Current concatenated sequence of the toehold design.
        toehold_struct (str): Current concatenated dot-bracket notation representing structure.
        seq (str): New nucleotide sequence to append.
        struct (str): Corresponding dot-bracket symbols to append for seq.

    Returns:
        tuple: Updated (toehold_seq, toehold_struct) after appending seq and struct.
    """
    toehold_seq += seq
    toehold_struct += struct
    return toehold_seq, toehold_struct


def check_anti_aug(anti_aug):
    """
    Checks whether the anti-AUG region violates specific nucleotide constraints
    derived from forward-engineering rules (Green et al., 2014).

    The rules being checked:
    - The first nucleotide in anti_aug should NOT be 'C'
    - The second nucleotide should NOT be 'A' or 'G'
    - The third nucleotide should NOT be 'U'

    Parameters:
        anti_aug (str): A 3-nucleotide sequence complementary to the start codon region.

    Returns:
        bool: True if any constraint is violated (sequence is disallowed),
              False otherwise (sequence passes the check).
    """
    if anti_aug[0] == "C" or anti_aug[1] == "A" or anti_aug[1] == "G" or anti_aug[2] == "U":
        return True
    return False

reporter_seq = list(SeqIO.parse(reporter_file, "fasta"))[0].seq.replace("T", "U")
transcript_seq = list(SeqIO.parse(transcript_file, "fasta"))[0].seq.replace("T", "U")
best_score = 0
toehold = ''

for i in range(len(kmer_lengths)):

    klen = kmer_lengths[i] # k-mer length
    comp, stem = pair_lengths[i] # complementary-to-transcript part and stem part (legnths)

    stem_partition = split_with_fixed_middle(stem) # partition of toehold's stem part into pre-NON-AUG; NON-AUG; post-NON-AUG

    transcript_kmers = [transcript_seq[i:i+klen] for i in range(0, len(transcript_seq) - len(transcript_seq) % klen, klen)] # k-mers of length klen to check

    # complex
    for kmer in transcript_kmers:
        tseq = ""
        tstruct = ""

        toehold = Seq(kmer).reverse_complement().replace("T", "U") # toehold sequence

        th_1, th_2 = toehold[:comp], toehold[comp:] 

        # check Green et. al (2014) forward-engineering rules 
        if th_2[0] not in W or th_2[1] not in S or th_2[2] not in S:
            continue

        up, mid, down = stem_partition

        upstream_aug = th_2[up[0]:up[1]]
        anti_aug = th_2[mid[0]:mid[1]] 
        downstream_aug = "n" * len(th_2[down[0]:down[1]])
        
        comp_to_kmer = th_1 + upstream_aug + anti_aug

        # next_anti_aug = th_2[mid[0]+1:mid[1]+1] 
        # before_anti_aug = th_2[mid[0]-1:mid[0]-1] 

        if check_anti_aug(anti_aug):
            continue

        comp_th_2 = kmer[:stem]
        rc_th_2_up = "n" * len(comp_th_2[stem - down[1]:stem - down[0]]) # rev-comp of th_2[down[0]:down[1]]
        rc_th_2_down = comp_th_2[stem - up[1]:stem - up[0]] # rev of th_2[up[0]:up[1]]

        to_model = KOZAK + "UAU" + rc_th_2_up + "AUG" + rc_th_2_down

        result = subprocess.run(
            ["RNAfold", "--noPS"],
            input=str(to_model),
            text=True,
            capture_output=True,
            check=True
        )

        down_struct = result.stdout.split("\n")[1].split(" ")[0]

        # if ")" in down_struct or ")" in down_struct: # omit secondary structures
        #     continue


        tseq, tstruct = update_toehold(tseq, tstruct, kmer + "&", "(" * len(kmer) + "&") # transcript part
        tseq, tstruct = update_toehold(tseq, tstruct, comp_to_kmer, ")" * len(comp_to_kmer))
        tseq, tstruct = update_toehold(tseq, tstruct, downstream_aug[:-3] + "AUA", "(" * len(downstream_aug[:-3] + "AUA"))
        tseq, tstruct = update_toehold(tseq, tstruct, "AAACAA", "." * 6)
        tseq, tstruct = update_toehold(tseq, tstruct, KOZAK + "UAU" +  rc_th_2_up[3:] + "AUG" + rc_th_2_down, down_struct) # descending
        tseq, tstruct = update_toehold(tseq, tstruct, "AACCUGGCGGCAGCGCAAAAG", "......(((....))).....") # linker
        # tseq, tstruct = update_toehold(tseq, tstruct, reporter_seq, rep_sec_str) # reporter

        start = tseq.find("n")
        end = tseq.rfind("n") + 1

        num_single_n = int(tseq.count("n")/2)
        struct = "(" * num_single_n + ")" * num_single_n
        n_seq = "n" * num_single_n*2

        input_data = f"{struct}\n{n_seq}\n"

        RNAinverse = subprocess.Popen(
            ['RNAinverse', "Fp", "--noGU", "-R10"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True 
        )

        stdout, stderr = RNAinverse.communicate(input=input_data)

        out_list = stdout.strip().split("\n")

        already_checked = []

        for result in out_list:
            inv_seq = result.replace(" ", "$").split("$")[0]
            if inv_seq not in already_checked:
                already_checked.append(inv_seq)
                inv_seq_1 = inv_seq[:len(inv_seq)//2]
                inv_seq_2 = inv_seq[len(inv_seq)//2:]
                tseq = str(tseq)
                to_check_1 = re.sub(r"n+", inv_seq_1, tseq, count=1)
                complex_to_check = re.sub(r"n+", inv_seq_2, to_check_1, count=1)

                result_complex_rnaup = subprocess.run(
                    ["RNAup", "--interaction_pairwise", f"--ulength={len(th_1)}", f"--window={len(comp_to_kmer)}", 
                     "--include_both"],
                    input=str(complex_to_check),
                    text=True,
                    capture_output=True,
                    check=True,
                    cwd=tmp_dir  
                )

                rnaup_score = float(result_complex_rnaup.stdout.split("\n")[0].split("  ")[-1].split(" ")[0][1:])


                result_complex_rnacf = subprocess.run(
                    ["RNAcofold", "--noPS", '-p'], 
                    input=str(complex_to_check),
                    text=True,
                    capture_output=True,
                    check=True,
                    cwd=tmp_dir
                )

                rnacf_score = float(result_complex_rnacf.stdout.split("\n")[1].split(" ")[-1][1:-1])
                
                current_score = rnaup_score + rnacf_score

                if current_score < best_score:
                    best_score = current_score
                    best_complex = complex_to_check

print(f"Best complex: {best_complex}")
print(f"Best toehold (with linker): {best_complex.split('&')[1]}")
shutil.rmtree(tmp_dir)