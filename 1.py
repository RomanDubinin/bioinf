import re
from collections import defaultdict
import math
import itertools

rna_to_protein_dict = {"UUU": "F",    "CUU": "L", "AUU": "I", "GUU": "V",
                       "UUC": "F",    "CUC": "L", "AUC": "I", "GUC": "V",
                       "UUA": "L",    "CUA": "L", "AUA": "I", "GUA": "V",
                       "UUG": "L",    "CUG": "L", "AUG": "M", "GUG": "V",
                       "UCU": "S",    "CCU": "P", "ACU": "T", "GCU": "A",
                       "UCC": "S",    "CCC": "P", "ACC": "T", "GCC": "A",
                       "UCA": "S",    "CCA": "P", "ACA": "T", "GCA": "A",
                       "UCG": "S",    "CCG": "P", "ACG": "T", "GCG": "A",
                       "UAU": "Y",    "CAU": "H", "AAU": "N", "GAU": "D",
                       "UAC": "Y",    "CAC": "H", "AAC": "N", "GAC": "D",
                       "UAA": "Stop", "CAA": "Q", "AAA": "K", "GAA": "E",
                       "UAG": "Stop", "CAG": "Q", "AAG": "K", "GAG": "E",
                       "UGU": "C",    "CGU": "R", "AGU": "S", "GGU": "G",
                       "UGC": "C",    "CGC": "R", "AGC": "S", "GGC": "G",
                       "UGA": "Stop", "CGA": "R", "AGA": "R", "GGA": "G",
                       "UGG": "W",    "CGG": "R", "AGG": "R", "GGG": "G"}

def cg_content(str):
    dict = {"A":0,"C":0, "G":0, "T":0}

    for nuc in str:
        dict[nuc] += 1

    return (dict["C"] + dict["G"])/len(str) * 100

def fasta_to_list(data):
    res = []
    for line in data.split("\n"):
        if (re.findall(r'Rosalind', line)):
            res.append("")
        else:
            res[len(res)-1] += line
    return res

def fasta_to_dict(data):
    dict = {}
    last_line = ""
    for line in data.split("\n"):
        if (re.findall(r'Rosalind', line)):
            dict[line] = ""
            last_line = line
        else:
            dict[last_line] += line
    return dict

def create_most_likely_common_ancestor_table(strings):
    m = len(strings[0])
    table = [{"A":0,"C":0, "G":0, "T":0} for x in range(m)]

    for string in strings:
        for i in range(len(string)):
            table[i][string[i]] += 1


    return table


def most_likely_common_ancestor(strings):
    table = create_most_likely_common_ancestor_table(strings)
    ancestor = ""
    for i in range(len(table)):
        max_val = max(table[i].values())
        ancestor += [x for x,y in table[i].items() if y == max_val][0]
    return ancestor

############################################### superstring

def overlap(a, b):
    best = 0
    for i in range(1, min(len(a), len(b))+1):
        if b.startswith(a[-i:]):
            best = i
    return best

def merge(first_word, second_word, overlap):
    return first_word + second_word[overlap:]

def is_first_word(word, words):
    for j in range(len(words)):
        if word != words[j]:
            over = overlap(words[j], word)
            if over > len(word)/2 and over > len(words[j])/2:
                return False
    return True

def get_first_word(words):
    for i in range(len(words)):
        if is_first_word(words[i], words):
            return words[i]

def get_next_word(word, words):
    for j in range(len(words)):
        over = overlap(word, words[j])
        if over > len(word)/2 and over > len(words[j])/2:

            return words[j]

def min_super_string(words):
    first = get_first_word(words)
    super_string = first
    words.remove(first)
    last_added = first
    while len(words) > 0:
        next_word = get_next_word(last_added, words)
        words.remove(next_word)
        super_string = merge(super_string, next_word, overlap(last_added, next_word))
        last_added = next_word

    return super_string

def protein_mass(protein):
    weights = {
        'A':   71.03711,
        'C':   103.00919,
        'D':   115.02694,
        'E':   129.04259,
        'F':   147.06841,
        'G':   57.02146,
        'H':   137.05891,
        'I':   113.08406,
        'K':   128.09496,
        'L':   113.08406,
        'M':   131.04049,
        'N':   114.04293,
        'P':   97.05276,
        'Q':   128.05858,
        'R':   156.10111,
        'S':   87.03203,
        'T':   101.04768,
        'V':   99.06841,
        'W':   186.07931,
        'Y':   163.06333,
    }

    weight = 0
    for amino in protein:
        weight += weights[amino]
    return weight

def get_complement(dna):
    complement = { "A" : "T", "T" : "A", "C" : "G", "G" : "C"}

    result = ""
    for nuc in dna:
        result += complement[nuc]

    return result[::-1]


def get_correct_reads(all_reads):
    correct_reads = []

    for i in range(len(all_reads)):
        for j in range(i+1, len(all_reads)):
            if all_reads[i] == all_reads[j]:
                correct_reads.append(all_reads[i])

            if get_complement(all_reads[i]) == all_reads[j]:
                correct_reads.append(all_reads[i])
                correct_reads.append(all_reads[j])

    return correct_reads

def hamming_disstance(s1, s2):
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def get_error_correction(reas):
    correct = get_correct_reads(reads)
    incorrect_reads = set(reads) - set(correct)
    result = []

    for incorrect_read in incorrect_reads:
        for read in correct:
            if hamming_disstance(incorrect_read, read) == 1:
                result.append((incorrect_read, read))
                break
            if hamming_disstance(incorrect_read, get_complement(read)) == 1:
                result.append((incorrect_read, get_complement(read)))
                break

    return result

def dna_to_rna(dna):
    return dna.replace("T", "U")

def rna_to_protein(rna):
    res = ""
    for i in range(0, len(rna), 3):
        if len(rna[i:i+3]) != 3:
            return
        prot = rna_to_protein_dict[rna[i:i+3]]
        if prot == "Stop":
            return res
        else:
            res += prot

    return


def rna_splicing(dna, introns):
    new_dna = dna
    for intron in introns:
        new_dna = new_dna.replace(intron, "")

    rna = dna_to_rna(new_dna)
    proteins = rna_to_protein(rna).replace("Stop", "")
    return proteins

def number_of_different_rna_from_protien(protein):
    res = 1

    count_of_different_codons = defaultdict(int)
    for codon in rna_to_protein_dict:
        count_of_different_codons[rna_to_protein_dict[codon]] += 1

    for letter in protein:
        res = res * count_of_different_codons[letter]

    return res * 3 # 4 Stop codons

def nCr(n,r):
    f = math.factorial
    return f(n) / f(n-r)

def get_all_proteins_from_dna(dna):
    rna = dna_to_rna(dna)
    start_codon = "AUG"

    protein_starts = [a.start() for a in list(re.finditer(start_codon, rna))]

    proteins = []
    for start in protein_starts:
        prot = rna_to_protein(rna[start:])
        if prot is not None:
            proteins.append(prot)

    return proteins


def k_mers_lexicographically(alphabet, k):
    product = itertools.product(alphabet, repeat = k)
    for word in product:
        print("".join(letter for letter in word))

k = 3
alphabet = "O V Y G R K J".split(" ")

k_mers_lexicographically(alphabet, k)