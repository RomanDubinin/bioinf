import re

def cg_content(str):
    dict = {"A":0,"C":0, "G":0, "T":0}

    for nuc in str:
        dict[nuc] += 1

    return (dict["C"] + dict["G"])/len(str) * 100

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


data = '''>Rosalind_6404
CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
TCCCACTAATAATTCTGAGG
>Rosalind_5959
CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
ATATCCATTTGTCAGCAGACACGC
>Rosalind_0808
CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
TGGGAACCTGCGGGCAGTAGGTGGAAT'''

data_dict = fasta_to_dict(data)
content_dict = {}
for key in data_dict:
    content_dict[key] = cg_content(data_dict[key])


max_val = max(content_dict.values())             #finds the max value
keys = [x for x,y in content_dict.items() if y == max_val]
for key in keys:
    print(key)
    print(content_dict[key])