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

data = '''>Rosalind_1
ATCCAGCT
>Rosalind_2
GGGCAACT
>Rosalind_3
ATGGATCT
>Rosalind_4
AAGCAACC
>Rosalind_5
TTGGAACT
>Rosalind_6
ATGCCATT
>Rosalind_7
ATGGCACT'''

data_dict = fasta_to_dict(data)
table = create_most_likely_common_ancestor_table([data_dict[key] for key in data_dict])
ancestor = most_likely_common_ancestor([data_dict[key] for key in data_dict])

print(ancestor)

for nuc in ["A", "C", "G", "T"]:
    string = ""
    for i in range(len(table)):
        string += str(table[i][nuc]) + " "

    print(nuc + ":", string)


