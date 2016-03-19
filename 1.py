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


data = ''''''

data_dict = fasta_to_dict(data)
data = [data_dict[key] for key in data_dict]

super_str = min_super_string(data)

print(super_str)

