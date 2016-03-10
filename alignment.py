import numpy as np

gap_price = -1

class Cell(object):
    """docstring for Cell"""

    first_string_shift = False
    second_string_shift = False
    no_shift = False
    value = 0

    def __init__(self, val):
        value = val

def print_table(table):
    for row in table:
        print([cell.value for cell in row])


def similarity(nuc1, nuc2):
    if (nuc1 == nuc2):
        return 1
    return -1

def create_empty_alignment_table(n, m):
    table = [[Cell(0) for x in range(m)] for x in range(n)]
    for i in range(1, n):
        table[i][0].second_string_shift = True
    for j in range(1,m):
        table[0][j].first_string_shift = True
    return table


def create_alignment_table(seq1, seq2):
    seq1 = " " + seq1
    seq2 = " " + seq2
    table = create_empty_alignment_table(len(seq1), len(seq2))
    for j in range(1, len(seq2)):
        for i in range(1, len(seq1)):
            match = table[i-1][j-1].value + similarity(seq1[i], seq2[j])
            delete = table[i-1][j].value + gap_price
            incert = table[i][j-1].value + gap_price

            value = max(match, delete, incert)
            table[i][j].value = value

            if (value == match):
                table[i][j].no_shift = True
            if (value == delete):
                table[i][j].second_string_shift = True
            if (value == incert):
                table[i][j].first_string_shift = True

    return table

def get_alignments(sequense1, sequense2):
    table = create_alignment_table(sequense1, sequense2)
    queue = []
    queue.append((len(table)-1, len(table[0])-1, "", ""))
    res = []

    while len(queue) != 0:
        i, j, seq1, seq2 = queue.pop()
        if i == 0 and j == 0:
            res.append((seq1, seq2))
            continue
        if table[i][j].first_string_shift:
            seq1_copy = "-" + seq1
            seq2_copy = sequense2[j-1] + seq2
            queue.append((i, j-1, seq1_copy, seq2_copy))
        if table[i][j].second_string_shift:
            seq2_copy = "-" + seq2
            seq1_copy = sequense1[i-1] + seq1
            queue.append((i-1, j, seq1_copy, seq2_copy))
        if table[i][j].no_shift:
            seq1_copy = sequense1[i-1] + seq1
            seq2_copy = sequense2[j-1] + seq2
            queue.append((i-1, j-1, seq1_copy, seq2_copy))

    return res

#t = create_alignment_table("cat", "atactgcat")
#print_table(t)

alignments = get_alignments("cat", "atactctat")

for a in alignments:
    print(a[0])
    print(a[1])
    print("_______")
