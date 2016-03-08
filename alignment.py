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
    return 0

def create_empty_alignment_table(n, m):
    table = [[Cell(0) for x in range(m)] for x in range(n)]
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
                table[i][j].first_string_shift = True
            if (value == incert):
                table[i][j].second_string_shift = True

    return table




t = create_alignment_table("ATA", "AGTA")
print_table(t)
