# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: YOUR NAME HERE

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq





def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    else:
        return ''

    pass


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    rDna = ''
    for i in range(len(dna)-1, -1, -1):
        rDna = rDna + get_complement(dna[i])
    return rDna


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    returnDna = ''
    for i in range(0, len(dna), 3):
        if aa_table[dna[i:i+3]] == '|':
            break
        returnDna = returnDna + dna[i:i+3]
    return returnDna


def find_all_ORFs_oneframe(dna, i):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC",0)
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """

    seq = []
    temp = ''
    reading = False
    while True:

        if aa_table[dna[i:i+3]] == 'M':
            reading = True
        elif aa_table[dna[i:i+3]] == '|' and reading:
            reading = False
            seq.append(temp)
            temp = ''
        if reading:
            temp = temp + dna[i:i+3]

        i = i + 3
        if (len(dna) - 3) < i:
            if reading:
                temp = temp + dna[i:]
                seq.append(temp)
            break

    return seq


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    seq = []
    for i in range(3):
        seq = seq + find_all_ORFs_oneframe(dna, i)
    return seq


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    seq = []
    seq = seq + find_all_ORFs(dna)
    seq = seq + find_all_ORFs(get_reverse_complement(dna))
    return seq


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    longest = ''
    for i in find_all_ORFs_both_strands(dna):
        if len(i) > len(longest):
            longest = i
    return longest


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    longest = ''
    for i in range(num_trials):
        dna = shuffle_string(dna)
        temp = longest_ORF(dna, 0)
        if len(longest) < len(temp):
            longest = temp
    return longest


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    protein = ''
    for i in range(0, len(dna), 3):
        if len(dna) <= i + 2:
            break
        protein += aa_table[dna[i:i+3]]

    return protein


def gene_finder():
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.

        >>>from load import load_seq
        >>>dna = load_seq("/home/gsteelman/Desktop/SoftDes/GeneFinder/data/X73525.fa")
        >>>print(gene_finder(dna))

    """
    #threshold = len(longest_ORF_noncoding(dna, 1500))
    #ans = longest_ORF(dna)
    #proteins = []
    #for i in ans:
        #proteins.append(coding_strand_to_AA(i))
    dna = load_seq("./data/X73525.fa")
    return coding_strand_to_AA(longest_ORF(dna))


if __name__ == "__main__":
    gene_finder(dna)
