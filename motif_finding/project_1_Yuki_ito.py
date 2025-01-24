#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 19:54:22 2024

@author: yuki
"""
# 1: Evolution as a Sequence of Mistakes
def HammingDistance(s, t):
    # Edge case
    if not isinstance(s, str) or not isinstance(t, str): # When inputs are not string
        raise ValueError("Both inputs must be strings.")
    if len(s) != len(t): # When inputs are not the same length
        raise ValueError("Strings must be of the same length.")

    count = 0
    for i in range(len(s)):
        if s[i] != t[i]:
            count += 1 # When ith nucleotide is different between two sequences, increase the count
    return count

# 2: Translating RNA into Protein
CODON_TABLE = {
    'AUG': 'M', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'UCU': 'S',
    'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'UAU': 'Y', 'UAC': 'Y', 'UGU': 'C',
    'UGC': 'C', 'UGG': 'W', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAU': 'H', 'CAC': 'H',
    'CAA': 'Q', 'CAG': 'Q', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T',
    'ACG': 'T', 'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGU': 'S',
    'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V',
    'GUG': 'V', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAU': 'D',
    'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'
}
def Translation(s):
    # Edge case: input must be a string
    if not isinstance(s, str):
        raise ValueError("Inputs must be strings.")
    # Find the start codon
    for i in range(len(s)):
        if s[i:i+3] == 'AUG':
            start = i
            break
    # else:
    #     raise ValueError("There is no start codon.")
    # Translate the protein
    prot = ""
    for j in range(start, len(s) - len(s) % 3, 3):  # Exclude leftover bases
        codon = s[j:j+3]
        amino_acid = CODON_TABLE.get(codon)
        if amino_acid == 'Stop':
            break
        elif amino_acid:
            prot += amino_acid
        else:
            # Ignore invalid codons (optional: raise an error here if needed)
            raise ValueError("Invalid codon found in the sequence.")
    return prot

# 3: Finding a Motif in DNA
def FindingMotif(s, t):
    # Edge case
    if not isinstance(s, str): # When inputs are not string
        raise ValueError("Inputs must be strings.")
    if len(s) < len(t):
        raise ValueError("t must be no longer than s")

    out = []
    for i in range(len(s)-len(t)+1): # Limit range to avoid index out of range
        for j in range(len(t)):
            if s[i+j] != t[j]: 
                break
            elif j+1 == len(t): # When t is a substring of s
                out.append(i+1)
    return out

# 4: RNA Splicing
# Function to parse fasta file
def parse_fasta(fasta_str):
    sequences = []
    current_seq = []
    for line in fasta_str.strip().splitlines():
        if line.startswith(">"):
            # Save the current sequence and reset
            if current_seq:
                sequences.append("".join(current_seq))
                current_seq = []
        else:
            current_seq.append(line.strip())
    if current_seq:
        sequences.append("".join(current_seq))  # Add the last sequence
    return sequences

def RNASplicing(s, introns):
    s = parse_fasta(s)[0]
    introns = parse_fasta(introns)
    # Edge case
    if not isinstance(s, str): # When inputs are not string
        raise ValueError("Inputs must be strings.")
    if not all(isinstance(intr, str) for intr in introns):
        raise ValueError("Introns must be provided as a list of strings.")
        
    # Delete introns from rna
    for intr in introns:
        s = s.replace(intr, "_")
    s = s.replace("_", "")
    # T is replaced by U
    mrna = s.replace("T", "U")
    # Translation
    prot = Translation(mrna)
    return prot

# 5: Finding a Shared Motif
def LongestCommonSubstring(k):
    k = parse_fasta(k)
    # Edge case
    if not all(isinstance(seq, str) for seq in k):
        raise ValueError("Input must be provided as a list of strings.")    
    k.sort(key=len)
    ref = k[0] # Shortest string as the reference
    out = "" # initiate the lengest common string
    # Loop over all substring of the reference string
    for start in range(len(ref)): 
       for end in range(len(ref), start, -1): # Loop from end to start+1
           substring = ref[start:end]
           if all(substring in seq for seq in k[1:]): # Check if the substring is in all other strings
               if len(substring) > len(out): # Update the lengest common string if substring is longer
                   out = substring                     
    return out

# 6: Finding a Spliced Motif
def FindingSubsequence(s, t):
    s = parse_fasta(s)
    t = parse_fasta(t)
    # Edge case
    if not isinstance(s, str) or not isinstance(t, str):
        raise ValueError("Input must be provided as a list of strings.")    

    out = []
    start = 0
    for nuc in t:
        for i in range(start, len(s)): # Loop from the last index to the end
            if nuc == s[i]:
                out.append(i+1) # Use 1-indexed positions
                start = i+1 # Update the start of the next search
                break
    return out

# 7: Finding a Shared Spliced Motif
def LongestCommonSubsequence(s, t):
    s = parse_fasta(s)
    t = parse_fasta(t)
    # Edge case
    if not isinstance(s, str) or not isinstance(t, str):
        raise ValueError("Input must be provided as a list of strings.")    

    # Initialize DP table
    dp = [[0] * (len(t) + 1) for _ in range(len(s) + 1)]
    # Update DP table
    # Calculate LCS length for every substring pairs
    for i in range(1, len(s) + 1):
        for j in range(1, len(t) + 1):
            if s[i - 1] == t[j - 1]: 
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])
        # Restore LCS
    lcs = []
    i, j = len(s), len(t)
    while i > 0 and j > 0:
        if s[i - 1] == t[j - 1]:
            lcs.append(s[i - 1]) # Append the string as a LCS component
            i -= 1
            j -= 1
        elif dp[i - 1][j] >= dp[i][j - 1]:
            i -= 1  # skip s to maintain LCS length
        else:
            j -= 1  # skip t to maintain LCS length
    
    # Reverse LCS list to the original order
    return ''.join(reversed(lcs))

# 8: Two Motifs, One Gene
def ShortestCommonSupersequence(s, t):
    # Edge case
    if not isinstance(s, str) or not isinstance(t, str):
        raise ValueError("Input must be provided as a list of strings.") 
        
    lcs = LongestCommonSubsequence(s, t)
    result = []
    i = j = 0  # Index
    # Insert strings that are not in LCS
    for char in lcs:
        while i < len(s) and s[i] != char:
            result.append(s[i])
            i += 1
        while j < len(t) and t[j] != char:
            result.append(t[j])
            j += 1
        result.append(char)  # Append LCS string
        i += 1
        j += 1
    # Add remained strings
    result.extend(s[i:])
    result.extend(t[j:])
    return ''.join(result)