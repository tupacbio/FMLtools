#!/usr/bin/env python3

# given a specific sequence (e.g. restriction enzyme's motif), make a BED file of all perfect matches in a reference genome
# allow target sequence or reference genome to contain ambiguous bases in IUPAC format
# for efficiency, precomputes a large rainbow table of possible matches including ambiguous and lowercase bases
# best to run with pypy

import itertools, sys, genomerator

N_LIMIT = 0 # don't allow hits with more than this many Ns
REVCOMP = {
	'A': 'T',
	'C': 'G',
	'G': 'C',
	'T': 'A'
}
IUPAC = { # each key is a base code and each value is all the base codes that it can refer to
	'A': {'A'},
	'C': {'C'},
	'G': {'G'},
	'T': {'T'},
	'R': {'A', 'G'},
	'Y': {'C', 'T'},
	'S': {'G', 'C'},
	'W': {'A', 'T'},
	'K': {'G', 'T'},
	'M': {'A', 'C'},
	'B': {'C', 'G', 'T'},
	'D': {'A', 'G', 'T'},
	'H': {'A', 'C', 'T'},
	'V': {'A', 'C', 'G'},
	'N': {'A', 'C', 'G', 'T'}
}
REVERSE_IUPAC = {} # each key is one of the four standard bases and the value is all the codes that can refer to it (including lowercase)
for base in REVCOMP.keys():
	codes = {code for code, bases in IUPAC.items() if base in bases}
	REVERSE_IUPAC[base] = codes.union({code.lower() for code in codes})

def revcomp (seq):
	return ''.join(REVCOMP[base] for base in seq[::-1])


def expand_to_bases (seq):
	'''
	expand a possibly ambiguous sequence code into all possible sequences of the four standard bases
	returns a generator	
	'''
	return (''.join(seq_list) for seq_list in itertools.product(*(IUPAC[base] for base in seq)))

def expand_to_matches (seq):
	'''
	expand a sequence of the standard four bases into all possible ambiguous sequences that match it
	include possible lowercase bases
	but exclude sequences with more than the number of allowed N's
	returns a generator
	'''
	return (''.join(code_list) for code_list in itertools.product(*(REVERSE_IUPAC[base] for base in seq)) if code_list.count('N') + code_list.count('n') <= N_LIMIT)


if len(sys.argv) != 2: sys.exit('usage: %s [target sequence] < [reference FASTA] > [BED file]' % sys.argv[0])

target_seq = sys.argv[1].upper()
target_seqs_forward = set(expand_to_bases(target_seq))
target_seqs_reverse = set(revcomp(seq) for seq in target_seqs_forward)
matches_forward = set(itertools.chain(*map(expand_to_matches, target_seqs_forward)))
matches_reverse = set(itertools.chain(*map(expand_to_matches, target_seqs_reverse)))

fasta_stream = genomerator.FastaStream(
	source = sys.stdin,
	span = len(target_seq),
	overlap = True,
	include_partial = False
)

for query_feature in fasta_stream:
	query_seq = query_feature.data
	if query_seq in matches_forward: print(query_feature.bed(fasta_stream.references), file = sys.stdout)
	if query_seq in matches_reverse: print(query_feature.switched_strand().bed(fasta_stream.references), file = sys.stdout)

