from bm_preproc import BoyerMoore



def naive_count(p,t):
    """Matches the occurences of pattern 'p' in text 't'."""
    occurrences = []
    alignment_counts = 0
    character_comparisions = 0

    for i in range(len(t)-len(p)+1):
        alignment_counts += 1
        match = True
        for j in range(len(p)):
            character_comparisions += 1
            if not t[i+j] == p[j]:
                match = False
                break
        if match:
            occurrences.append(i)

    print(f"Number of alignment(s) tried by naive exact match algorithm is/are {alignment_counts}.")
    print(f"Number of character(s) compared by naive exact match algorithm is/are {character_comparisions}.")
    print(f"number of occurences of pattern in the text are {occurrences}.")
    print("\n")


def boyer_moore_count(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    alignment_counts = 0
    character_comparisions = 0

    while i < len(t) - len(p) + 1:
        alignment_counts += 1
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            character_comparisions += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    
    print(f"Number of alignment(s) tried by Boyer moore algorithm is/are {alignment_counts}.")
    print(f"Number of character(s) compared by Boyer moore algorithm is/are {character_comparisions}.")
    print(f"number of occurences of pattern in the text are {occurrences}.")



def readGenome(filename):
    genome = ""
    with open(filename,'r') as f:
        for line in f:
            if not line[0] == ">":
                genome += line.rstrip()
    return genome



pattern = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
p = pattern.lower()
text = readGenome(r"F:\Coursera\Algorithms for DNA sequences\Module_2\chr1.GRCh38.excerpt.fasta")
t = text.lower()

results_1 = naive_count(p,t)
lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
p_bm = BoyerMoore(p, lowercase_alphabet)
result_2 = boyer_moore_count(p, p_bm, t)
