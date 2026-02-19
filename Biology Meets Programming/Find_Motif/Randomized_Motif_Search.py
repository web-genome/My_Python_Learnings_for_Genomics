import  random

# Input:  Integers k, t, and N, followed by a collection of strings Dna
# Output: GibbsSampler(Dna, k, t, N)
def GibbsSampler(Dna, k, t, N):
    """ Input:  Integers k, t, and N, followed by a collection of strings Dna
        Output: GibbsSampler(Dna, k, t, N)
        Description: Run GibbsSampler(Dna, k, t, N) on the given input. Remember to call RandomMotifs(Dna, k, t) to obtain your initial motifs. 
        Note: You should use the pseudocounts version of Profile in your implementation of GibbsSampler. """
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = [] # output variable
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = Motifs[:]  # Copy to track the best motifs
    for j in range(N):
        # Randomly select one motif to remove
        i = random.randint(0, t - 1)  # Corrected index range
        removed_motif = Motifs.pop(i)
        profile = ProfileWithPseudocounts(Motifs) # Create profile from remaining motifs
        new_motif = ProfileGeneratedString(Dna[i], profile, k) # Sample a new motif from the removed sequence
        Motifs.insert(i, new_motif)# Insert new motif into motifs list
        if Score(Motifs) < Score(BestMotifs):# Update BestMotifs if current motifs are better
            BestMotifs = Motifs[:]
    return BestMotifs



def RandomMotifs(Dna, k, t):
    """ Input:  A list of strings Dna, and integers k and t
        Output: A list of strings Motifs, where Motifs[i] is a randomly selected k-mer from Dna[i]. """
    Motifs = []
    for i in range(t):
        b = random.randint(0, len(Dna[i]) - k)
        Motifs.append(Dna[i][b:b + k])
    return Motifs


def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = {}
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)

    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


def ProfileWithPseudocounts(Motifs):
    t = len(Motifs) + len(Motifs)
    k = len(Motifs[0])
    profile = {} # output variable
    count = CountWithPseudocounts(Motifs)
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            profile[symbol].append(count[symbol][j]/ t)
    return profile


def Score(Motifs):
    k = len(Motifs[0])
    t = len(Motifs)
    consensus = Consensus(Motifs)
    total = 0
    for i in range(t):
        for j in range(k):
            
            if Motifs[i][j] != consensus[j]:
                total += 1
    return total


def Consensus(Motifs):
    k = len(Motifs[0])
    count = CountWithPseudocounts(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus


def Pr(Text, Profile):

    p = 1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]
    return p


def ProfileMostProbableKmer(text, k, profile):
    
    most_probable_kmer = ""
    m = 0
    for i in range(len(text) - k + 1):
        k_mer = text[i:i+k]
        p = Pr(k_mer, profile)
        if p > m:
            m = p
            most_probable_kmer = k_mer
    return most_probable_kmer


def Normalize(Probabilities):
    NP = {}
    sum_of_all = sum(Probabilities.values())
    for i in Probabilities:
        new_p_value = Probabilities[i] / sum_of_all
        NP[i] = new_p_value
    return NP



import random
# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    kmer = '' # output variable
    # your code here
    p = random.uniform(0,1)
    cumulative_prob = 0.0
    for i,prob in Probabilities.items():
        cumulative_prob += prob
        if p < cumulative_prob:
            kmer = i
            break
    return kmer


def ProfileGeneratedString(Text, profile, k):
    # your code here
    n = len(Text)
    probabilities = {}
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)
