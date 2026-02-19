# Input:  A set of kmers Motifs
# Output: Count(Motifs)
def Count(Motifs):
    
    count = {} # initializing the count dictionary
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count



def Profile(Motifs):
    """ Input:  A set of kmers Motifs
        Output: Profile(Motifs) """
    k = len(Motifs[0])
    t = len(Motifs)
    Profile = {}
    count = Count(Motifs)
    for symbol in "ACGT":
        Profile[symbol] = []
        for j in range(k):
            Profile[symbol].append(count[symbol][j]/ t)
    return Profile



# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus(Motifs):

    k = len(Motifs[0])
    count = Count(Motifs)
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


# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
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


# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)
def Pr(Text, Profile):   
    p = 1
    rows = "ACGT"
    for i in range(len(Text)):
        for j in range(len(rows)):
            if Text[i] == rows[j]:
                p *= Profile[rows[j]][i]
                break
    return p


def ProfileMostProbableKmer(text, k, profile):
    """ Input:  String Text, integer k, and profile matrix Profile
        Output: A profile-most probable k-mer in Text with respect to Profile. 
        If there are multiple most probable k-mers, then return the first one appearing in Text."""
    most_probable_kmer = text[0:0+k]
    m = 0
    for i in range(len(text) - k + 1):
        k_mer = text[i:i+k]
        p = Pr(k_mer, profile)
        if p > m:
            m = p
            most_probable_kmer = k_mer
    return most_probable_kmer
    

def GreedyMotifSearch(Dna, k, t):
    """ Input:  A list of strings Dna, and integers k and t
        Output: A list of strings BestMotifs resulting from applying GreedyMotifSearch(Dna, k, t). 
        If there are multiple sets of most probable motifs, then return the first one found."""
    BestMotifs = []
    n = len(Dna[0])
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])  
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
                BestMotifs = Motifs
    return BestMotifs
