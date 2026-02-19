def PatternCount(Text, Pattern):
    """Counts the number of times a pattern appears in a text"""
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i:i + len(Pattern)] == Pattern:
            count = count + 1
    return count

def FrequentWords(Text, k):
    """Finds the most frequent k-mers in a string"""
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for x in freq:
        if freq[x] == m:
            words.append(x)
    return words

def FrequencyMap(Text, k):
    """Creates a frequency map of all k-mers in a string"""
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        if Pattern == Pattern:
           freq[Pattern] = freq[Pattern] + 1
    return freq

def PatternMatching(Pattern_1, __Pattern_2, Genome):
    """Finds all starting positions where Pattern_1 and __Pattern_2 appear in Genome""" 
    positions_1  = []
    positions_2 = [] # output variable
    
    for i in range(len(Genome) - len(Pattern_1) + 1):
        if Genome[i:i + len(Pattern_1)] == Pattern_1:
            positions_1.append(i)
    print("positions of Pattern_1 is " + str(positions_1))
    
    
    for x in range(len(Genome) - len(__Pattern_2) + 1): 
        if Genome[x:x + len(__Pattern_2)] == __Pattern_2: 
            positions_2.append(x)
    print("positions of Pattern_2 is " + str(positions_2))    
    return positions_1, positions_2 
           

# Input:  Strings Genome and symbol
# Output: SymbolArray(Genome, symbol)
def SymbolArray(Genome, symbol):
    """Computes the symbol array for symbol in Genome"""
    
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(ExtendedGenome[i:i+(n//2)], symbol)
    return array
#This algorithym is inefficient while processing large data set like above here
# Therefore, we need to optimize the algorithm so it do not go through entire n//2 length while counting the symbols at i position in array every time it slide the window to next position.


# Input:  Strings Genome and symbol
# Output: FasterSymbolArray(Genome, symbol)
def FasterSymbolArray(Genome, symbol):
    """Computes the symbol array for symbol in Genome in a faster way"""
    array = {}
    
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(Genome[0:n//2], symbol)
    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]
        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i] - 1            
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i] + 1            
    return array

import sys
lines = sys.stdin.read().splitlines()
print(FasterSymbolArray(lines[0],lines[1]))



def SkewArray(Genome):
    """Computes the skew array of Genome"""
    skew = {} # output variable

    skew[0] = 0

    for i in range(len(Genome)):

        if Genome[i] == "A":

            skew[i+1] = skew[i]

        if Genome[i] == "T":

            skew[i+1] = skew[i]

        if Genome[i] == "G":

            skew[i+1] = skew[i] + 1

        if Genome[i] == "C":

            skew[i+1] = skew[i] - 1

    return skew

def MinimumSkew(Genome):
    """Finds all positions in Genome where SkewArray is minimum"""
    positions = [] # output variable
    skew = SkewArray(Genome)

    min_v = min(skew.values())

    

    for i in skew:

        if skew[i] == min_v:

            positions.append(i)

    return positions