
def HammingDistance(p, q):
    """Computes the Hamming distance between two strings"""
    D = 0 
    for i in range(len(p)):

        if p[i] != q[i]:

            D = D + 1

    return D


# Input:  Strings Pattern and Text, and an integer d
# Output: The number of times Pattern appears in Text with at most d mismatches
def ApproximatePatternCount(Text, Pattern, d):
    """Counts the number of times Pattern appears in Text with at most d mismatches"""
    count = 0 # initialize count variable
    for i in range(len(Text)-len(Pattern)+1):

        if HammingDistance(Text[i:i + len(Pattern)], Pattern) <= d:

            count += 1

    return count

### DO NOT MODIFY THE CODE BELOW THIS LINE ###
import sys
lines = sys.stdin.read().splitlines()
print((ApproximatePatternCount(lines[0],lines[1],int(lines[2]))))