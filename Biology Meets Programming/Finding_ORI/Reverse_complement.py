
def Reverse(Pattern):
    """Reverses a string"""
    rev = ""
    for char in Pattern:
        rev = char + rev
    return rev


def Complement(Pattern):
    """Computes the complement of a DNA string"""
    comp = ""
    for char in Pattern:
        if char == "A":
            comp = comp + "T"
        if char == "T":
            comp = comp + "A"
        if char == "G":
            comp = comp + "C"
        if char == "C":
            comp = comp + "G"
    return comp

def ReverseComplement(Pattern):
    """Computes the reverse complement of a DNA string"""
    Pattern = Reverse(Pattern) # reverse all letters in a string
    Pattern = Complement(Pattern) # complement each letter in a string
    return Pattern
