# Define console output colors
class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

#Below code is for Week 1 inputs- some minor differences i.e. capitalisation
#Text is the Vibrio Cholerae ori
#k is the length of the pattern being searched

Text='ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC'

#1a. Count the occurences of a pattern in a given text
def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

print (color.RED + '''1. PatternCount(Pattern, Text):
Count the occurrences of a pattern in a given text'''+ color.END)
print(PatternCount('AGCA',Text))
print('')

#2. Map the occurrences of each k-mer in a given text
def CountDict(Text, k):
    Count = {}
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        Count[i] = PatternCount(Pattern, Text)
    return Count

print(color.RED + '''2. CountDict(Text, k):
Map the occurrences of each k-mer in a given text'''+color.END)
print (CountDict(Text,6))
print('')

# Canceled code (dont remember why)
# def FrequentWords(Text, k):
#     FrequentPatterns = []
#     Count = CountDict(Text, k)
#     m = max(Count.values())
#     for i in Count:
#         if Count[i] == m:
#             FrequentPatterns.append(Text[i:i+k])
#     return FrequentPatterns

#3a. Return only those kmer(s) with max occurrences in a given text
def FrequentWords(Text, k):
    FrequentPatterns = []
    Count = CountDict(Text, k)
    m = max(Count.values())
    for i in Count:
        if Count[i] == m:
            FrequentPatterns.append(Text[i:i+k])
    FrequentPatternsNoDuplicates = remove_duplicates(FrequentPatterns)
    return FrequentPatternsNoDuplicates

#3b. Subroutine: Remove duplicates from a given list
def remove_duplicates(Items):
    ItemsNoDuplicates=[]
    for i in Items:
        if i not in ItemsNoDuplicates:
            ItemsNoDuplicates.append(i)
    return ItemsNoDuplicates

print (color.RED + '''3a,b. FrequentWords(Text, k), remove_duplicates(Items):
Return only those kmer(s) with max occurrences in a given text'''+color.END)
print (FrequentWords(Text,6))
print('')

#3c. MY EXTENSION: Return the top hits for k-mers with occurrence in a given range AND their respective counts

def TopHits(Text, k):
    FrequentPatterns = {}
    Count = CountDict(Text, k)
    m = max(Count.values())
    for i in Count:
        if Count[i] >= m-4:
            y=Text[i:i+k]
            FrequentPatterns[y]=Count[i]
    return FrequentPatterns

print (color.BLUE + '''3c. Top_Hits(Text, k):
Return the top hits for k-mers with occurrence in a given range AND their respective counts'''+color.END)
print(TopHits(Text,6))
print('')

#3d. MY EXTENSION: Write a subroutine that orders Top_Hits from highest to lowest occurrence
print (color.BLUE + '''3d. Order(Items):
Write a subroutine that orders Top_Hits from highest to lowest occurrence'''+color.END)
print('')

#4a. Compute reverse complement for a given text
def Reverse_Complement(Pattern):
    revComp=complement(reverse(Pattern))
    return revComp

#4b. Substitute each nucleotide for their complement
def complement(Nucleotide):
    comp=''
    for i in Nucleotide:
        if i=='A':
            comp+='T'
        elif i=='T':
            comp += 'A'
        elif i=='G':
            comp += 'C'
        elif i=='C':
            comp += 'G'
        else:
            comp += 'n'
    return comp

#4c. Reverse given text string
def reverse(Pattern):
    result=''
    index=len(Pattern)-1
    while index>=0:
        result+=Pattern[index]
        index-=1
    return result

# Below code collapses reverse complement into 2 functions instead of 3.
# def reverse_complement(x):
#     revComp=''
#     for i in reverse(x):
#         if i=='A':
#             revComp+='T'
#         elif i=='T':
#             revComp += 'A'
#         elif i=='G':
#             revComp += 'C'
#         elif i=='C':
#             revComp += 'G'
#         else:
#             revComp += 'n'
#     return revComp

print (color.RED + '''4a,b,c. reverse_complement(x),complement(Nucleotide),reverse(x):
Reverse text string and return its complement'''+color.END)
print(Reverse_Complement('ATGCTxGCATCAGAGTCA'))
print('')


#4d. MY EXTENSION: Modify the PatternCount function to include the reverse complement of the given pattern

def PatternCountPlusrevComp(Pattern, Text):
    count = 0
    revComp=Reverse_Complement(Pattern)
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
        if Text[i:i+len(Pattern)] == revComp:
            count = count+1
    return count

print (color.BLUE + '''4d. PatternCountPlusrevComp(Pattern,Text):
Modify PatternCount to include the reverse complement of the given pattern'''+color.END)
print(PatternCountPlusrevComp('AGCA',Text))
print('')

#5. Modify PatternCount to return the starting index positions of each match in a string
def PatternMatching(Pattern, Genome):
    positions = []
    for i in range(len(Genome) - len(Pattern) + 1):
       if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(i)
    return positions

print (color.RED + '''5. PatternMatching(Pattern, Genome):
Modifies PatternCount to return the starting index positions of each match in a string'''+color.END)
print(PatternMatching('ATAT','GATATATGCATATACTT'))
print('')

#6a. Identify the frequency of a nucleotide in a shifting window of spanning half the genomic sequence

def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array

print (color.RED + '''6a. SymbolArray(Genome,symbol):
Identify the frequency of a nucleotide in a shifting window of spanning half the genomic sequence'''+color.END)
print(SymbolArray('AAAAGGGG','A'))
print('')

#6b. Identify the frequency of a nucleotide in a shifting window of spanning half the genomic sequence (faster algorithm)
#After scanning the starting window this function will check for the presence of symbol only in new characters as the windows shifts/iterates

def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    array[0] = PatternCount(symbol, Genome[0:n//2])
    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

print (color.RED + '''6b. FasterSymbolArray(Genome,symbol):
Identify the frequency of a nucleotide in a shifting window of spanning half the genomic sequence'''+color.END)
print(SymbolArray('AAAAGGGG','A'))
print('')

#6c. Generates a map with a running total whereby C and G are replaced with -1 or +1 as a sequence is read:

def Skew(Genome):
    skew={}
    n=len(Genome)
    skew[0]=0
    for i in range(1,(n+1)):
        skew[i]=skew[i-1]
        if Genome[i-1] == 'G':
            skew[i]=skew[i]+1
        elif Genome[i-1] == 'C':
            skew[i]=skew[i]-1
        else:
            skew[i]=skew[i]
    return skew

print (color.RED + '''6c. Skew(Genome):
Generates a map with a running total whereby C and G are replaced with -1 or +1 as a sequence is read:'''+color.END)
print (Skew('CATGGGCATCGGCCATACGCC'))
print('')

#6d. Returns the candidate ori position for a genome, where ori is defined by the skew minimum:

def MinimumSkew(Genome):
    positions=[]
    SkewMap=Skew(Genome)
    n=min(SkewMap.values())
    for i in SkewMap:
        if SkewMap[i]==n:
            positions.append(i)
    return positions

print (color.RED + '''6d. MinimumSkew(Genome):
Returns the candidate ori position for a genome, where ori is defined by the skew minimum:'''+color.END)
print (MinimumSkew('CATTCCAGTACTTCGATGATGGCGTGAAGA'))
print('')

#6e. Returns the candidate ori position for a genome, where ori is defined by the skew maximum:

def MaximumSkew(Genome):
    positions=[]
    SkewMap=Skew(Genome)
    n=max(SkewMap.values())
    for i in SkewMap:
        if SkewMap[i]==n:
            positions.append(i)
    return positions

print (color.RED + '''6e. MaximumSkew(Genome):
Returns the candidate replication terminus for a genome, where ori is defined by the skew maximum:'''+color.END)
print (MaximumSkew('GATACACTTCCCGAGTAGGTACTG'))
print('')

#7a. Returns the Hamming Distance separating two sequences of length p

def HammingDistance(p,q):
    HamDist=0
    for i in range(0,len(p)):
        if p[i] is not q[i]:
            HamDist+=1
    return HamDist

print (color.RED + '''7a. HammingDistance(p,q):
Returns the Hamming Distance separating two sequences of length p:'''+color.END)
print (HammingDistance('CTACAGCAATACGATCATATGCGGATCCGCAGTGGCCGGTAGACACACGT','CTACCCCGCTGCTCAATGACCGGGACTAAAGAGGCGAAGATTATGGTGTG'))
print('')

#7b. Returns the index positions of patterns within mismatch tolerance range 'd' (Hamming Distance)

def ApproximatePatternMatching(Pattern, Text, d):
    positions = []
    for i in range(len(Text) - len(Pattern) + 1):
        x= Pattern
        y = Text[i:i + len(Pattern)]
        if HammingDistance(x,y)<=d:
            positions.append(i)
    return positions

hampattern='GTGCCG'
hamtext='AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT'

print (color.RED + '''7b. ApproximatePatternMatching(Pattern, Text, d):
Returns the index positions of patterns within mismatch tolerance range '''+color.END)
print(ApproximatePatternMatching(hampattern, hamtext, 3))
print('')

#7c. Returns the PatternCount for a given pattern within mismatch tolerance range 'd' (Hamming Distance)

def ApproximatePatternCount(Pattern, Text, d):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        x=Pattern
        y = Text[i:i + len(Pattern)]
        if HammingDistance(x,y)<=d:
            count = count + 1
    return count

print (color.RED + '''7c. ApproximatePatternCount(Pattern, Text, d):
Returns the PatternCount for a given pattern within mismatch tolerance range 'd' (Hamming Distance) '''+color.END)
print(ApproximatePatternCount('GAGG', 'TTTAGAGCCTTCAGAGG', 2))
print('')

#8a. Returns a dictionary of count lists (count matrix) for each nucleotide in a given series of motifs.

def InitialiseMatrix(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
    return count

def Count(Motifs):
    count=InitialiseMatrix(Motifs)
    t = len(Motifs)
    k = len(Motifs[0])
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1     #for key[base] look at list index [j] in each row and add +1
    return count

Motifs=['AACGTA','CCCGTT','CACCTT','GGATTA','TTCCGG']

print (color.RED + '''8a. Count(Motifs):
Returns a dictionary of count lists (count matrix) for each nucleotide in a given series of motifs.'''+color.END)
print (Count(Motifs))
print('')

#8b. Returns a dictionary of frequency lists (frequency matrix) for each nucleotide in a given series of motifs.

def Profile(Motifs):
    profile=InitialiseMatrix(Motifs)
    t = len(Motifs)
    k = len(Motifs[0])
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            profile[symbol][j] += 1/t
    return profile

print (color.RED + '''8b. Profile(Motifs):
Returns a dictionary of frequency lists (frequency matrix) for each nucleotide in a given series of motifs.'''+color.END)
print (Profile(Motifs))
print('')

#8c. Returns a candidate motif modeled on the highest frequency bases at each index position from a series of sequences.

def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]            #I don't understand how the highest frequency symbol is being returned
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

print (color.RED + '''8c. Consensus(Motifs):
Returns a candidate motif modeled on the highest frequency bases at each index position from a series of sequences.'''+color.END)
print (Consensus(Motifs))
print('')

#8d. Returns a score for the consensus sequence. Score = sum(number of times the consensus base does NOT appear at each index)

def Score(Motifs):
    consensus=Consensus(Motifs)
    score=0
    t = len(Motifs)
    k = len(Motifs[0])
    for i in range(t):
        for j in range(k):
            if consensus[j] is not Motifs[i][j]:
                score+=1
    return score

print (color.RED + '''8d. Score(Motifs):
Returns a score for the consensus sequence. Score = sum(number of times the consensus base does NOT appear at each index)'''+color.END)
print (Score(Motifs))
print('')

#9a. Returns the probability a motif will appear given the frequency of its nucleotide profile matrix

def Pr(Text, Profile):
    probability=1
    t = len(Text)
    for i in range(t):
        probability*=Profile[Text[i]][i]
    return probability

Profile1={'A':[0.2,0.2,0.0,0.0,0.0,0.0,0.9,0.1,0.1,0.1,0.3,0.0],'C':[0.1,0.6,0.0,0.0,0.0,0.0,0.0,0.4,0.1,0.2,0.4,0.6],'G':[0.0,0.0,1.0,1.0,0.9,0.9,0.1,0.0,0.0,0.0,0.0,0.0],'T':[0.7,0.2,0.0,0.0,0.1,0.1,0.0,0.5,0.8,0.7,0.3,0.4]}
Text1='ACGGGGATTACC'
QuizProfile={'A':[.4,.3,0,.1,0,.9],'C':[.2,.3,0,.4,0,.1],'G':[.1,.3,1,.1,.5,0],'T':[.3,.1,0,.4,.5,0]}

print (color.RED + '''9a. Pr(Text, Profile):
Returns the probability a motif will appear given the frequency of its nucleotide profile matrix'''+color.END)
print(Pr(Text1,Profile1))
print(Pr('TCGGTA',QuizProfile))
print('')

#9b. Returns a dictionary with index and probability for each kmer in given text.
def ProfilePatternDict(Text, k, Profile):
    PatternDict={}
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        PatternDict[i]=Pr(Pattern,Profile)
    return PatternDict

Profile2={'A':[.2,.2,.3,.2,.3],'C':[.4,.3,.1,.5,.1],'G':[0.3,0.3,0.5,0.2,0.4],'T':[.1,.2,.1,.1,.2]}
Text2='ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCTCCGGG'

print (color.BLUE + '''9b. ProfilePatternDict(Text, k, Profile):
Returns a dictionary with index and probability for each kmer in given text.'''+color.END)
print(ProfilePatternDict(Text2,5,Profile2))
print('')

#9c. Lists the ALL the most probable motifs in text for a given k-mer and profile matrix (It is better than 9e).
def PMPPList(Text, k, Profile):
    PrPatterns = []
    PatternDict=ProfilePatternDict(Text, k, Profile)
    m = max(PatternDict.values())
    for i in PatternDict:
        if PatternDict[i] == m:
            PrPatterns.append(Text[i:i + k])
    PrPatternsNoDuplicates = remove_duplicates(PrPatterns)
    return PrPatternsNoDuplicates

print (color.BLUE + '''9c. PMPPList(Text, k, Profile):
Lists the most probable motifs in text for a given k-mer and profile matrix.'''+color.END)
print(PMPPList(Text2, 5, Profile2))
print('')

#9d. Return the x most probable motifs AND their respective probabilities

print (color.BLUE + '''9d. Top_Motifs(Text, k, Profile):
Return the x most probable motifs AND their respective probabilities.'''+color.END)
print('')

#9e. Returns only the FIRST encountered most probable motif in text for a given k-mer and profile matrix.
def ProfileMostProbablePattern(Text, k, Profile):
    motif =''
    P=-1
    for i in range(len(Text) - k + 1):
        Pattern = Text[i:i + k]
        if Pr(Pattern,Profile) >P:
            motif=Pattern
            P = Pr(Pattern, Profile)
    return motif

print (color.RED + '''9e. ProfileMostProbablePattern(Text, k, Profile):
Returns the most probable motif in text for a given k-mer and profile matrix.'''+color.END)
print(ProfileMostProbablePattern(Text2, 5, Profile2))
print('')

#9f. Returns a list of k-mers that DO NOT appear in text according to mismatch tolerance 'x'.

print (color.BLUE + '''9f. Absent_Motifs():
Returns a list of k-mers that DO NOT appear in text according to mismatch tolerance 'x'.'''+color.END)
print('')

#10a. Returns a list of the most probable motifs of length 'k' in a list of given sequences 'Dna'. 't' = len(Dna) and can be removed (it is left here because the exercise stipulates it as argument).

Dna=['GGCGTTCAGGCA','AAGAATCAGTCA','CAAGGAGTTCGC','CACGTCAATCAC','CAATAATATTCG']
DosR=['GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC','CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG','ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC','GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC','GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG','CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA','GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA','GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG','GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG','TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC']

def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])

    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):   #This code must be indented to access the Motifs[] list. The BestMotifs[] list is also accessible because it is outside all the loops.
            BestMotifs = Motifs

    return BestMotifs

print (color.RED + '''10a. GreedyMotifSearch(Dna, k, t):
Returns a list of the most probable motifs of length 'k' in a list of given sequences 'Dna'.
't' = len(Dna) and can be removed (it is left here because the exercise stipulates it as argument.'''+color.END)
print (GreedyMotifSearch(Dna, 3,5))
print (GreedyMotifSearch(DosR, 15,10))
print(Score(GreedyMotifSearch(DosR, 15,10)))
print('')

#11a. CountWithPseudocounts(Motifs):
#Input: List of strings 'Motifs'. Output: Returns a count Matrix with pseudocounts.

def InitialisePseudoMatrix(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(1)
    return count

def CountWithPseudocounts(Motifs):
    count=InitialisePseudoMatrix(Motifs)
    t = len(Motifs)
    k = len(Motifs[0])
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1     #for key[base] look at list index [j] in each row and add +1
    return count

print (color.RED + '''11a. CountWithPseudocounts(Motifs):
Input: List of strings 'Motifs'. Output: Returns a count Matrix with pseudocounts.'''+color.END)
print(InitialisePseudoMatrix(Motifs))
print('%s: %s' %('Count(Motifs)', Count(Motifs)))
print('%s: %s' %('CountWithPseudocounts(Motifs)', CountWithPseudocounts(Motifs)))
print('')

#11b. ProfileWithPseudocounts(Motifs):
#Input: List of strings 'Motifs'. Output: A dictionary of frequency lists (matrix) for each nucleotide at each position in Motifs.

def ProfileWithPseudocounts(Motifs):
    count=CountWithPseudocounts(Motifs)
    profile={}
    t=len(Motifs)
    for i in 'ACGT':
        ls=[]
        list=count[i]
        for k in range(len(list)):
            ls.append(list[k] / (t+4))
        profile[i]=ls
    return profile

print (color.RED + '''11b. ProfileWithPseudocounts(Motifs):
Input: List of strings 'Motifs'. Output: A dictionary of frequency lists (matrix) for each nucleotide at each position in Motifs.'''+color.END)
print (ProfileWithPseudocounts(Motifs))
print('')

#11c. GreedyMotifSearchWithPseudocounts(Dna, k, t)
# Input: A Dna string, a k-mer of length k, and t=len(Dna). Output: A list of probable motifs of lenght 'k' where each profile matrix was generated with pseudocounts

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])

    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):   #This code must be indented to access the Motifs[] list. The BestMotifs[] list is also accessible because it is outside all the loops.
            BestMotifs = Motifs

    return BestMotifs

print (color.RED + '''11c. GreedyMotifSearchWithPseudocounts(Dna, k, t):
Input: A Dna string, a k-mer of length k, and t=len(Dna). Output: A list of probable motifs of lenght 'k' where each profile matrix was generated with pseudocounts'''+color.END)
print('%s: %s' %('GreedyMotifSearchWithPseudocounts(Dna, 3,5)', GreedyMotifSearchWithPseudocounts(Dna, 3,5)))
print('%s: %s' %('GreedyMotifSearchWithPseudocounts(DosR, 15,10)', GreedyMotifSearchWithPseudocounts(DosR, 15,10)))
print('')

# 12a.Motifs(Profile, Dna):
# Input: A profile matrix Profile corresponding to a list of strings Dna. Output: A list of the Profile-most probable k-mers in each string from Dna.

def Motifs(Profile,Dna):
    Prlist=[]
    t=len(Dna)
    k=len(Profile['A'])
    for i in range(t):
        j=ProfileMostProbablePattern(Dna[i],k,Profile)
        Prlist.append(j)
    return Prlist

Profile3={'A':[.8,0,0,.2], 'C':[0,.6,.2,.0],'G':[.2,.2,.8,0],'T':[0,.2,0,.8]}
Dna1=['TTACCTTAAC','GATGTCTGTC','ACGGCGTTAG','CCCTAACGAG','CGTCAGAGGT']

print (color.RED + '''12a.Motifs(Profile, Dna):
# Input: A profile matrix Profile corresponding to a list of strings Dna. Output: A list of the Profile-most probable k-mers in each string from Dna.'''+color.END)
print(Motifs(Profile3,Dna1))
print('')

#12b. RandomMotifs(Dna, k):
#Input: A list of len(Dna) strings from Dna and a value k for length kmer. Output: A list of len(Dna) random kmers generated from Dna.

import random
def RandomMotifs(Dna,k,t):
    t=len(Dna)
    ranmotif=[]
    for i in range(t):
        y = random.randint(0, len(Dna[0]) - k)
        slice=Dna[i][y:y+k]
        ranmotif.append(slice)
    return ranmotif

# def RandomMotifs(Dna,k):
#     ranmotif=[]
#     for i in range(len(Dna)):
#         y = random.randint(0, len(Dna[0]) - k)
#         slice=Dna[i][y:y+k]
#         ranmotif.append(slice)
#     return ranmotif

print (color.RED + '''12b. RandomMotifs(Dna, k,t):
Input: A list of len(Dna)=t strings from Dna and a value k for length kmer. Output: A list of t random kmers generated from Dna.'''+color.END)
print (RandomMotifs(Dna,3,5))
print('')

#12c. RandomizedMotifSearch(Dna, k,t):
# Input: A list of len(Dna)=t strings from Dna and a value k for length kmer. N= the number of times the function is run.
# Output: The best scoring motifs as generated from looping (N) and scoring random slices.

def RandomizedMotifSearch(Dna,k,t,N):
    for i in range(N):
        M = RandomMotifs(Dna, k, t)
        BestMotifs=M
        while True:
            Profile = ProfileWithPseudocounts(M)
            M = Motifs(Profile, Dna)
            if Score(M) < Score(BestMotifs):
                BestMotifs = M
            else:
                return BestMotifs

print (color.RED + '''12c. RandomizedMotifSearch(Dna, k,t):
Input: A list of len(Dna)=t strings from Dna and a value k for length kmer. N= the number of times the function is run.
Output: The best scoring motifs as generated from looping (N) and scoring random slices.'''+color.END)
Dna2=['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA','GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG','TAGTACCGAGACCGAAAGAAGTATACAGGCGT','TAGATCAAGTTTCAGGTGCACGTCGGTGAACC','AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
print(RandomizedMotifSearch(Dna2,8,5,15))
print('')

#12d. Identifying the consensus sequence
print (color.BLUE + '''12d. Identifying the consensus sequence:
Use the functions generated so far to identify and score the consensus sequence for a given DNA string'''+color.END)
print('Use RandomizedMotifSearch to identify the DosR consensus')
print ('Consensus: %s   Score: %s' %((Consensus(RandomizedMotifSearch(DosR, 15,10,100))),Score(Consensus(RandomizedMotifSearch(DosR, 15,10,100)))))
print('Use GreedyMotifSearch to identify the DosR consensus')
print ('Consensus: %s   Score: %s' %((Consensus(GreedyMotifSearch(DosR, 15,10))),Score(Consensus(GreedyMotifSearch(DosR, 15,10)))))
print('')

#13a. Normalize(Probabilities):
#Input: A dictionary 'Probabilities' whose keys are k-mers and whose values are the probabilities of these k-mers (which do not necessarily sum to 1)
#Output: A dictionary in which each value in Probabilities has been divided by the sum of all values in Probabilities.

def SumDict(Probabilities):
    y=0
    for i in Probabilities.keys():
        y+=Probabilities[i]
    return y

def Normalize(Probabilities):
    norm={}
    y=SumDict(Probabilities)
    for j in Probabilities.keys():
        norm[j]=Probabilities[j]/y
    return norm

print (color.RED + '''13a. Normalize(Probabilities):
Input: A dictionary 'Probabilities' whose keys are k-mers and whose values are the probabilities of these k-mers (which do not necessarily sum to 1)
Output: A dictionary in which each value in Probabilities has been divided by the sum of all values in Probabilities.'''+color.END)
Probabilities={'A': 0.1, 'C': 0.2, 'G': 0.05, 'T': 0.1}
print('%s %s' %('SumDict is',SumDict(Probabilities)))
print('%s %s' %('Normalize is',Normalize(Probabilities)))
print('')

#13b. WeightedDie(Probabilities):
#Input: A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
#Output: A randomly chosen k-mer with respect to the values in Probabilities

def WeightedDie(Probabilities):
    Pr=0
    kmer=''
    for k in Probabilities.keys():
        p = random.uniform(0, 1)
        Pr+=Probabilities[k]
        if p<=Pr:
            kmer+=k
            break
    return kmer

print (color.RED + '''13b. WeightedDie(Probabilities):
Input: A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
Output: A randomly chosen k-mer with respect to the values in Probabilities:'''+color.END)
print('%s %s' %('Weighted Die is',WeightedDie(Normalize(Probabilities))))
print('')

#13c. WeightedDie_profilewithlist(Probabilities):
#Rebuild WeightedDie, and the necessary subroutines, to function for a dictionary Probabilities whose keys each reference a list of probabilities
#Input: A dictionary Probabilities whose keys are nucleotides and whose values contained by a list of probabilities that each represent their index in a string.
#Output: A randomly chosen k-mer with respect to the values in Probabilities whose total length matches that of the given probabilities:

def SumDict_profilenolist(Probabilities):
    y=0
    for i in Probabilities.keys():
        y+=Probabilities[i]
    return y

def SumDict_profilewithlist(Probabilities):
    ls=genlist(Probabilities)
    x=ls[0]
    sumdict={}
    for i in range(len(Probabilities[x])):
        y = 0
        for k in Probabilities.keys():
            y += Probabilities[k][i]
            sumdict[i]=y
    return sumdict

def Normalize_profilewithlist(Probabilities):
    sum_map=SumDict_profilewithlist(Probabilities)
    ls=genlist(Probabilities)
    y=ls[0]
    norm=initialisenorm(Probabilities)
    for i in range(len(Probabilities[y])):
        for j in Probabilities.keys():
            norm[j].append(Probabilities[j][i]/sum_map[i])
    return norm

def initialisenorm(Profile):
    ls=genlist(Profile)
    norm={}
    for i in range(len(ls)):
        norm[ls[i]]=[]
    return norm

def genlist(Probabilities):
    kmer=[]
    for i in Probabilities.keys():
        kmer.append(i)
    return kmer

def WeightedDie_profilewithlist(Probabilities):        # functions where the dictionary 'Probabilities' does contain lists of values.
    kmer=''
    ls=genlist(Probabilities)
    y=ls[0]
    for i in range(len(Probabilities[y])):
        pr = 0
        p = random.uniform(0, 1)
        for k in Probabilities.keys():
            pr+=Probabilities[k][i]
            if p<=pr:
                kmer+=k
                break
    return kmer

print (color.BLUE + '''13c. WeightedDie_profilewithlist(Probabilities):
Rebuild WeightedDie, and the necessary subroutines, to function for a dictionary Probabilities whose keys each reference a list of probabilities
Input: A dictionary Probabilities whose keys are nucleotides and whose values contained by a list of probabilities that each represent their index in a string.
Output: A randomly chosen k-mer with respect to the values in Probabilities whose total length matches that of the given probabilities:'''+color.END)
profile4={'A':[0.5,0.1,.1],'C':[0.3,0.2,.1],'G':[0.2,0.4,.1],'T':[0,0.3,.1]}
print('%s %s' %('WeightedDie_profilewithlist is',WeightedDie_profilewithlist(Normalize_profilewithlist(profile4))))
print('')

#13d. dictype(Probabilities):
#Create a function that identifies if an dictionary input Probabilities does, or does not, contain lists of probabilites and executes the appropriate weighted die function.
#Input: A dictionary Probabilities whose keys reference either lists or single values.
#Output: A randomly chosen k-mer with respect to the values in Probabilities whose total length matches that of the given probabilities:

def dictype(Probabilities):
    ls=genlist(Probabilities)
    y=ls[0]
    if type(Probabilities[y]) is list:
        x=WeightedDie_profilewithlist(Probabilities)
    else:
        x=WeightedDie(Probabilities)
    return x

print (color.BLUE + '''13d. dictype(Probabilities):
Create a function that identifies if an dictionary input Probabilities does, or does not, contain lists of probabilites and executes the appropriate weighted die function.
Input: A dictionary Probabilities whose keys reference either lists or single values.
Output: A randomly chosen k-mer with respect to the values in Probabilities whose total length matches that of the given probabilities:'''+color.END)
print(dictype(profile4))
print(dictype(Probabilities))
print('')

# 13e. ProfileGeneratedString(Text, profile, k):
# Assemble this code into a function ProfileGeneratedString(Text, profile, k) that takes a string Text, a profile matrix profile , and an integer k as input.
# It should then return a randomly generated k-mer from Text whose probabilities are generated from profile, as described above.
# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)

def PGS1(Text, profile, k):
    n=len(Text)
    probabilities={}
    for i in range(0,n-k+1):
        probabilities[Text[i:i + k]] = Pr(Text[i:i + k], profile)
    return probabilities

def ProfileGeneratedString(Text, profile, k):
    step1 = PGS1(Text, profile, k)
    step2 = Normalize(step1)
    return WeightedDie(step2)

print (color.RED + '''13e. ProfileGeneratedString(Text, profile, k):
Assemble this code into a function ProfileGeneratedString(Text, profile, k) that takes a string Text, a profile matrix profile , and an integer k as input.
It should then return a randomly generated k-mer from Text whose probabilities are generated from profile, as described above.
Input:  A string Text, a profile matrix Profile, and an integer k
Output: ProfileGeneratedString(Text, profile, k)'''+color.END)
teststring= 'ATTCAGGTTCAATCGGGAATGTTATGTAGC'
testprofile= {'G': [0.21428571428571427, 0.21428571428571427, 0.35714285714285715, 0.2857142857142857, 0.21428571428571427, 0.21428571428571427, 0.35714285714285715, 0.07142857142857142, 0.35714285714285715, 0.14285714285714285], 'T': [0.21428571428571427, 0.21428571428571427, 0.21428571428571427, 0.2857142857142857, 0.2857142857142857, 0.35714285714285715, 0.2857142857142857, 0.07142857142857142, 0.21428571428571427, 0.2857142857142857], 'C': [0.2857142857142857, 0.21428571428571427, 0.14285714285714285, 0.2857142857142857, 0.35714285714285715, 0.35714285714285715, 0.21428571428571427, 0.5714285714285714, 0.35714285714285715, 0.35714285714285715], 'A': [0.2857142857142857, 0.35714285714285715, 0.2857142857142857, 0.14285714285714285, 0.14285714285714285, 0.07142857142857142, 0.14285714285714285, 0.2857142857142857, 0.07142857142857142, 0.21428571428571427]}
testk= 10
print(ProfileGeneratedString(teststring,testprofile,testk))
print('')

print (color.BLUE + '''13f. randomstringtally(List):
Take a list of string inputs and return a frequency count for each unique string. No mismatch tolerance.
Input:  A list of strings
Output: A matrix with keys (unique strings) with values (frequencies)'''+color.END)

def genstringlist(x):
    stringlist=[]
    for i in range(x):
        stringlist.append(ProfileGeneratedString(teststring, testprofile, testk))
    return stringlist

pretallylist=genstringlist(40)
print(pretallylist)
print(pretallylist[3])


def randomstringtally(list):
    tally={}
    for i in range(len(list)):
        y=list[i]
        tally[i]=0
        for x in list:
            if y==x:
                tally[i]+=1
    return tally

print(randomstringtally(pretallylist))