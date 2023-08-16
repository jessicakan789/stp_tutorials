# 1.1 Find the square of the hypotenuse

def c_sq(a, b):
    c_sq = a**2 + b**2
    print("The square of the hypotenuse is {}".format(c_sq))
    return c_sq


c_sq(5, 6)

c_sq(10,12)  ## Answer: 244

#####################################################################################################

# 1.2 Slice string with two number pairs

def slice(string, i, j, k, l):
    """
    This function takes input of string, start index, stop index, start index, stop index
    """
    
    for num in [i, j, k, l]:
        if num > len(string):
            raise ValueError("The slicing index is bigger than the length of the string")
    
    if i > j or k > l:
            print("Warning: slicing indexes may be in the wrong order")
            
    tmp = ""
    tmp += string[i:j] + " " + string[k:l]
    return tmp


slice("TheUniversityOfManchesterFacultyofBiologyMedicineAndHealth", 3, 13, 15, 25)

slice("GCTGAGACTTCCTGGACGGGGGACAGGCTGTGGGGTTTCTCAGATAACTGGGCCCCTGCGCTCAG\
GAGGCCTTCACCCTCTGCTCTGGGTAAAGTTCATTGGAACAGAAAGAAATGGATTTATCTGCTCT\
TCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTC\
TGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTG\
AAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAG\
GAGCCTACAAGAAAGTACGAGATTTGAT", 137, 143, 274, 280)  ## Answer: 'GAAGAA GAAGAA'

#####################################################################################################

# 1.3 Return the sum of all odd integers between two numbers

def sum_odd(a, b):
    count = 0
    
    for num in range(a, b):
        if num % 2 != 0:
            count += num
            
    return count


sum_odd(50, 100)

sum_odd(10, 25)  ## Answer: 119

#####################################################################################################

# 1.4 

def blocks(string, num):
    tmp = ""
    i = 0
    j = num
    
    for a in range(int(len(string)/num)):  # More efficient way using range(0, len(dna), num)
        tmp += string[i:j] + " "
        i += num
        j += num
    
    return tmp


blocks("aggagtaagcccttgcaactggaaatacacccattg", 3)

blocks("GCTGAGACTTCCTGGACGGGGGACAGGCTGTGGGGTTTCTCAGATAACTGGGCCCCTGCGCTCAG\
GAGGCCTTCACCCTCTGCTCTGGGTAAAGTTCATTGGAACAGAAAGAAATGGATTTATCTGCTCT\
TCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTC\
TGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTG\
AAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAG\
GAGCCTACAAGAAAGTACGAGATTTGAT", 10)  # Answer: 'GCTGAGACTT CCTGGACGGG GGACAGGCTG TGGGGTTTCT CAGATAACTG GGCCCCTGCG CTCAGGAGGC CTTCACCCTC TGCTCTGGGT AAAGTTCATT GGAACAGAAA GAAATGGATT TATCTGCTCT TCGCGTTGAA GAAGTACAAA ATGTCATTAA TGCTATGCAG AAAATCTTAG AGTGTCCCAT CTGTCTGGAG TTGATCAAGG AACCTGTCTC CACAAAGTGT GACCACATAT TTTGCAAATT TTGCATGCTG AAACTTCTCA ACCAGAAGAA AGGGCCTTCA CAGTGTCCTT TATGTAAGAA TGATATAACC AAAAGGAGCC TACAAGAAAG TACGAGATTT '

#####################################################################################################

# 1.5 Transcribe DNA sequence to tRNA

def transcribe(dna):
    dna = dna.upper()
    tmp = ""
    
    for letter in dna:
        if letter == "T":
            tmp += "U"
        else:
            tmp += letter
            
    return tmp


transcribe("aggagtaagcccttgcaactggaaatacacccattg")

transcribe("GCTGAGACTTCCTGGACGGGGGACAGGCTGTGGGGTTTCTCAGATAACTGGGCCCCTGCGCTCAGGAG\
GCCTTCACCCTCTGCTCTGGGTAAAGTTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGT\
TGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGA\
TCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAAC\
CAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGAGCCTACAAGAAAG\
TACGAGATTTGAT")  ## Answer: 'GCUGAGACUUCCUGGACGGGGGACAGGCUGUGGGGUUUCUCAGAUAACUGGGCCCCUGCGCUCAGGAGGCCUUCACCCUCUGCUCUGGGUAAAGUUCAUUGGAACAGAAAGAAAUGGAUUUAUCUGCUCUUCGCGUUGAAGAAGUACAAAAUGUCAUUAAUGCUAUGCAGAAAAUCUUAGAGUGUCCCAUCUGUCUGGAGUUGAUCAAGGAACCUGUCUCCACAAAGUGUGACCACAUAUUUUGCAAAUUUUGCAUGCUGAAACUUCUCAACCAGAAGAAAGGGCCUUCACAGUGUCCUUUAUGUAAGAAUGAUAUAACCAAAAGGAGCCUACAAGAAAGUACGAGAUUUGAU'

#####################################################################################################

# 2.1 

def genbank(dna):
    dna = dna.lower()
    tmp = ""
    count = 1
    i = 0
    j = 10
    
    for a in range(int(len(dna)/60)): # More efficient way using range(0, len(dna), 60)
        tmp += str(count) + " "
        count += 60
        
        for b in range(6): # More efficient way using range(0, len(dna), 6)
            tmp += dna[i:j] + " "
            i += 3
            j += 3
            
        tmp += ' \n '
    
    result = print(tmp)
    return result


genbank("GCTGAGACTTCCTGGACGGGGGACAGGCTGTGGGGTTTCTCAGATAACTGGGCCCCTGCGCTCAGGAGGC\
CTTCACCCTCTGCTCTGGGTAAAGTTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAA\
GAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGG\
AACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAA\
AGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGAGCCTACAAGAAAGTACGAGATTT\
AGTCAACTTGTTGAAGAGCTATTGAAAATCATTTGTGCTTTTCAGCTTGACACAGGTTTGGAGTATGCAA\
ACAGCTATAATTTTGCAAAAAAGGAAAATAACTCTCCTGAACATCTAAAAGATGAAGTTTCTATCATCCA\
AAGTATGGGCTACAGAAACCGTGCCAAAAGACTTCTACAGAGTGAACCCGAAAATCCTTCCTTGCAGGAA\
ACCAGTCTCAGTGTCCAACTCTCTAACCTTGGAACTGTGAGAACTCTGAGGACAAAGCAGCGGATACAAC\
CTCAAAAGACGTCTGTCTACATTGAATTGGGATCTGATTCTTCTGAAGATACCGTTAATAAGGCAACTTA\
TTGCAGTGTGGGAGATCAAG")  

# Answer for above
"""
1 gctgagactt gagacttcct acttcctgga tcctggacgg tggacggggg acgggggaca 
 61 ggggacaggc gacaggctgt aggctgtggg ctgtggggtt tggggtttct ggtttctcag 
 121 ttctcagata tcagataact gataactggg aactgggccc tgggcccctg gcccctgcgc  
 181 cctgcgctca gcgctcagga ctcaggaggc aggaggcctt aggccttcac ccttcaccct  
 241 tcaccctctg ccctctgctc tctgctctgg gctctgggta ctgggtaaag ggtaaagttc  
 301 aaagttcatt gttcattgga cattggaaca tggaacagaa aacagaaaga agaaagaaat  
 361 aagaaatgga aaatggattt tggatttatc atttatctgc tatctgctct ctgctcttcg  
 421 ctcttcgcgt ttcgcgttga gcgttgaaga ttgaagaagt aagaagtaca aagtacaaaa  
 481 tacaaaatgt aaaatgtcat atgtcattaa tcattaatgc ttaatgctat atgctatgca  
 541 ctatgcagaa tgcagaaaat agaaaatctt aaatcttaga tcttagagtg tagagtgtcc  
 601 agtgtcccat gtcccatctg ccatctgtct tctgtctgga gtctggagtt tggagttgat  
 661 agttgatcaa tgatcaagga tcaaggaacc aggaacctgt aacctgtctc ctgtctccac 
 """
 
#####################################################################################################
 
# 2.2 Translate DNA into protein

# Table from https://gist.github.com/juanfal/09d7fb53bd367742127e17284b9c47bf
codontab = {
    'TCA': 'S',    # Serina
    'TCC': 'S',    # Serina
    'TCG': 'S',    # Serina
    'TCT': 'S',    # Serina
    'TTC': 'F',    # Fenilalanina
    'TTT': 'F',    # Fenilalanina
    'TTA': 'L',    # Leucina
    'TTG': 'L',    # Leucina
    'TAC': 'Y',    # Tirosina
    'TAT': 'Y',    # Tirosina
    'TAA': '*',    # Stop
    'TAG': '*',    # Stop
    'TGC': 'C',    # Cisteina
    'TGT': 'C',    # Cisteina
    'TGA': '*',    # Stop
    'TGG': 'W',    # Triptofano
    'CTA': 'L',    # Leucina
    'CTC': 'L',    # Leucina
    'CTG': 'L',    # Leucina
    'CTT': 'L',    # Leucina
    'CCA': 'P',    # Prolina
    'CCC': 'P',    # Prolina
    'CCG': 'P',    # Prolina
    'CCT': 'P',    # Prolina
    'CAC': 'H',    # Histidina
    'CAT': 'H',    # Histidina
    'CAA': 'Q',    # Glutamina
    'CAG': 'Q',    # Glutamina
    'CGA': 'R',    # Arginina
    'CGC': 'R',    # Arginina
    'CGG': 'R',    # Arginina
    'CGT': 'R',    # Arginina
    'ATA': 'I',    # Isoleucina
    'ATC': 'I',    # Isoleucina
    'ATT': 'I',    # Isoleucina
    'ATG': 'M',    # Methionina
    'ACA': 'T',    # Treonina
    'ACC': 'T',    # Treonina
    'ACG': 'T',    # Treonina
    'ACT': 'T',    # Treonina
    'AAC': 'N',    # Asparagina
    'AAT': 'N',    # Asparagina
    'AAA': 'K',    # Lisina
    'AAG': 'K',    # Lisina
    'AGC': 'S',    # Serina
    'AGT': 'S',    # Serina
    'AGA': 'R',    # Arginina
    'AGG': 'R',    # Arginina
    'GTA': 'V',    # Valina
    'GTC': 'V',    # Valina
    'GTG': 'V',    # Valina
    'GTT': 'V',    # Valina
    'GCA': 'A',    # Alanina
    'GCC': 'A',    # Alanina
    'GCG': 'A',    # Alanina
    'GCT': 'A',    # Alanina
    'GAC': 'D',    # Acido Aspartico
    'GAT': 'D',    # Acido Aspartico
    'GAA': 'E',    # Acido Glutamico
    'GAG': 'E',    # Acido Glutamico
    'GGA': 'G',    # Glicina
    'GGC': 'G',    # Glicina
    'GGG': 'G',    # Glicina
    'GGT': 'G'     # Glicina
}

 
def translate(dna):
    dna = dna.upper()
    protein = ""
    cursor = 0
    
    for letter in dna:
        if letter not in "ATGC":
            raise TypeError("Sequence must contain letters only of ATGC")
            
    if len(dna) % 3 != 0:
        raise TypeError("Sequence must be divisible by 3")
    
    for i in range(int(len(dna)/3)):  # More efficient way using range(0, len(dna), 3)
        codon = dna[cursor:cursor+3]
        # Get the list of dictionary values and find the index of codon
        # Get dictionary key of index (amino acid)
        amino = codontab.get(codon) 
        protein += amino
        cursor += 3
  
    return protein

 
translate("aggagtaagcccttgcaactggaaatacacccattg")
 
translate("ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCC\
CATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGA\
AACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAA")  # Answer: 'MDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLKLLNQKKGPSQCPLCKNDITK'

#####################################################################################################

# 2.3 Reverse complement DNA
## Adapted from https://stackoverflow.com/questions/68445715/how-to-transcribe-a-list-of-dna-sequences-to-rna-in-python

def complement(dna):
    dna = dna.upper()
    
    for letter in dna:
        if letter not in "ATGC":
            raise TypeError("Sequence must contain letters only of ATGC")
    
    mapping = str.maketrans("GCTA", "CGAT") # replace G with C, C with G etc.
    
    compl = dna.translate(mapping)[::-1]
    
    return compl


complement("aggagtaagcccttgcaactggaaatacacccattg")

complement("GCTGAGACTTCCTGGACGGGGGACAGGCTGTGGGGTTTCTCAGATAACTGGGCCCCTGCGCTCAGGAG\
GCCTTCACCCTCTGCTCTGGGTAAAGTTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGT\
TGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGA\
TCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAAC\
CAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGAGCCTACAAGAAAG\
TACGAGATTTGAT")

# Answer for above
'ATCAAATCTCGTACTTTCTTGTAGGCTCCTTTTGGTTATATCATTCTTACATAAAGGACACTGTGAAGGCCCTTTCTTCTGGTTGAGAAGTTTCAGCATGCAAAATTTGCAAAATATGTGGTCACACTTTGTGGAGACAGGTTCCTTGATCAACTCCAGACAGATGGGACACTCTAAGATTTTCTGCATAGCATTAATGACATTTTGTACTTCTTCAACGCGAAGAGCAGATAAATCCATTTCTTTCTGTTCCAATGAACTTTACCCAGAGCAGAGGGTGAAGGCCTCCTGAGCGCAGGGGCCCAGTTATCTGAGAAACCCCACAGCCTGTCCCCCGTCCAGGAAGTCTCAGC'

#####################################################################################################

# 2.4 Forward/Reverse translation starting on 1st/2nd/3rd bases

codontab = {
    'TCA': 'S',    # Serina
    'TCC': 'S',    # Serina
    'TCG': 'S',    # Serina
    'TCT': 'S',    # Serina
    'TTC': 'F',    # Fenilalanina
    'TTT': 'F',    # Fenilalanina
    'TTA': 'L',    # Leucina
    'TTG': 'L',    # Leucina
    'TAC': 'Y',    # Tirosina
    'TAT': 'Y',    # Tirosina
    'TAA': '*',    # Stop
    'TAG': '*',    # Stop
    'TGC': 'C',    # Cisteina
    'TGT': 'C',    # Cisteina
    'TGA': '*',    # Stop
    'TGG': 'W',    # Triptofano
    'CTA': 'L',    # Leucina
    'CTC': 'L',    # Leucina
    'CTG': 'L',    # Leucina
    'CTT': 'L',    # Leucina
    'CCA': 'P',    # Prolina
    'CCC': 'P',    # Prolina
    'CCG': 'P',    # Prolina
    'CCT': 'P',    # Prolina
    'CAC': 'H',    # Histidina
    'CAT': 'H',    # Histidina
    'CAA': 'Q',    # Glutamina
    'CAG': 'Q',    # Glutamina
    'CGA': 'R',    # Arginina
    'CGC': 'R',    # Arginina
    'CGG': 'R',    # Arginina
    'CGT': 'R',    # Arginina
    'ATA': 'I',    # Isoleucina
    'ATC': 'I',    # Isoleucina
    'ATT': 'I',    # Isoleucina
    'ATG': 'M',    # Methionina
    'ACA': 'T',    # Treonina
    'ACC': 'T',    # Treonina
    'ACG': 'T',    # Treonina
    'ACT': 'T',    # Treonina
    'AAC': 'N',    # Asparagina
    'AAT': 'N',    # Asparagina
    'AAA': 'K',    # Lisina
    'AAG': 'K',    # Lisina
    'AGC': 'S',    # Serina
    'AGT': 'S',    # Serina
    'AGA': 'R',    # Arginina
    'AGG': 'R',    # Arginina
    'GTA': 'V',    # Valina
    'GTC': 'V',    # Valina
    'GTG': 'V',    # Valina
    'GTT': 'V',    # Valina
    'GCA': 'A',    # Alanina
    'GCC': 'A',    # Alanina
    'GCG': 'A',    # Alanina
    'GCT': 'A',    # Alanina
    'GAC': 'D',    # Acido Aspartico
    'GAT': 'D',    # Acido Aspartico
    'GAA': 'E',    # Acido Glutamico
    'GAG': 'E',    # Acido Glutamico
    'GGA': 'G',    # Glicina
    'GGC': 'G',    # Glicina
    'GGG': 'G',    # Glicina
    'GGT': 'G'     # Glicina
}


def translate(dna):
    dna = dna.upper()
    protein = ""
    cursor = 0
    
    for letter in dna:
        if letter not in "ATGC":
            raise TypeError("Sequence must contain letters only of ATGC")
            
    #if len(dna) % 3 != 0:
        #raise TypeError("Sequence must be divisible by 3")
    
    for i in range(int(len(dna)/3)):
        codon = dna[cursor:cursor+3]
        # Get the list of dictionary values and find the index of codon
        # Get dictionary key of index (amino acid)
        amino = codontab.get(codon) 
        protein += amino
        cursor += 3
  
    return protein


def complement(dna):
    dna = dna.upper()
    
    for letter in dna:
        if letter not in "ATGC":
            raise TypeError("Sequence must contain letters only of ATGC")
    
    mapping = str.maketrans("GCTA", "CGAT") # replace G with C, C with G etc.
    
    compl = dna.translate(mapping)[::-1]
    
    return compl 


def combo(dna):
    dna_2 = dna[1:]
    dna_3 = dna[2:]
    dna_4 = dna[:len(dna)-1]
    dna_5 = dna[:len(dna)-2]
    
    print("Forward")
    print(translate(dna))
    print(translate(dna_2))   
    print(translate(dna_3))
            
    print("Reverse")
    print(translate(complement(dna)))
    print(translate(complement(dna_4)))
    print(translate(complement(dna_5)))
    
    return None



combo("aggagtaagcccttgcaactggaaatacacccattg")

combo("GCTGAGACTTCCTGGACGGGGGACAGGCTGTGGGGTTTCTCAGATAACTGGGCCCCTGCGCTCAGGAGGCCTTCACCC")

# Answer to above
"""
Forward
AETSWTGDRLWGFSDNWAPALRRPSP
LRLPGRGTGCGVSQITGPLRSGGLH
*DFLDGGQAVGFLR*LGPCAQEAFT
Reverse
G*RPPERRGPVI*ETPQPVPRPGSLS
GEGLLSAGAQLSEKPHSLSPVQEVS
VKAS*AQGPSYLRNPTACPPSRKSQ
"""

#####################################################################################################

# 2.5 Count single/di/tri-nucleotides

def nuc_count(dna):
    single = dict()
    di = dict()
    tri = dict()
    
    for letter in dna:
        if letter not in single.keys():
            single[letter] = 1
        else:
            single[letter] += 1
        
    for i in range(0, len(dna), 2):
        tmp = dna[i:i+2]
        if tmp not in di.keys():
            di[tmp] = 1
        else:
            di[tmp] += 1
            
    for i in range(0, len(dna), 3):
        tmp = dna[i:i+3]
        if tmp not in tri.keys():
            tri[tmp] = 1
        else:
            tri[tmp] += 1
    
    return single, di, tri


nuc_count("aggagtaagcccttgcaactggaaatacacccattg")

nuc_count("GAACCCGAAAATCCTTCCTTGCAGGAAACCAGTCTCAGTGTCCAACTCTCTAACCTTGGAACTGTGAGAA\
CTCTGAGGACAAAGCAGCGGATACAACCTCAAAAGACGTCTGTCTACATTGAATTGGGATCTGATTCTTC\
TGAAGATACCGTTAATAAGGCAACTTATTGCAGTGTGGGAGATCAAGAATTGTTACAAATCACCCCTCAA\
GGAACCAGGGATGAAATCAGTTTGGATTCTGCAAAAAAGGCTGCTTGTGAATTTTCTGAGACGGATGTAA")

# Answer to above
"""
({'G': 62, 'A': 89, 'C': 58, 'T': 71}, {'GA': 15, 'AC': 9, 'CC': 9, 'AA': 18, 'AT': 11, 'TT': 10, 'GC': 4, 'AG': 8, 'TC': 13, 'TG': 11, 'CA': 6, 'TA': 4, 'CT': 7, 'GG': 9, 'CG': 1, 'GT': 5}, {'GAA': 8, 'CCC': 1, 'AAT': 2, 'CCT': 3, 'TCC': 1, 'TTG': 4, 'CAG': 2, 'ACC': 4, 'AGT': 4, 'CTC': 2, 'GTC': 2, 'CAA': 6, 'TCT': 7, 'AAC': 1, 'CTT': 1, 'GGA': 4, 'ACT': 3, 'GTG': 2, 'AGA': 1, 'CTG': 1, 'AGG': 2, 'ACA': 1, 'AAG': 4, 'CGG': 1, 'ATA': 1, 'ACG': 2, 'TAC': 1, 'ATT': 1, 'GAT': 6, 'GTT': 1, 'GCA': 2, 'TAT': 1, 'TGC': 1, 'TTA': 1, 'ATC': 2, 'AAA': 1, 'GCT': 2, 'TGT': 1, 'TTT': 1, 'GAG': 1, 'GTA': 1, 'A': 1})
"""

#####################################################################################################

# 2.6 GC content

def gc_content(dna):
    dna = dna.upper()
    single = dict()
    
    for letter in dna:
        if letter not in single.keys():
            single[letter] = 1
        else:
            single[letter] += 1
    
    count = round((single['G'] + single['C']) / sum(single.values()) * 100, 2)
    result = "The GC content is {}%".format(count)
    
    return result


gc_content("aggagtaagcccttgcaactggaaatacacccattg")

gc_content("GAACCCGAAAATCCTTCCTTGCAGGAAACCAGTCTCAGTGTCCAACTCTCTAACCTTGGAACTGTGAGAA\
CTCTGAGGACAAAGCAGCGGATACAACCTCAAAAGACGTCTGTCTACATTGAATTGGGATCTGATTCTTC\
TGAAGATACCGTTAATAAGGCAACTTATTGCAGTGTGGGAGATCAAGAATTGTTACAAATCACCCCTCAA\
GGAACCAGGGATGAAATCAGTTTGGATTCTGCAAAAAAGGCTGCTTGTGAATTTTCTGAGACGGATGTAA")  # The GC content is 42.86%

#####################################################################################################

# 3 Analyse amino acid count in file

def analysis():
    file = input('Enter filename: ')
    amino = input('Enter 3-letter code of amino acid: ')
    amino = amino.upper()
    
    # Check input
    amino_list = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", 
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", 
    "TYR", "VAL"]
    
    if amino not in amino_list:
        raise ValueError("Amino acid does not exist or is not in 3-letter format. Please try again.")
    
    # Open file and count lines: 141
    try:
        with open(file) as f:
            num_lines = len(f.readlines())
            print("Number of lines: {}".format(num_lines))
    except FileNotFoundError:
        print("Sorry file not found")
        return None
        
    # Open file and count lines with 'ALA': 19
    with open(file) as f:
        count = 0
        lines = f.readlines()
        
        for line in lines:
            if amino in line:
                count += 1
                
        result = "Number of lines with {}: {}".format(amino, count)
        return result


analysis()

# Input
## 1abc.sec
## ALA

# Output
## Number of lines: 141
## 'Number of lines with ALA: 19'



import pandas as pd


def auto_analysis():
    file = input('Enter filename: ')
    result = dict()
    
    try:
        with open("1abc.sec") as file:
            pdFile = pd.read_csv(file, delim_whitespace=True, header=None)
            amino_acids = pdFile[3]
            
            for aa in amino_acids:
                if aa in result.keys():
                    result[aa] += 1
                else:
                    result[aa] = 1
            
            return result
            
    except FileNotFoundError:
        print("Sorry file not found")
        return None
        

# Input
## 1abc.sec

# Output
## {'VAL': 15, 'LEU': 14, 'SER': 8, 'ALA': 19, 'ASP': 7, 'LYS': 12, 'THR': 8, 'ASN': 3, 'GLY': 8, 'PHE': 8, 'ILE': 5, 'HIS': 10, 'GLU': 6, 'TYR': 4, 'ARG': 3, 'MET': 1, 'PRO': 4, 'GLN': 4, 'CYS': 2}


#####################################################################################################

# 4 Identify ORFs in FASTA file
## https://stackoverflow.com/questions/7654971/parsing-a-fasta-file-using-a-generator-python/7655072#7655072
## https://stackoverflow.com/questions/31757876/python-find-longest-orf-in-dna-sequence#:~:text=You%20should%20look%20into%20regular%20expressions%3A%20import%20re,%28re.findall%20%28r%27ATG%20%28%3F%3A%20%28%3F%21TAA%7CTAG%7CTGA%29...%29%2A%20%28%3F%3ATAA%7CTAG%7CTGA%29%27%2Cs%29%2C%20key%20%3D%20len%29

import re

codontab = {
    'TCA': 'S',    # Serina
    'TCC': 'S',    # Serina
    'TCG': 'S',    # Serina
    'TCT': 'S',    # Serina
    'TTC': 'F',    # Fenilalanina
    'TTT': 'F',    # Fenilalanina
    'TTA': 'L',    # Leucina
    'TTG': 'L',    # Leucina
    'TAC': 'Y',    # Tirosina
    'TAT': 'Y',    # Tirosina
    'TAA': '*',    # Stop
    'TAG': '*',    # Stop
    'TGC': 'C',    # Cisteina
    'TGT': 'C',    # Cisteina
    'TGA': '*',    # Stop
    'TGG': 'W',    # Triptofano
    'CTA': 'L',    # Leucina
    'CTC': 'L',    # Leucina
    'CTG': 'L',    # Leucina
    'CTT': 'L',    # Leucina
    'CCA': 'P',    # Prolina
    'CCC': 'P',    # Prolina
    'CCG': 'P',    # Prolina
    'CCT': 'P',    # Prolina
    'CAC': 'H',    # Histidina
    'CAT': 'H',    # Histidina
    'CAA': 'Q',    # Glutamina
    'CAG': 'Q',    # Glutamina
    'CGA': 'R',    # Arginina
    'CGC': 'R',    # Arginina
    'CGG': 'R',    # Arginina
    'CGT': 'R',    # Arginina
    'ATA': 'I',    # Isoleucina
    'ATC': 'I',    # Isoleucina
    'ATT': 'I',    # Isoleucina
    'ATG': 'M',    # Methionina
    'ACA': 'T',    # Treonina
    'ACC': 'T',    # Treonina
    'ACG': 'T',    # Treonina
    'ACT': 'T',    # Treonina
    'AAC': 'N',    # Asparagina
    'AAT': 'N',    # Asparagina
    'AAA': 'K',    # Lisina
    'AAG': 'K',    # Lisina
    'AGC': 'S',    # Serina
    'AGT': 'S',    # Serina
    'AGA': 'R',    # Arginina
    'AGG': 'R',    # Arginina
    'GTA': 'V',    # Valina
    'GTC': 'V',    # Valina
    'GTG': 'V',    # Valina
    'GTT': 'V',    # Valina
    'GCA': 'A',    # Alanina
    'GCC': 'A',    # Alanina
    'GCG': 'A',    # Alanina
    'GCT': 'A',    # Alanina
    'GAC': 'D',    # Acido Aspartico
    'GAT': 'D',    # Acido Aspartico
    'GAA': 'E',    # Acido Glutamico
    'GAG': 'E',    # Acido Glutamico
    'GGA': 'G',    # Glicina
    'GGC': 'G',    # Glicina
    'GGG': 'G',    # Glicina
    'GGT': 'G'     # Glicina
}


def read_fasta(file):
    name, seq = None, []  # Create an iterable where name is None and seq is a list
    
    for line in file:
        line = line.rstrip()  # Get rid of white space at the end of each line
        if line.startswith(">"):  # Find header
            if name:  # if name not empty i.e. new header found and name = old header
                yield (name, ''.join(seq))  # Generator: Add old header and associated seq to a tuple
            name, seq = line, []  # Set name = new header
        else:
            seq.append(line)  # Add sequence lines to seq
    
    if name:  # Loop for last tuple
        yield (name, ''.join(seq))  # Add last header and associated seq to a tuple
    
    return name, seq  # Gives multiple tuples (name_1, seq_1) (name_2, seq_2)


def translate(dna):
    dna = dna.upper()
    protein = ""
    cursor = 0
    
    for letter in dna:
        if letter not in "ATGC":
            raise TypeError("Sequence must contain letters only of ATGC")
    
    for i in range(int(len(dna)/3)):
        codon = dna[cursor:cursor+3]
        # Get the list of dictionary values and find the index of codon
        # Get dictionary key of index (amino acid)
        amino = codontab.get(codon) 
        protein += amino
        cursor += 3
  
    return protein


def orf(file):
    
    try:
        with open(file) as f:
            for name, seq in read_fasta(f):
                all_orfs = re.findall(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)', seq)  # Get all sequences starting with ATG, anything that isn't a stop codon in between, and finishes with stop codon
                longest = max(all_orfs, key = len)  # Get the orf with the maximum length
                aa = translate(longest)
                
                return name, aa
    
    except FileNotFoundError:
        print("Sorry file not found")
        return None



orf("U15422.fasta")


# Output
## ('>U15422.1 Homo sapiens protamine 1 (PRM1), protamine 2 (PRM2), and transition protein 2 (TNP2) genes, complete cds', 'MWCPPPTSVTSSPTNPHLSLLQPCWPCCSLGTAEGLCIPTPGFPLFPPCPHPIMFPLPGALFPDVLMTNSLTAFKICLNVTFSTRPALSTLNCHLSQLPPSGPRVLTVPYFFCVPFNIKYRLSIYHSYCFAHLLH*')




