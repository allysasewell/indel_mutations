import csv


def GetReverseComplement(str):
    if str == 'NA':
        return 
    if str == '-':
        return '-'
    dict = {}
    dict['A'] = 'T'
    dict['C'] = 'G'
    dict['G'] = 'C'
    dict['T'] = 'A'
    dict['-'] = '-'
    dict['*'] = '*'
    dict['N'] = 'N'
    str2 = ''
    for i in range(0, len(str)):
        str2 = str2 + dict[str[i]]
    return str2[::-1]



def GetComplement(str):
    if str == 'NA':
        return 
    if str == '-':
        return '-'
    dict = {}
    dict['A'] = 'T'
    dict['C'] = 'G'
    dict['G'] = 'C'
    dict['T'] = 'A'
    dict['-'] = '-'
    dict['*'] = '*'
    str2 = ''
    for i in range(0, len(str)):
        str2 = str2 + dict[str[i]]
    return str2




def GetSequence(file):
    f1 = open(file)
    sequence1 = f1.read()
    sequence1 = str(sequence1.strip().upper())
    sequence1 = ''.join(sequence1.split('\n'))
    return sequence1


sequences = []
sequence1 = GetSequence('chr1.txt')
sequences.append(sequence1)
sequence2 = GetSequence('chr2.txt')
sequences.append(sequence2)
sequence3 = GetSequence('chr3.txt')
sequences.append(sequence3)
sequence4 = GetSequence('chr4.txt')
sequences.append(sequence4)
sequence5 = GetSequence('chr5.txt')
sequences.append(sequence5)
sequence6 = GetSequence('chr6.txt')
sequences.append(sequence6)
sequence7 = GetSequence('chr7.txt')
sequences.append(sequence7)
sequence8 = GetSequence('chr8.txt')
sequences.append(sequence8)
sequence9 = GetSequence('chr9.txt')
sequences.append(sequence9)
sequence10 = GetSequence('chr10.txt')
sequences.append(sequence10)
sequence11 = GetSequence('chr11.txt')
sequences.append(sequence11)
sequence12 = GetSequence('chr12.txt')
sequences.append(sequence12)
sequence13 = GetSequence('chr13.txt')
sequences.append(sequence13)
sequence14 = GetSequence('chr14.txt')
sequences.append(sequence14)
sequence15 = GetSequence('chr15.txt')
sequences.append(sequence15)
sequence16 = GetSequence('chr16.txt')
sequences.append(sequence16)



#Returns the start and end positions and the length of the homopolymer containing a deleted base
def GetHomopolymerLength(sequence, mut_pos1, mut_pos, mut_position):
    #mut_position = int(len(sequence)/2) + 1
    mut_pos1 = int(mut_pos1) #0 index position
    mut_pos = int(mut_pos) #1 index position
    difference = mut_pos - mut_pos1 - 1
    mut_position2 = mut_position + difference
    end_position = -1
    lengths = {}
    for i in range(0, len(sequence)):
        flag = True
        length = 1
        n = 0
        while flag == True and (i + n + 1) < len(sequence):
            if sequence[i + n] == sequence[i + n + 1]:
                end_position = i + n + 2
                length = length + 1 
                n = n + 1
            else:
                flag = False
               
                if mut_position <= end_position  and mut_position >= end_position - (length - 1 ) or mut_position2 <= end_position  and mut_position2 >= end_position - (length - 1 ):
                     # or [(mut_position2 <= end_position)  and (mut_position2 >= end_position - (length - 1 ))]:
                     lengths[(end_position  + (mut_pos - mut_position) - (length - 1 ) , end_position  + (mut_pos - mut_position))] = length
    
        if mut_position <= end_position  and mut_position >= end_position - (length - 1 ) or mut_position2 <= end_position  and mut_position2 >= end_position - (length - 1 ):
    
            lengths[(end_position  + (mut_pos - mut_position) - (length - 1 ) , end_position  + (mut_pos - mut_position) )] = length
           
    if len(lengths) == 0:
        return 1, (mut_pos, mut_pos)
    for key in lengths:
        if lengths[key] == max(lengths.values()):
            return lengths[key], key

#Returns the start and end positions and the length of the homopolymer containing an inserted base
def GetHomopolymerLengthIns(sequence, mut_pos1, mut_pos, mut_position, strand):
    #mut_position = int(len(sequence)/2) + 1
    mut_position = mut_position 
    mut_pos1 = int(mut_pos1)  #0 index position
    mut_pos = int(mut_pos) #1 index position
    end_position = -1
    lengths = {}
    for i in range(0, len(sequence)):
        flag = True
        length = 1
        n = 0
        while flag == True and (i + n + 1) < len(sequence):
            if sequence[i + n] == sequence[i + n + 1]:
                end_position = i + n + 2
                length = length + 1 
                n = n + 1
            else:
                flag = False
                
                if mut_position <= end_position + 1   and (mut_position) >= end_position - (length - 1):
                        lengths[(end_position  + (mut_pos - mut_position) - (length - 1 ) , end_position  + (mut_pos - mut_position))] = length
              
    
       
            if mut_position <= end_position + 1   and mut_position >= end_position - (length - 1):
                lengths[(end_position  + (mut_pos - mut_position) - (length - 1 ) , end_position  + (mut_pos - mut_position))] = length
         
    
           
    if len(lengths) == 0:
        return 1, (mut_pos, mut_pos)
    for key in lengths:
        if lengths[key] == max(lengths.values()):
            return lengths[key], key
        


with open("_ROB_PL3_1-A01_COMPILED_mlprocess(Sheet1).csv", 'r', newline = '') as f:  #Input csv of photolyase data
    data = csv.reader(f)
    header = next(data)
#lines = [s.strip().split() for s in data.splitlines()]
   
    genotype = []
    chromosome = []
    position = []
    mutation_type = []
    length = []
    allele = []
    mutation = []
    frequency = []
    trinucleotide = []
    strand = []
    isolate = []
    isolate_list = {}
#lines.pop(0)
    for line in data:
    
        #if 'Mito' not in line[6] and 'ySR128' not in line[3] and line[0] != 'dRad16_25':#and line[10] == "DIP" or line[10] == "Complex DIP":
        isolate.append(line[0])
        if 'XVI' in line[1]:
            chromosome.append('chrXVI')
        elif 'XV' in line[1]:
            chromosome.append('chrXV')
        elif 'XIV' in line[1]:
            chromosome.append('chrXIV')
        elif 'XIII' in line[1]:
            chromosome.append('chrXIII')
        elif 'XII' in line[1]:
            chromosome.append('chrXII')
        elif 'XI' in line[1]:
            chromosome.append('chrXI')
        elif 'IX' in line[1]:
            chromosome.append('chrIX')
        elif 'X' in line[1]:
            chromosome.append('chrX')
        elif 'VIII' in line[1]:
            chromosome.append('chrVIII')
        elif 'VII' in line[1]:
            chromosome.append('chrVII')
        elif 'VI' in line[1]:
            chromosome.append('chrVI')
        elif 'IV' in line[1]:
            chromosome.append('chrIV')
        elif 'V' in line[1]:
            chromosome.append('chrV')
        elif 'III' in line[1]:
            chromosome.append('chrIII')
        elif 'II' in line[1]:
            chromosome.append('chrII')
        elif 'I' in line[1]:
            chromosome.append('chrI')
        else:
            chromosome.append('NA')

        position.append(line[2])
        
        mutation_type.append(line[3])
     

        length.append(line[4])

        allele.append(line[5])

        mutation.append(line[6])

        frequency.append(line[7])
        


for pos in range(0, len(position)):  #Get chromosome sequences
    #if mutation_type_2[pos] == 'Insertion' or mutation_type_2[pos] == 'Deletion':# or mutation_type_2[pos] == 'Deletion' and len(allele2[pos]) == 1:
        if chromosome[pos] == 'chrI':
            sequence = GetSequence('chr1.txt')
           
            if sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'C' or sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position[pos]) - 2: int(position[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position[pos]) - 2: int(position[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrII':
            sequence = GetSequence('chr2.txt')
            if sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'C' or sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position[pos]) - 2: int(position[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position[pos]) - 2: int(position[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrIII':
            sequence = GetSequence('chr3.txt')
            if sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'C' or sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position[pos]) - 2: int(position[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position[pos]) - 2: int(position[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrIV':
            sequence = GetSequence('chr4.txt')
            if sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'C' or sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position[pos]) - 2: int(position[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position[pos]) - 2: int(position[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrV':
            sequence = GetSequence('chr5.txt')
            if sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'C' or sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position[pos]) - 2: int(position[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position[pos]) - 2: int(position[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrVI':
            sequence = GetSequence('chr6.txt')
            if sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'C' or sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position[pos]) - 2: int(position[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position[pos]) - 2: int(position[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrVII':
            sequence = GetSequence('chr7.txt')
            if sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'C' or sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position[pos]) - 2: int(position[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position[pos]) - 2: int(position[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrVIII':
            sequence = GetSequence('chr8.txt')
            if sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'C' or sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position[pos]) - 2: int(position[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position[pos]) - 2: int(position[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrIX':
            sequence = GetSequence('chr9.txt')
            if sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'C' or sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position[pos]) - 2: int(position[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position[pos]) - 2: int(position[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrX':
            sequence = GetSequence('chr10.txt')
            if sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'C' or sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position[pos]) - 2: int(position[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position[pos]) - 2: int(position[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrXI':
            sequence = GetSequence('chr11.txt')
            if sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'C' or sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position[pos]) - 2: int(position[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position[pos]) - 2: int(position[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrXII':
            sequence = GetSequence('chr12.txt')
            if sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'C' or sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position[pos]) - 2: int(position[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position[pos]) - 2: int(position[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrXIII':
            sequence = GetSequence('chr13.txt')
            if sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'C' or sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position[pos]) - 2: int(position[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position[pos]) - 2: int(position[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrXIV':
            sequence = GetSequence('chr14.txt')
            if sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'C' or sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position[pos]) - 2: int(position[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position[pos]) - 2: int(position[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrXV':
            sequence = GetSequence('chr15.txt')
            if sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'C' or sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position[pos]) - 2: int(position[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position[pos]) - 2: int(position[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrXVI':
            sequence = GetSequence('chr16.txt')
            if sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'C' or sequence[int(position[pos]) - 2: int(position[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position[pos]) - 2: int(position[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position[pos]) - 2: int(position[pos]) + 1]))
                strand.append('-')
        else:
            trinucleotide.append('NA')
            strand.append('+')
for m in range(0, len(mutation)):
    if '-' in allele[m]:
        if mutation[m] == 'A' or mutation[m] == 'G':
            strand[m] = '-'
        else:
            strand[m] = '+'
    elif strand[m] == '-':
        mutation[m] = GetReverseComplement(mutation[m])

for i1 in range(0, len(isolate)):

    if isolate[i1] not in isolate_list.keys():
            isolate_list[isolate[i1]] = [(position[i1], str(chromosome[i1]), trinucleotide[i1], mutation[i1], strand[i1], mutation_type[i1], isolate[i1], int(length[i1]))]
    else:
        isolate_list[isolate[i1]] .append((position[i1], str(chromosome[i1]), trinucleotide[i1], mutation[i1], strand[i1], mutation_type[i1], isolate[i1], int(length[i1])))


def MutationsByIsolate(isolate_list):  #Number of indels or substitutions per isolate
    mutation_numbers = {}
    indel_numbers = {}
    substitution_numbers = {}
    for key in isolate_list.keys():
        mutation_number = 0
        indel_number = 0
        sub_number = 0
        for i in range(0, len(isolate_list[key])):
            mutation_number = mutation_number + 1
            if isolate_list[key][i][5] == 'Deletion':#and isolate_list[key][i][7] == 1:
                indel_number = indel_number + 1
            if isolate_list[key][i][5] == 'Insertion':# and isolate_list[key][i][7] == 1:
                indel_number = indel_number + 1
            if isolate_list[key][i][5] == 'SNV':
                sub_number = sub_number + 1
            
            

        mutation_numbers[key] = mutation_number
        indel_numbers[key] = indel_number
        substitution_numbers[key] = sub_number
    return mutation_numbers, indel_numbers, substitution_numbers

mutation_numbers = MutationsByIsolate(isolate_list)[0]
indel_numbers = MutationsByIsolate(isolate_list)[1]
substitution_numbers = MutationsByIsolate(isolate_list)[2]

def IndelRatio(mutation_numbers, indel_numbers):  #fraction of indels for each isolate
    ratios = {}
    for key in indel_numbers.keys():
        ratios[key] = indel_numbers[key]/ mutation_numbers[key]
    return ratios

#IndelRatio(mutation_numbers, indel_numbers)


def TotalMutationCount(mutation_numbers, indel_numbers): #counts total mutations and indels for each isolate, separated by photolyase treatment group
    CPD_photolyases = {'PL3_21':0, 'PL3_22':0, 'PL3_23':0, 'PL3_24': 0, 'PL3_25': 0, 'PL3_26':0, 'PL3_27':0, 'PL3_28':0, 'PL3_29':0
    , 'PL3_30':0, 'PL3_31':0,  'PL3_33':0, 'PL3_34':0, 'PL3_35':0, 'PL3_36':0}
    photolyase64 = {'PL3_57':0, 'PL3_58':0, 'PL3_59':0, 'PL3_60':0, 'PL3_61':0, 'PL3_62':0, 'PL3_63':0, 'PL3_64':0, 'PL3_65':0, 'PL3_66':0, 'PL3_67':0, 'PL3_68':0, 'PL3_69':0,
            'PL3_70':0, 'PL3_71':0, 'PL3_72':0, 'PL3_73':0, 'PL3_74':0, 'PL3_75':0}
    both = {'PL3_77':0, 'PL3_78':0, 'PL3_79':0, 'PL3_80':0, 'PL3_81':0, 'PL3_82':0, 'PL3_83':0, 'PL3_84':0, 'PL3_85':0, 'PL3_86':0, 'PL3_87':0, 'PL3_88':0, 'PL3_89':0,
                  'PL3_90':0, 'PL3_91':0, 'PL3_92':0, 'PL3_93':0, 'PL3_94':0, 'PL3_95':0, 'PL3_96':0}
    none = {'PL3_37':0, 'PL3_38':0, 'PL3_39':0, 'PL3_40':0, 'PL3_41':0, 'PL3_42':0, 'PL3_43':0, 'PL3_44':0, 'PL3_45':0, 'PL3_46':0, 'PL3_47':0 , 'PL3_49':0,
                  'PL3_50':0, 'PL3_51':0, 'PL3_52':0, 'PL3_54':0, 'PL3_55':0, 'PL3_56':0}
    
    CPD_photolyases2 = {'PL3_21':0, 'PL3_22':0, 'PL3_23':0, 'PL3_24': 0, 'PL3_25': 0, 'PL3_26':0, 'PL3_27':0, 'PL3_28':0, 'PL3_29':0
    , 'PL3_30':0, 'PL3_31':0, 'PL3_33':0, 'PL3_34':0, 'PL3_35':0, 'PL3_36':0}
    photolyase64_2 = {'PL3_57':0, 'PL3_58':0, 'PL3_59':0, 'PL3_60':0, 'PL3_61':0, 'PL3_62':0, 'PL3_63':0, 'PL3_64':0, 'PL3_65':0, 'PL3_66':0, 'PL3_67':0, 'PL3_68':0, 'PL3_69':0,
            'PL3_70':0, 'PL3_71':0, 'PL3_72':0, 'PL3_73':0, 'PL3_74':0, 'PL3_75':0}
    both2 = {'PL3_77':0, 'PL3_78':0, 'PL3_79':0, 'PL3_80':0, 'PL3_81':0, 'PL3_82':0, 'PL3_83':0, 'PL3_84':0, 'PL3_85':0, 'PL3_86':0, 'PL3_87':0, 'PL3_88':0, 'PL3_89':0,
                  'PL3_90':0, 'PL3_91':0, 'PL3_92':0, 'PL3_93':0, 'PL3_94':0, 'PL3_95':0, 'PL3_96':0}
    none2 = {'PL3_37':0, 'PL3_38':0, 'PL3_39':0, 'PL3_40':0, 'PL3_41':0, 'PL3_42':0, 'PL3_43':0, 'PL3_44':0, 'PL3_45':0, 'PL3_46':0, 'PL3_47':0,  'PL3_49':0,
                  'PL3_50':0, 'PL3_51':0, 'PL3_52':0,  'PL3_54':0, 'PL3_55':0, 'PL3_56':0}
    CPD_mutation_number = 0
    for key in CPD_photolyases.keys():
        CPD_photolyases[key] = mutation_numbers[key]
        CPD_mutation_number = CPD_mutation_number + mutation_numbers[key]
    CPD_average = CPD_mutation_number / len(CPD_photolyases)
    P64_mutation_number = 0
    for key in photolyase64.keys():
        photolyase64[key] = mutation_numbers[key]
        P64_mutation_number = P64_mutation_number + mutation_numbers[key]
    P64_average = P64_mutation_number / len(photolyase64)
    both_mutation_number = 0
    for key in both.keys():
        both[key] = mutation_numbers[key]
        both_mutation_number = both_mutation_number + mutation_numbers[key]
    both_average = both_mutation_number / len(both)

    none_number = 0
    for key in none.keys():
        none[key] = mutation_numbers[key]
        none_number = none_number + mutation_numbers[key]
    none_average = none_number / len(none)


    CPD_indel_number = 0
    for key in CPD_photolyases2.keys():
        CPD_photolyases2[key] = indel_numbers[key]
        CPD_indel_number = CPD_indel_number + indel_numbers[key]
    CPD_indel_average = CPD_indel_number / len(CPD_photolyases2)

    P64_indel_number = 0
    for key in photolyase64_2.keys():
        photolyase64_2[key] = indel_numbers[key]
        P64_indel_number = P64_indel_number + indel_numbers[key]
    P64_indel_average = P64_indel_number / len(photolyase64_2)

    both_indel_number = 0
    for key in both2.keys():
        both2[key] = indel_numbers[key]
        both_indel_number = both_indel_number + indel_numbers[key]
    both_indel_average = both_indel_number / len(both2)

    none_indel_number = 0
    for key in none2.keys():
        none2[key] = indel_numbers[key]
        none_indel_number = none_indel_number + indel_numbers[key]
    none_indel_average = none_indel_number / len(none2)
    print(none_indel_average)
   
   
    
    CPD_isolates = {}
    isolates_P64 = {}
    both_isolates = {}
    none_isolates = {}

    for key in CPD_photolyases.keys():
        CPD_isolates[key] = substitution_numbers[key]
    for key in photolyase64.keys():
        isolates_P64[key] = substitution_numbers[key]
    for key in both.keys():
        both_isolates[key] = substitution_numbers[key]
    for key in none.keys():
        none_isolates[key] = substitution_numbers[key]

    print(CPD_isolates)
    print(isolates_P64)
    print(both_isolates)
    print(none_isolates)




        



TotalMutationCount(mutation_numbers, indel_numbers)


CPD_photolyases = {'PL3_21':0, 'PL3_22':0, 'PL3_23':0, 'PL3_24': 0, 'PL3_25': 0, 'PL3_26':0, 'PL3_27':0, 'PL3_28':0, 'PL3_29':0
    , 'PL3_30':0, 'PL3_31':0,  'PL3_33':0, 'PL3_34':0, 'PL3_35':0, 'PL3_36':0}
photolyase64 = {'PL3_57':0, 'PL3_58':0, 'PL3_59':0, 'PL3_60':0, 'PL3_61':0, 'PL3_62':0, 'PL3_63':0, 'PL3_64':0, 'PL3_65':0, 'PL3_66':0, 'PL3_67':0, 'PL3_68':0, 'PL3_69':0,
            'PL3_70':0, 'PL3_71':0, 'PL3_72':0, 'PL3_73':0, 'PL3_74':0, 'PL3_75':0}
both = {'PL3_77':0, 'PL3_78':0, 'PL3_79':0, 'PL3_80':0, 'PL3_81':0, 'PL3_82':0, 'PL3_83':0, 'PL3_84':0, 'PL3_85':0, 'PL3_86':0, 'PL3_87':0, 'PL3_88':0, 'PL3_89':0,
                  'PL3_90':0, 'PL3_91':0, 'PL3_92':0, 'PL3_93':0, 'PL3_94':0, 'PL3_95':0, 'PL3_96':0}
none = {'PL3_37':0, 'PL3_38':0, 'PL3_39':0, 'PL3_40':0, 'PL3_41':0, 'PL3_42':0, 'PL3_43':0, 'PL3_44':0, 'PL3_45':0, 'PL3_46':0, 'PL3_47':0 , 'PL3_49':0,
                  'PL3_50':0, 'PL3_51':0, 'PL3_52':0,  'PL3_54':0, 'PL3_55':0, 'PL3_56':0}


def FindSpectrum(isolate_keys, mut_type): #Counts number of dipyrimidine mutations and finds mutation spectra for trinucleotide contexts
    spectrum = {}
    dipyrimidine_count = 0
    non_dipyrimidine_count = 0
    for i in range(0, len(isolate)):
        if isolate[i] in isolate_keys.keys():
            if mutation_type[i] == mut_type and int(length[i]) == 1:
                if (trinucleotide[i], mutation[i]) in spectrum.keys():
                    spectrum[(trinucleotide[i], mutation[i])] = spectrum[(trinucleotide[i], mutation[i])] + 1
                else:
                    spectrum[(trinucleotide[i], mutation[i]) ] = 1
                if 'CC' in trinucleotide[i] or 'CT' in trinucleotide[i] or 'TC' in trinucleotide[i] or 'TT' in trinucleotide[i]:
                    dipyrimidine_count = dipyrimidine_count + 1
                else:
                    non_dipyrimidine_count = non_dipyrimidine_count + 1
    return spectrum, dipyrimidine_count, non_dipyrimidine_count

FindSpectrum(none, 'SNV')
FindSpectrum(CPD_photolyases, 'SNV')
FindSpectrum(photolyase64, 'SNV')
FindSpectrum(both, 'SNV')

FindSpectrum(none, 'Insertion')
FindSpectrum(CPD_photolyases, 'Insertion')
FindSpectrum(photolyase64, 'Insertion')
FindSpectrum(both, 'Insertion')

FindSpectrum(none, 'Deletion')
FindSpectrum(CPD_photolyases, 'Deletion')
FindSpectrum(photolyase64, 'Deletion')
FindSpectrum(both, 'Deletion')







f1a = open("CPD_deletions.bed", 'w+')
f1b = open("6-4PP_deletions.bed", 'w+')
f1c = open("both_deletions.bed", 'w+')
f1d = open("none_deletions.bed", 'w+')
f2a = open("CPD_insertions.bed", 'w+')
f2b = open("6-4PP_insertions.bed", 'w+')
f2c = open("both_insertions.bed", 'w+')
f2d = open("none_insertions.bed", 'w+')


#Makes bed files for deletions and insertions in each treatment group
def WriteChr(chr_name, file1a, file1b, file1c, file1d,file2a, file2b, file2c, file2d, chromosome, position,
             trinucleotide, allele, mutation, mutation_type ,isolate, strand, length):
    file_data = []
    chromosome1_dict = {}
    chromosome1 =[]
    for x in range(0, len(mutation_type)):
        if mutation_type[x] == 'Insertion' or mutation_type[x] == 'Deletion':
        #if mutation_type[x] == 'TDIP':
                
                if trinucleotide[x].isalpha():
                    if chromosome[x] == chr_name: 
                        chromosome1_dict[int(position[x]), isolate[x]] = x
                


    for key1 in chromosome1_dict.keys():
        chromosome1.append((int(key1[0]), key1[1]))
    sorted_chromosome1= sorted(chromosome1)
    #sorted_chromosome1.sort()
    for y in range(0, len(sorted_chromosome1)):
        file_data.append(chromosome[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + str(int(position[chromosome1_dict[sorted_chromosome1[y]]]) - 1) + '\t' + str(position[chromosome1_dict[sorted_chromosome1[y]]]) + '\t' + trinucleotide[chromosome1_dict[sorted_chromosome1[y]]]  + '\t' + mutation[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + strand[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + allele[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + isolate[chromosome1_dict[sorted_chromosome1[y]]])
    for fd in range(0, len(file_data)):
            if  mutation_type[chromosome1_dict[sorted_chromosome1[fd]]] == 'Deletion':# and length[chromosome1_dict[sorted_chromosome1[fd]]] == 1 :
                if isolate[chromosome1_dict[sorted_chromosome1[fd]]] in CPD_photolyases.keys():
                    file1a.write(file_data[fd])
                    file1a.write('\n')
                if isolate[chromosome1_dict[sorted_chromosome1[fd]]] in photolyase64.keys():
                    file1b.write(file_data[fd]) 
                    file1b.write('\n')
                if isolate[chromosome1_dict[sorted_chromosome1[fd]]] in both.keys():
                    file1c.write(file_data[fd])
                    file1c.write('\n')
                if isolate[chromosome1_dict[sorted_chromosome1[fd]]] in none.keys():
                    file1d.write(file_data[fd]) 
                    file1d.write('\n')
            
            if mutation_type[chromosome1_dict[sorted_chromosome1[fd]]] == 'Insertion':# and length[chromosome1_dict[sorted_chromosome1[fd]]] == 1 :
                if isolate[chromosome1_dict[sorted_chromosome1[fd]]] in CPD_photolyases.keys():
                    file2a.write(file_data[fd])
                    file2a.write('\n')
                if isolate[chromosome1_dict[sorted_chromosome1[fd]]] in photolyase64.keys():
                    file2b.write(file_data[fd]) 
                    file2b.write('\n')
                if isolate[chromosome1_dict[sorted_chromosome1[fd]]] in both.keys():
                    file2c.write(file_data[fd])
                    file2c.write('\n')
                if isolate[chromosome1_dict[sorted_chromosome1[fd]]] in none.keys():
                    file2d.write(file_data[fd]) 
                    file2d.write('\n')
    
            
        
       
            
WriteChr('chrI',f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, chromosome, position,
trinucleotide, allele, mutation, mutation_type, isolate, strand, length)
WriteChr('chrII', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, chromosome, position,
trinucleotide, allele, mutation, mutation_type, isolate, strand, length)
WriteChr('chrIII',  f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, chromosome, position,
trinucleotide, allele, mutation, mutation_type, isolate, strand, length)
WriteChr('chrIV', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, chromosome, position,
trinucleotide, allele, mutation, mutation_type, isolate, strand, length)
WriteChr('chrIX', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, chromosome, position,
trinucleotide, allele, mutation, mutation_type, isolate, strand, length)
WriteChr('chrV', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, chromosome, position,
trinucleotide, allele, mutation, mutation_type, isolate, strand, length)
WriteChr('chrVI', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, chromosome, position,
trinucleotide, allele, mutation, mutation_type, isolate, strand, length)
WriteChr('chrVII', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, chromosome, position,
trinucleotide, allele, mutation, mutation_type, isolate, strand, length)
WriteChr('chrVIII', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, chromosome, position,
trinucleotide, allele, mutation, mutation_type, isolate, strand, length)
WriteChr('chrX', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, chromosome, position,
trinucleotide, allele, mutation, mutation_type, isolate, strand, length)
WriteChr('chrXI', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, chromosome, position,
trinucleotide, allele, mutation, mutation_type, isolate, strand, length)
WriteChr('chrXII', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, chromosome, position,
trinucleotide, allele, mutation, mutation_type, isolate, strand, length)
WriteChr('chrXIII', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, chromosome, position,
trinucleotide, allele, mutation, mutation_type, isolate, strand, length)
WriteChr('chrXIV', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, chromosome, position,
trinucleotide, allele, mutation, mutation_type, isolate, strand, length) 
WriteChr('chrXV', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, chromosome, position,
trinucleotide, allele, mutation, mutation_type, isolate, strand, length)
WriteChr('chrXVI', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d,chromosome, position,
trinucleotide, allele, mutation, mutation_type, isolate, strand, length)


f1a.close()
f1b.close()
f1c.close()
f1d.close()
f2a.close()
f2b.close()
f2c.close()
f2d.close()


f1a = open("CPD_deletions.bed")
f1b = open("6-4PP_deletions.bed")
f1c = open("both_deletions.bed")
f1d = open("none_deletions.bed")
f2a = open("CPD_insertions.bed")
f2b = open("6-4PP_insertions.bed")
f2c = open("both_insertions.bed")
f2d = open("none_insertions.bed")


f3a = open("CPD_deletions_multinucleotide.bed", 'w+')
f3b = open("6-4PP_deletions_mumltinucleotide.bed", 'w+')
f3c = open("both_deletions_multinucleotide.bed", 'w+')
f3d = open("none_deletions_multinucleotide.bed", 'w+')
f4a = open("CPD_insertions_multinucleotide.bed", 'w+')
f4b = open("6-4PP_insertions_multinucleotide.bed", 'w+')
f4c = open("both_insertions_multinucleotide.bed", 'w+')
f4d = open("none_insertions_multinucleotide.bed", 'w+')






def MakePentaNucFile(file_1, file_2):
    new_data = file_1.read()
    new_chromosome = []
    new_position = []
    new_position1 = []
    new_mutation = []
    new_allele = []
    new_trinucleotide = []
    new_strand = []
    new_isolate = []
    lengths = []
    new_lines = [s.strip().split() for s in new_data.splitlines()]
    for line in new_lines:
        new_chromosome.append(line[0])
        new_position1.append(int(line[1]))
        new_position.append(int(line[2]))
        new_mutation.append(line[4])
        new_allele.append(line[6])
        new_trinucleotide.append(line[3])
        new_strand.append(line[5])
        new_isolate.append(line[7])

    #Add more nucleotides with the indel in the middle
    pentanucleotide = []
    for pos in range(0, len(new_position)):
        if new_chromosome[pos] == 'chrI':
            sequence = GetSequence('chr1.txt')
            #if strand[pos] == '+':
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrII':
            sequence = GetSequence('chr2.txt')
            #if strand[pos] == '+':
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrIII':
            sequence = GetSequence('chr3.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
           
            


        elif new_chromosome[pos] == 'chrIV':
            sequence = GetSequence('chr4.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrV':
            sequence = GetSequence('chr5.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrVI':
            sequence = GetSequence('chr6.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrVII':
            sequence = GetSequence('chr7.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrVIII':
            sequence = GetSequence('chr8.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrIX':
            sequence = GetSequence('chr9.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq , int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrX':
            sequence = GetSequence('chr10.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq , int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrXI':
            sequence = GetSequence('chr11.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq , int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq , int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrXII':
            sequence = GetSequence('chr12.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq , int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq , int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrXIII':
            sequence = GetSequence('chr13.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq , int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq , int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrXIV':
            sequence = GetSequence('chr14.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq , int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq , int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrXV':
            sequence = GetSequence('chr15.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq , int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq , int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrXVI':
            sequence = GetSequence('chr16.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    pentanucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    pentanucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq , int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    pentanucleotide.append(multinucleotide_seq)
                else:
                     pentanucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                     lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq , int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            
            
    
    for i in range(0, len(new_position)):
        #if int(lengths) >= 4:
        #if int(lengths[i]) >= 4 and new_allele[i] == '-' and len(new_mutation[i]) == 1:
        if int(lengths[i]) < 4 and new_allele[i] == '-' and len(new_mutation[i]) == 1:
        #if len(new_allele[i]) == 1:
            if new_strand[i] == '-':
                file_2.write(new_chromosome[i] + '\t' + str(int(new_position[i]) - 1) + '\t'  +str(new_position[i]) + '\t' + new_trinucleotide[i] + '\t' + new_mutation[i] + '\t' + new_strand[i]  + '\t'+ pentanucleotide[i] + '\t' + str(lengths[i]) + '\t' + new_isolate[i])
                file_2.write('\n')
            if new_strand[i] == '+':
                file_2.write(new_chromosome[i] + '\t' + str(int(new_position[i]) - 1)  + '\t' + str(new_position[i]) + '\t' + new_trinucleotide[i] + '\t' +  new_mutation[i] + '\t' + new_strand[i] + '\t' + pentanucleotide[i] + '\t' + str(lengths[i]) + '\t' + new_isolate[i])
                file_2.write('\n')

    
    
    return pentanucleotide, new_strand, new_chromosome, new_position, new_isolate, new_trinucleotide, new_mutation, new_allele, lengths, new_position1

sequencefile1 = open("CPD_single_deletion_sequences.txt", 'w+')
sequencefile2 = open("CPD_single_insertion_sequences.txt", 'w+')
sequencefile3 = open("6-4PP_single_deletion_sequences.txt", 'w+')
sequencefile4 = open("6-4PP_single_insertion_sequences.txt", 'w+') 
sequencefile5 = open("both_single_deletion_sequences.txt", 'w+')
sequencefile6 = open("both_single_insertion_sequences.txt", 'w+') 
sequencefile7 = open("none_single_deletion_sequences.txt", 'w+')
sequencefile8 = open("none_single_insertion_seqnces.txt", 'w+') 



def MakeSequenceFile(file1, file2, pentanucleotide, pentanucleotide2, strand, strand2):
    for p in range(0, len(pentanucleotide)):
            if strand[p] == '-':
                file1.write(GetComplement(GetReverseComplement(pentanucleotide[p])))
                file1.write('\n')
            else:
                file1.write(pentanucleotide[p])
    for q in range(0, len(pentanucleotide2)):
            if strand2[q] == '-':
                file2.write(GetComplement(GetReverseComplement(pentanucleotide2[q])))
                file2.write('\n')
            else:
                file2.write(pentanucleotide2[q])
                file2.write('\n')

            

def GetPolymerCounts(pentanucleotide,  new_chromosome, new_position, new_isolate, new_trinucleotide, new_mutation, new_allele, lengths):
    homopolymers = {}
    homopolymer_list = {}
    homopolymer_deletions = 0
    homopolymer_insertions = 0
    polymer_counts = []
    for j in range(0, 4):
        polymer_counts.append({})
    for i in range (0, len(lengths)):
        #if len(new_mutation[i]) >= 2:
            #print(pentanucleotide[i])
            #print(new_trinucleotide[i])
            key = lengths[i]
            
            #elif strand[i] == '-':
                #key = GetHomopolymerLength(GetReverseComplement(sequence), )[1]
            homopolymers[(new_isolate[i], new_chromosome[i], new_position[i])] = key
            if int(key) >= 4:
                #homopolymers[(new_isolate[i], new_chromosome[i], new_position[i])] = key
                if '-' in new_mutation[i]:
                #if new_mutation[i] == '-':
                    homopolymer_deletions = homopolymer_deletions + 1
                if '-' in new_allele[i]:
                #if new_allele[i] == '-':
                    homopolymer_insertions = homopolymer_insertions + 1
                #homopolymer_list[key[1]] = [(new_chromosome[i], key[0][0], key[0][1])]
                print(new_trinucleotide[i])
                print(pentanucleotide[i])
                print(key)

            if new_allele[i] == 'C' or new_allele[i] == 'G':

                if key in polymer_counts[0].keys():
                    polymer_counts[0][key] = polymer_counts[0][key] + 1
                else:
                    polymer_counts[0][key] = 1
            if new_allele[i] == 'A' or new_allele[i] == 'T':

                if key in polymer_counts[1].keys():
                    polymer_counts[1][key] = polymer_counts[1][key] + 1
                else:
                    polymer_counts[1][key] = 1
            if new_mutation[i] == 'C' or new_mutation[i] == 'G':

                if key in polymer_counts[2].keys():
                    polymer_counts[2][key] = polymer_counts[2][key] + 1
                else:
                    polymer_counts[2][key] = 1
            if new_mutation[i] == 'A' or new_mutation[i] == 'T':

                if key in polymer_counts[3].keys():
                    polymer_counts[3][key] = polymer_counts[3][key] + 1
                else:
                    polymer_counts[3][key] = 1
    return polymer_counts, homopolymers, homopolymer_list


CPD_deletion_homopolymer_data = MakePentaNucFile(f1a, f3a)
CPD_deletion_polymer_counts = GetPolymerCounts(CPD_deletion_homopolymer_data[0],  CPD_deletion_homopolymer_data[2],CPD_deletion_homopolymer_data[3], CPD_deletion_homopolymer_data[4], CPD_deletion_homopolymer_data[5], CPD_deletion_homopolymer_data[6], CPD_deletion_homopolymer_data[7], CPD_deletion_homopolymer_data[8])[0]
homopolymer_data_64PP_deletion = MakePentaNucFile(f1b, f3b)
polymer_counts_64PP_deletion = GetPolymerCounts(homopolymer_data_64PP_deletion[0],  homopolymer_data_64PP_deletion[2], homopolymer_data_64PP_deletion[3], homopolymer_data_64PP_deletion[4], homopolymer_data_64PP_deletion[5], homopolymer_data_64PP_deletion[6], homopolymer_data_64PP_deletion[7], homopolymer_data_64PP_deletion[8])[0]
both_deletion_homopolymer_data = MakePentaNucFile(f1c, f3c)
both_deletion_polymer_counts = GetPolymerCounts(both_deletion_homopolymer_data[0],  both_deletion_homopolymer_data[2], both_deletion_homopolymer_data[3], both_deletion_homopolymer_data[4], both_deletion_homopolymer_data[5], both_deletion_homopolymer_data[6], both_deletion_homopolymer_data[7], both_deletion_homopolymer_data[8])[0]
none_deletion_homopolymer_data = MakePentaNucFile(f1d, f3d)
none_deletion_polymer_counts = GetPolymerCounts(none_deletion_homopolymer_data[0],  none_deletion_homopolymer_data[2], none_deletion_homopolymer_data[3], none_deletion_homopolymer_data[4], none_deletion_homopolymer_data[5], none_deletion_homopolymer_data[6], none_deletion_homopolymer_data[7], none_deletion_homopolymer_data[8])[0]


CPD_insertion_homopolymer_data = MakePentaNucFile(f2a, f4a)
CPD_insertioon_polymer_counts = GetPolymerCounts(CPD_insertion_homopolymer_data[0],  CPD_insertion_homopolymer_data[2], CPD_insertion_homopolymer_data[3], CPD_insertion_homopolymer_data[4], CPD_insertion_homopolymer_data[5], CPD_insertion_homopolymer_data[6], CPD_insertion_homopolymer_data[7], CPD_insertion_homopolymer_data[8])[0]
homopolymer_data_64PP_insertion = MakePentaNucFile(f2b, f4b)
polymer_counts_64PP_insertion = GetPolymerCounts(homopolymer_data_64PP_insertion[0],  homopolymer_data_64PP_insertion[2], homopolymer_data_64PP_insertion[3], homopolymer_data_64PP_insertion[4], homopolymer_data_64PP_insertion[5], homopolymer_data_64PP_insertion[6], homopolymer_data_64PP_insertion[7], homopolymer_data_64PP_insertion[8])[0]
both_insertion_homopolymer_data = MakePentaNucFile(f2c, f4c)
both_insertion_polymer_counts = GetPolymerCounts(both_insertion_homopolymer_data[0],  both_insertion_homopolymer_data[2], both_insertion_homopolymer_data[3], both_insertion_homopolymer_data[4], both_insertion_homopolymer_data[5], both_insertion_homopolymer_data[6], both_insertion_homopolymer_data[7], both_insertion_homopolymer_data[8])[0]
none_insertion_homopolymer_data = MakePentaNucFile(f2d, f4d)
none_insertion_polymer_counts = GetPolymerCounts(none_insertion_homopolymer_data[0],  none_insertion_homopolymer_data[2], none_insertion_homopolymer_data[3], none_insertion_homopolymer_data[4], none_insertion_homopolymer_data[5], none_insertion_homopolymer_data[6], none_insertion_homopolymer_data[7], none_insertion_homopolymer_data[8])[0]






MakeSequenceFile(sequencefile1, sequencefile2, CPD_deletion_homopolymer_data[0], CPD_insertion_homopolymer_data[0], CPD_deletion_homopolymer_data[1], CPD_insertion_homopolymer_data[1])
MakeSequenceFile(sequencefile3, sequencefile4, homopolymer_data_64PP_deletion[0], homopolymer_data_64PP_insertion[0], homopolymer_data_64PP_deletion[1], homopolymer_data_64PP_insertion[1])
MakeSequenceFile(sequencefile5, sequencefile6, both_deletion_homopolymer_data[0], both_insertion_homopolymer_data[0], both_deletion_homopolymer_data[1], both_insertion_homopolymer_data[1])
MakeSequenceFile(sequencefile5, sequencefile6,  none_deletion_homopolymer_data[0], none_insertion_homopolymer_data[0], none_deletion_homopolymer_data[1], none_insertion_homopolymer_data[1])






