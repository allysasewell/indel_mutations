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




def GetSequence(file):  #Return reference chromosome sequence from the SacCer3 genome
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
    mut_pos1 = int(mut_pos1) #0 base position
    mut_pos = int(mut_pos) #1 base position
    difference = mut_pos - mut_pos1 - 1
    mut_position2 = mut_position + difference
    if difference > 0:
        print(difference)
        print(sequence[mut_position - 1:mut_position2])
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


#Returns the start and end positions and the length of the homopolymer containing a inserted base
def GetHomopolymerLengthIns(sequence, mut_pos1, mut_pos, mut_position, strand):
    #mut_position = int(len(sequence)/2) + 1
    mut_position = mut_position 
    mut_pos1 = int(mut_pos1) #0 base position
    mut_pos = int(mut_pos) #1 base position 
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
        


with open("41467_2023_38255_MOESM5_ESM.csv", 'r', newline = '') as f: #Input whole-genome data as a csv file
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
    isolate = []
    position1 = []
    isolate_list = {}
    WT_genotype = []
    WT_chromosome = []
    WT_position = []
    WT_mutation_type = []
    WT_length = []
    WT_allele = []
    WT_mutation = []
    WT_frequency = []
    WT_trinucleotide = []
    WT_strand = []
    WT_isolate = []
    WT_position1 = []
    WT_isolate_list = {}
#lines.pop(0)
    for line in data:
        if ("rad30" in line[0] or 'Rad30' in line[0]) and line[3] == "15" and line[4] != "Mito":
            isolate.append(line[0])
            if 'XVI' in line[4]:
                chromosome.append('chrXVI')
            elif 'XV' in line[4]:
                chromosome.append('chrXV')
            elif 'XIV' in line[4]:
                chromosome.append('chrXIV')
            elif 'XIII' in line[4]:
                chromosome.append('chrXIII')
            elif 'XII' in line[4]:
                chromosome.append('chrXII')
            elif 'XI' in line[4]:
                chromosome.append('chrXI')
            elif 'IX' in line[4]:
                chromosome.append('chrIX')
            elif 'X' in line[4]:
                chromosome.append('chrX')
            elif 'VIII' in line[4]:
                chromosome.append('chrVIII')
            elif 'VII' in line[4]:
                chromosome.append('chrVII')
            elif 'VI' in line[4]:
                chromosome.append('chrVI')
            elif 'IV' in line[4]:
                chromosome.append('chrIV')
            elif 'V' in line[4]:
                chromosome.append('chrV')
            elif 'III' in line[4]:
                chromosome.append('chrIII')
            elif 'II' in line[4]:
                chromosome.append('chrII')
            elif 'I' in line[4]:
                chromosome.append('chrI')
            else:
                chromosome.append('NA')
            position.append(int(line[5]))

            mutation_type.append(line[7])

            length.append(int(line[6]))

            allele.append(line[8])

            mutation.append(line[9])

            position1.append(position[len(position) - 1] - length[len(length) - 1])

        
        if "WT" in line[1] and line[3] == "15" and line[4] != "Mito":
            WT_isolate.append(line[0])
            if 'XVI' in line[4]:
                WT_chromosome.append('chrXVI')
            elif 'XV' in line[4]:
                WT_chromosome.append('chrXV')
            elif 'XIV' in line[4]:
                WT_chromosome.append('chrXIV')
            elif 'XIII' in line[4]:
                WT_chromosome.append('chrXIII')
            elif 'XII' in line[4]:
                WT_chromosome.append('chrXII')
            elif 'XI' in line[4]:
                WT_chromosome.append('chrXI')
            elif 'IX' in line[4]:
                WT_chromosome.append('chrIX')
            elif 'X' in line[4]:
                WT_chromosome.append('chrX')
            elif 'VIII' in line[4]:
                WT_chromosome.append('chrVIII')
            elif 'VII' in line[4]:
                WT_chromosome.append('chrVII')
            elif 'VI' in line[4]:
                WT_chromosome.append('chrVI')
            elif 'IV' in line[4]:
                WT_chromosome.append('chrIV')
            elif 'V' in line[4]:
                WT_chromosome.append('chrV')
            elif 'III' in line[4]:
                WT_chromosome.append('chrIII')
            elif 'II' in line[4]:
                WT_chromosome.append('chrII')
            elif 'I' in line[4]:
                WT_chromosome.append('chrI')
            else:
                WT_chromosome.append('NA')
            
            WT_position.append(int(line[5]))

            WT_mutation_type.append(line[7])

            WT_length.append(int(line[6]))

            WT_allele.append(line[8])

            WT_mutation.append(line[9])

            WT_position1.append(WT_position[len(WT_position) - 1] - WT_length[len(WT_length) - 1])

def inputChromosome(position, chromosome): #Find trinucleotide sequence contexts
    trinucleotide = []
    strand = []
    for pos in range(0, len(position)):
        #if mutation_type_2[pos] == 'Insertion' or mutation_type_2[pos] == 'Deletion':# or mutation_type_2[pos] == 'Deletion' and len(allele2[pos]) == 1:
            if chromosome[pos] == 'chrI':
                sequence = GetSequence('chr1.txt')
                print(sequence[int(position[pos]) - 2: int(position[pos]) + 1])
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
    return trinucleotide, strand

input = inputChromosome(position, chromosome)
trinucleotide = input[0]
strand = input[1]
input2 = inputChromosome(WT_position, WT_chromosome)
WT_trinucleotide = input2[0]
WT_strand = input2[1]

for m in range(0, len(mutation)):
    if '-' in allele[m]:
        if mutation[m] == 'A' or mutation[m] == 'G':
            strand[m] = '-'
        else:
            strand[m] = '+'
    elif strand[m] == '-':
        mutation[m] = GetReverseComplement(mutation[m])


for m1 in range(0, len(WT_mutation)):
    if '-' in WT_allele[m1]:
        if WT_mutation[m1] == 'A' or WT_mutation[m1] == 'G':
            WT_strand[m1] = '-'
        else:
            WT_strand[m1] = '+'
    elif WT_strand[m1] == '-':
        WT_mutation[m1] = GetReverseComplement(WT_mutation[m1])

for i1 in range(0, len(isolate)):

    if isolate[i1] not in isolate_list.keys():
            isolate_list[isolate[i1]] = [(position[i1], str(chromosome[i1]), trinucleotide[i1], mutation[i1], strand[i1], mutation_type[i1], isolate[i1], int(length[i1]), allele[i1])]
    else:
        isolate_list[isolate[i1]] .append((position[i1], str(chromosome[i1]), trinucleotide[i1], mutation[i1], strand[i1], mutation_type[i1], isolate[i1], int(length[i1]), allele[i1]))

for i2 in range(0, len(WT_isolate)):

    if WT_isolate[i2] not in WT_isolate_list.keys():
            WT_isolate_list[WT_isolate[i2]] = [(WT_position[i2], str(WT_chromosome[i2]), WT_trinucleotide[i2], WT_mutation[i2], WT_strand[i2], WT_mutation_type[i2], WT_isolate[i2], int(WT_length[i2]), WT_allele[i2])]
    else:
        WT_isolate_list[WT_isolate[i2]] .append((WT_position[i2], str(WT_chromosome[i2]), WT_trinucleotide[i2], WT_mutation[i2], WT_strand[i2], WT_mutation_type[i2], WT_isolate[i2], int(WT_length[i2]), WT_allele[i2]))

indel_count = 0
sub_count = 0
double_sub_count = 0
for u in range(0, len(mutation_type)):
    if mutation_type[u] == 'DIP':
        indel_count = indel_count + 1
    if mutation_type[u] == 'SNV' or mutation_type[u] == 'SNP':
        sub_count = sub_count + 1
    if mutation_type[u] == 'DNV' or mutation_type[u] == 'DNP':
        double_sub_count = double_sub_count + 1

WT_indel_count = 0
WT_sub_count = 0
WT_double_sub_count = 0
for u in range(0, len(WT_mutation_type)):
    if WT_mutation_type[u] == 'DIP':
        WT_indel_count = WT_indel_count + 1
    if WT_mutation_type[u] == 'SNV' or WT_mutation_type[u] == 'SNP':
        WT_sub_count = WT_sub_count + 1
    if WT_mutation_type[u] == 'DNV' or WT_mutation_type[u] == 'DNP':
        WT_double_sub_count = WT_double_sub_count + 1


def MutationsPerIsolate( isolate_list): #Counts the number of insertions and deletions for each isolate and finds the average
    Deletions = []
    Insertions = []
    average_rad30_deletion = 0
    average_rad30_insertion = 0

    for key in isolate_list.keys():
        #if key[0] in isolate_keys:

        d_count = 0
        i_count = 0
        for a in range(0, len(isolate_list[key])):
            if isolate_list[key][a][5] == 'DIP'or isolate_list[key][a][5] == 'Deletion' or isolate_list[key][a][5] == 'Insertion':


                if '-' in isolate_list[key][a][3]:# and  len(isolate_list[key][a][7]) == 1:
                    d_count = d_count + 1
                elif '-' not in isolate_list[key][a][3]:# and len(isolate_list[key][a][3]) == 1:
                    i_count = i_count + 1
      
        Deletions.append(d_count)
        Insertions.append(i_count)
    
    for e in Deletions:
        average_rad30_deletion = average_rad30_deletion+ e
    average_rad30_deletion = average_rad30_deletion/len(Deletions)

    for i in Insertions:
        average_rad30_insertion = average_rad30_insertion + i
    average_rad30_insertion = average_rad30_insertion/len(Insertions)
    print(Deletions)
    print(Insertions)
    
    return Deletions, Insertions
       

MutationsPerIsolate(isolate_list)


rad30_deletions = {}  
rad30_insertions = {}    
for a in range(0, len(allele)):
    if mutation_type[a] == 'DIP':
             
                if '-' in allele[a]:
                    if length[a] in rad30_insertions.keys():
                        rad30_insertions[length[a]] = rad30_insertions[length[a]] + 1
                    else: 
                        rad30_insertions[length[a]] = 1
                elif '-' in mutation[a]:
                    if length[a] in rad30_deletions.keys():
                        rad30_deletions[length[a]] = rad30_deletions[length[a]] + 1
                    else: 
                        rad30_deletions[length[a]] = 1

WT_deletions = {}  
WT_insertions = {}    
for a in range(0, len(WT_allele)):
    if WT_mutation_type[a] == 'DIP':
             
                if '-' in WT_allele[a]:
                    if WT_length[a] in WT_insertions.keys():
                        WT_insertions[WT_length[a]] = WT_insertions[WT_length[a]] + 1
                    else: 
                        WT_insertions[WT_length[a]] = 1
                elif '-' in WT_mutation[a]:
                    if WT_length[a] in WT_deletions.keys():
                        WT_deletions[WT_length[a]] = WT_deletions[WT_length[a]] + 1
                    else: 
                        WT_deletions[WT_length[a]] = 1



f1 = open("rad30_insertions_sorted_tandem.bed", 'w+')
f2 = open("rad30_deletions_sorted_tandem.bed", 'w+')
f3 = open("rad30_all_sorted_tandem.bed", 'w+')
f4 = open("WT_insertionsv2_sorted_tandem.bed", 'w+')
f5 = open("WT_deletionsv2_sorted_tandem.bed", 'w+')
f6 = open("WT_allv2_sorted_tandem.bed", 'w+')

f7 = open("rad30_allcombinedv2_sorted_tandem.bed", 'w+')
f8 = open("WT_allcombinedv2_sorted_tandem.bed", 'w+')


#Makes bed files
def WriteChr(chr_name, file1, file2,  file3, chromosome, position, position1,
             trinucleotide, allele, mutation, genotype, mutation_type ,isolate, strand):
    file_data = []
    chromosome1_dict = {}
    chromosome1 =[]
    for x in range(0, len(mutation_type)):
        if mutation_type[x] == 'DIP'  or mutation_type[x] == 'Insertion' or mutation_type[x] == 'Deletion':
        #if mutation_type[x] == 'DIP' or mutation_type[x] == 'Complex DIP'  or mutation_type[x] == 'Insertion' or mutation_type[x] == 'Deletion':
        #if mutation_type[x] == 'SNV' or mutation_type[x] == 'SNP':
                
                if trinucleotide[x].isalpha():# len(trinucleotide[x]) == 3 and trinucleotide[x].isalpha():
                    if chromosome[x] == chr_name: 
                        chromosome1_dict[int(position[x]), isolate[x]] = x
                
    for key1 in chromosome1_dict.keys():
        chromosome1.append((int(key1[0]), key1[1]))
    sorted_chromosome1= sorted(chromosome1)
    #sorted_chromosome1.sort()
    for y in range(0, len(sorted_chromosome1)):
        file_data.append(chromosome[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + str(int(position1[chromosome1_dict[sorted_chromosome1[y]]] )) + '\t' + str(position[chromosome1_dict[sorted_chromosome1[y]]]) + '\t' + trinucleotide[chromosome1_dict[sorted_chromosome1[y]]]  + '\t' + mutation[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + strand[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + allele[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + isolate[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + mutation_type[chromosome1_dict[sorted_chromosome1[y]]])
        #print(chromosome[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + str(int(position1[chromosome1_dict[sorted_chromosome1[y]]]) - 1) + '\t' + position2[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + trinucleotide[chromosome1_dict[sorted_chromosome1[y]]]  + '\t' + mutation[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + mutation_type[chromosome1_dict[sorted_chromosome1[y]]])
    for fd in range(0, len(file_data)):
        if '-' in allele[chromosome1_dict[sorted_chromosome1[fd]]] and '-' not in mutation[chromosome1_dict[sorted_chromosome1[fd]]]:
            file1.write(file_data[fd])
            file3.write(file_data[fd])
            file1.write('\n')
            file3.write('\n')
            

        elif '-' in mutation[chromosome1_dict[sorted_chromosome1[fd]]] and '-' not in allele[chromosome1_dict[sorted_chromosome1[fd]]]:# and len(trinucleotide[chromosome1_dict[sorted_chromosome1[fd]]]) == 3:
            file2.write(file_data[fd])
            file3.write(file_data[fd])
            file2.write('\n')
            file3.write('\n')
                
        else:
            file3.write(file_data[fd])
            file3.write('\n')

#Makes bed files with all the mutations    
def WriteAllChr(chr_name, file1a, chromosome, position2, position1,
             trinucleotide, allele, mutation, genotype, mutation_type ,isolate, strand):
    file_data = []
    chromosome1_dict = {}
    chromosome1 =[]
    for x in range(0, len(mutation_type)):
        #if mutation_type[x] == 'SNV' or mutation_type[x] == 'SNP':        
                #if trinucleotide[x].isalpha():# len(trinucleotide[x]) == 3 and trinucleotide[x].isalpha():
                    if chromosome[x] == chr_name: 
                        chromosome1_dict[int(position2[x]), isolate[x]] = x
                


    for key1 in chromosome1_dict.keys():
        chromosome1.append((int(key1[0]), key1[1]))
    sorted_chromosome1= sorted(chromosome1)
    #sorted_chromosome1.sort()
    for y in range(0, len(sorted_chromosome1)):
        file_data.append(chromosome[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + str(int(position1[chromosome1_dict[sorted_chromosome1[y]]])) + '\t' + str(position2[chromosome1_dict[sorted_chromosome1[y]]]) + '\t' + trinucleotide[chromosome1_dict[sorted_chromosome1[y]]]  + '\t' + mutation[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + strand[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + allele[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + isolate[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + mutation_type[chromosome1_dict[sorted_chromosome1[y]]])
        #print(chromosome[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + str(int(position1[chromosome1_dict[sorted_chromosome1[y]]]) - 1) + '\t' + position2[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + trinucleotide[chromosome1_dict[sorted_chromosome1[y]]]  + '\t' + mutation[chromosome1_dict[sorted_chromosome1[y]]])
    for fd in range(0, len(file_data)):
            #if genotype[chromosome1_dict[sorted_chromosome1[fd]]] == 'WT':
                file1a.write(file_data[fd])
               
                file1a.write('\n')
    
            
WriteChr('chrI',f1, f2, f3, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrII', f1, f2, f3, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrIII',  f1, f2, f3, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrIV', f1, f2, f3, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrIX', f1, f2, f3, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrV', f1, f2, f3, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrVI', f1, f2, f3, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrVII', f1, f2, f3, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrVIII', f1, f2, f3, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrX', f1, f2, f3, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXI', f1, f2, f3, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXII', f1, f2, f3, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXIII', f1, f2, f3, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXIV', f1, f2, f3, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand) 
WriteChr('chrXV', f1, f2, f3, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXVI', f1, f2, f3, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)


f1.close()
f2.close()
f3.close()

WriteAllChr('chrI',f7, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrII', f7, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrIII',  f7, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrIV', f7, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrIX', f7, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrV', f7, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrVI', f7, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrVII', f7, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrVIII', f7, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrX', f7, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrXI', f7, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrXII', f7, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrXIII', f7, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrXIV', f7, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand) 
WriteAllChr('chrXV', f7, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrXVI', f7, chromosome, position, position1, trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)

f7.close()

WriteChr('chrI',f4, f5, f6, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteChr('chrII', f4, f5, f6, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteChr('chrIII',  f4, f5, f6, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteChr('chrIV', f4, f5, f6, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteChr('chrIX', f4, f5, f6, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteChr('chrV', f4, f5, f6, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteChr('chrVI', f4, f5, f6, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteChr('chrVII', f4, f5, f6, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteChr('chrVIII', f4, f5, f6, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteChr('chrX', f4, f5, f6, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteChr('chrXI', f4, f5, f6, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteChr('chrXII', f4, f5, f6, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteChr('chrXIII', f4, f5, f6, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteChr('chrXIV', f4, f5, f6, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand) 
WriteChr('chrXV', f4, f5, f6, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteChr('chrXVI', f4, f5, f6, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)


f1.close()
f2.close()
f3.close()
f4.close()
f5.close()
f6.close()


WriteAllChr('chrI',f8, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteAllChr('chrII', f8, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteAllChr('chrIII',  f8, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteAllChr('chrIV', f8, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteAllChr('chrIX', f8, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteAllChr('chrV', f8, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteAllChr('chrVI', f8, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteAllChr('chrVII', f8, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteAllChr('chrVIII', f8, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteAllChr('chrX', f8, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteAllChr('chrXI', f8, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteAllChr('chrXII', f8, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteAllChr('chrXIII', f8, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteAllChr('chrXIV', f8, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand) 
WriteAllChr('chrXV', f8, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)
WriteAllChr('chrXVI', f8, WT_chromosome, WT_position, WT_position1, WT_trinucleotide, WT_allele, WT_mutation, WT_genotype, WT_mutation_type, WT_isolate, WT_strand)

f8.close()



#Makes bed files for single-base deletions, including just homopolymer mutations 
def MakeTriNucFile(file1, file2, file3, file4, lengths):
    file1.seek(0)
    new_data = file1.read()
    new_lines = [s.strip().split() for s in new_data.splitlines()]
    n = -1
    for line in new_lines:
            n = n + 1
            if len(line[6]) == 1:
                #if line[5] == '+':
                    file2.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' + line[5]  + '\t' + line[6])
                    file2.write('\n')
                    if lengths[n] >= 4:
                        file3.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] +  '\t' + line[4] + '\t' + line[5] + '\t' + line[6])
                        file3.write('\n')
                    else:
                        file4.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4]+ '\t' + line[5] + '\t' + line[6])
                        file4.write('\n')
               

#Makes bed files for single-base insertions, including just homopolymer mutations 
def MakeInsertionFile(file1, file2, file3, file4, lengths ):
    file1.seek(0)
    new_data = file1.read()
    new_lines = [s.strip().split() for s in new_data.splitlines()]
    n = -1
    for line in new_lines:
            n = n + 1
            if len(line[6]) == 1:
           
                file2.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' + line[5]  + '\t' + line[6])
                file2.write('\n')
              
                if lengths[n] >= 4:
                    file3.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' + line[5]  + '\t' + line[6])
                    file3.write('\n')
                else:
                    file4.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' + line[5]  + '\t' + line[6])
                    file4.write('\n')
               


f1a = open("rad30_insertions_sorted_tandem.bed")
f1b = open("rad30_singlebase_insertions_Trinuc.bed", 'w+')
f1c = open("rad30_singlebase_insertions_Trinuc_homopolymer.bed", 'w+')
f1d = open("rad30_singlebase_insertions_Trinuc_nonhomopolymer.bed", 'w+')


f2a = open("WT_insertionsv2_sorted_tandem.bed")
f2b = open("WT_singlebase_insertionsv2_Trinuc.bed", 'w+')
f2c = open("WT_singlebase_insertionsv2_Trinuc_homopolymer.bed", 'w+')
f2d = open("WT_singlebase_insertionsv2_Trinuc_nonhomopolymer.bed", 'w+')

#MakeInsertionFile(f1a, f1b)


#f1a.close()
#f1b.close()



a1 = open("rad30_deletions_sorted_tandem.bed")
a2 = open('rad30_singlebase_deletions_Trinuc.bed', 'w+')
a3 = open('rad30_singlebase_deletions_Trinuc_homopolymer.bed', 'w+')
a4 = open('rad30_singlebase_deletions_Trinuc_nonhomopolymer.bed', 'w+')
a5 = open('rad30_sorted_deletions_Pentanuc.bed', 'w+')
a6 = open('rad30_sorted_single_deletions_25nuc.bed', 'w+')
a7 = open("rad30_substitutionsv2_sorted_tandem.bed")
a8 = open('rad30_singlebase_substitutionsv2_Trinuc.bed', 'w+')
a9 = open('rad30_singlebase_substitutionsv2_Trinuc_homopolymer.bed', 'w+')
a10 = open('rad30_singlebase_substitutionsv2_Trinuc_nonhomopolymer.bed', 'w+')
a11 = open('rad30_sorted_substitutionsv2_Pentanuc.bed', 'w+')

# b1 = open("rad30_all_sorted_tandem.bed")
# b2 = open('rad30_singlebase_combined_Trinuc.bed', 'w+')
# b3 = open('rad30_singlebase_combined_Trinuc_dipyrimidine.bed', 'w+')
# b4 = open('rad30_singlebase_combined_Trinuc_nondipyrimidine.bed', 'w+')
# b5 = open('rad30_sorted_insertions_Pentanuc.bed', 'w+')
# b6 = open('rad30_sorted_singlebase_combined_25nuc.bed', 'w+')

b1 = open("rad30_allcombinedv2_sorted_tandem.bed")
b2 = open('rad30_singlebase_allcombined_Trinuc.bed', 'w+')
b3 = open('rad30_singlebase_allcombined_Trinuc_dipyrimidine.bed', 'w+')
b4 = open('rad30_singlebase_allcombined_Trinuc_nondipyrimidine.bed', 'w+')
b5 = open('rad30_sorted_insertions_Pentanuc.bed', 'w+')
b6 = open('rad30_sorted_singlebase_allcombined_25nuc.bed', 'w+')


c1 = open("WT_deletionsv2_sorted_tandem.bed")
c2 = open('WT_singlebase_deletionsv2_Trinuc.bed', 'w+')
c3 = open('WT_singlebase_deletionsv2_Trinuc_homopolymer.bed', 'w+')
c4 = open('WT_singlebase_deletionsv2_Trinuc_nonhomopolymer.bed', 'w+')

c6 = open('WT_sorted_deletionsv2_25nuc.bed', 'w+')
c7 = open("WT_substitutionsv2_sorted_tandem.bed")
c8 = open('WT_singlebase_substitutionsv2_Trinuc.bed', 'w+')
c9 = open('WT_singlebase_substitutionsv2_Trinuc_homopolymer.bed', 'w+')
c10 = open('WT_singlebase_substitutionsv2_Trinuc_nonhomopolymer.bed', 'w+')
c11 = open('WT_sorted_substitutionsv2_25nuc.bed', 'w+')

d1 = open("WT_insertionsv2_sorted_tandem.bed")
d6 = open('WT_sorted_insertionsv2_25nuc.bed', 'w+')
#e1 = open("WT_allv2_sorted_tandem.bed")
#e6 = open('WT_sorted_combinedv2_25nuc.bed', 'w+')

e1 = open("WT_allcombinedv2_sorted_tandem.bed")
e6 = open('WT_sorted_allcombinedv2_25nuc.bed', 'w+')





#MakeTriNucFile(a1, a2, a3, a4)
#MakeTriNucFile(b1, b2, b3, b4)   

#a1.close()
#a2.close()

#a3.close()
#a4.close()
#b1.close()
#b2.close()
b3.close()
b4.close()


def MakeMultiNucFile(file_1, file_2):  #Finds the 27 bases on either side of a mutation and the length of the homopolymer containing the mutation
    new_data = file_1.read()
    new_chromosome = []
    new_position = []
    new_position1 = []
    new_mutation = []
    new_allele = []
    new_trinucleotide = []
    new_strand = []
    new_isolate = []
    new_mutation_type = []
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
        new_mutation_type.append(line[8])

    #Add more nucleotides with the indel in the middle
    multinucleotide = []
    for pos in range(0, len(new_position)):
        if new_chromosome[pos] == 'chrI':
            sequence = GetSequence('chr1.txt')
            #if strand[pos] == '+':
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]: #Insertion
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                else: #Deletion
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    
            elif int(new_position[pos]) <= 27:   #Indels near telomeres 
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
        
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:  #Indels near telomeres
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrII':
            sequence = GetSequence('chr2.txt')
            
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[1]
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    
            elif int(new_position[pos]) <= 27:    
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
        
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                    
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrIII':
            sequence = GetSequence('chr3.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[1]
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    
            elif int(new_position[pos]) <= 27:    
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
        
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
            
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            


        elif new_chromosome[pos] == 'chrIV':
            sequence = GetSequence('chr4.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[1]
                       
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    
            elif int(new_position[pos]) <= 27:    
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
        
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrV':
            sequence = GetSequence('chr5.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[1]
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    
            elif int(new_position[pos]) <= 27:    
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
        
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrVI':
            sequence = GetSequence('chr6.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[1]
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    
            elif int(new_position[pos]) <= 27:    
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
        
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrVII':
            sequence = GetSequence('chr7.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[1]
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    
            elif int(new_position[pos]) <= 27:    
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
        
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chrVIII':
            sequence = GetSequence('chr8.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[1]
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    
            elif int(new_position[pos]) <= 27:    
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                       
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
        
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                       
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            
        elif new_chromosome[pos] == 'chrIX':
            sequence = GetSequence('chr9.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[1]
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    
            elif int(new_position[pos]) <= 27:    
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
        
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            
        elif new_chromosome[pos] == 'chrX':
            sequence = GetSequence('chr10.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[1]
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    
            elif int(new_position[pos]) <= 27:    
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                    
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
        
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            
        elif new_chromosome[pos] == 'chrXI':
            sequence = GetSequence('chr11.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[1]
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    
            elif int(new_position[pos]) <= 27:    
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
        
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            
        elif new_chromosome[pos] == 'chrXII':
            sequence = GetSequence('chr12.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[1]
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    
            elif int(new_position[pos]) <= 27:    
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
        
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            
        elif new_chromosome[pos] == 'chrXIII':
            sequence = GetSequence('chr13.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[1]
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    
            elif int(new_position[pos]) <= 27:    
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
        
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            
        elif new_chromosome[pos] == 'chrXIV':
            sequence = GetSequence('chr14.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[1]
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    
            elif int(new_position[pos]) <= 27:    
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
        
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            
        elif new_chromosome[pos] == 'chrXV':
            sequence = GetSequence('chr15.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[1]
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    
            elif int(new_position[pos]) <= 27:    
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                       
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
        
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            
        elif new_chromosome[pos] == 'chrXVI':
            sequence = GetSequence('chr16.txt')
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1, new_strand[pos])[1]
                       
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    
            elif int(new_position[pos]) <= 27:    
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
        
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, 27 - (len(sequence) - int(new_position[pos]))):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if '-' in new_allele[pos]:
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1, new_strand[pos])[1]
                    if lengths[len(lengths) - 1] >= 4:
                        
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                else:
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            
            
    
    for i in range(0, len(new_position)):
        #if int(lengths) >= 4:
        #if int(lengths[i]) >= 4 and new_allele[i] == '-' and len(new_mutation[i]) == 1:
        #if int(lengths[i]) < 4 and new_allele[i] == '-' and len(new_mutation[i]) == 1:
        #if len(new_allele[i]) == 1:
            if new_strand[i] == '-':
                file_2.write(new_chromosome[i] + '\t' + str(int(new_position[i]) - 1) + '\t'  +str(new_position[i]) + '\t' + new_trinucleotide[i] + '\t' + new_mutation[i] + '\t' + new_strand[i]  + '\t'+ multinucleotide[i] + '\t' + str(lengths[i]) + '\t' + new_isolate[i])
                file_2.write('\n')
            if new_strand[i] == '+':
                file_2.write(new_chromosome[i] + '\t' + str(int(new_position[i]) - 1)  + '\t' + str(new_position[i]) + '\t' + new_trinucleotide[i] + '\t' +  new_mutation[i] + '\t' + new_strand[i] + '\t' + multinucleotide[i] + '\t' + str(lengths[i]) + '\t' + new_isolate[i])
                file_2.write('\n')

    
    
    return multinucleotide, new_strand, new_chromosome, new_position, new_isolate, new_trinucleotide, new_mutation, new_allele, lengths, new_position1, new_mutation_type

sequencefile1 = open("rad30_single_deletion_sequences.txt", 'w+')
sequencefile2 = open("rad30_single_insertion_sequences.txt", 'w+')


def MakeSequenceFile(file1, file2, pentanucleotide, mutation, allele, strand):
    for p in range(0, len(pentanucleotide)):
        #if allele[p] == 'C' or allele[p] == 'G':
        #if allele[p] == 'T' or allele[p] == 'A':
            if '-' in mutation[p] and len(allele[p]) == 1:
                if strand[p] == '-':
                    file1.write(GetComplement(GetReverseComplement(pentanucleotide[p])))
                    file1.write('\n')
                else:
                    file1.write(pentanucleotide[p])
                    file1.write('\n')
        #elif mutation[p]  == 'C' or mutation[p] == 'G':
        #if mutation[p] == 'T' or mutation[p] == 'A':
            if '-' in allele[p]  and len(mutation[p]) == 1:
                if strand[p] == '-':
                    file2.write(GetComplement(GetReverseComplement(pentanucleotide[p])))
                    file2.write('\n')
                else:
                    file2.write(pentanucleotide[p])
                    file2.write('\n')

#Counts the numbers of mutations in homopolymers
def GetPolymerCounts(multinucleotide,  new_chromosome, new_position, new_isolate, new_trinucleotide, new_mutation, new_allele, lengths):
    homopolymers = {}
    homopolymer_list = {}
    deletion_numbers = {}
    insertion_numbers = {}
    homopolymer_deletions = 0
    homopolymer_insertions = 0
    polymer_counts = []
    for j in range(0, 4):
        polymer_counts.append({})
    for i in range (0, len(lengths)):
    
            key = lengths[i]
            
            #elif strand[i] == '-':
                #key = GetHomopolymerLength(GetReverseComplement(sequence), )[1]
            homopolymers[(new_isolate[i], new_chromosome[i], new_position[i])] = key
            if int(key) >= 4:
                #homopolymers[(new_isolate[i], new_chromosome[i], new_position[i])] = key
                if '-' in new_mutation[i]:
                #if new_mutation[i] == '-':
                    homopolymer_deletions = homopolymer_deletions + 1
                    if key in deletion_numbers.keys():
                        deletion_numbers[key] = deletion_numbers[key] + 1
                    else:
                        deletion_numbers[key] = 1
                elif '-' in new_allele[i]:
                #if new_allele[i] == '-':
                    homopolymer_insertions = homopolymer_insertions + 1
                    if key in insertion_numbers.keys():
                        insertion_numbers[key] = insertion_numbers[key] + 1
                    else:
                        insertion_numbers[key] = 1
                #homopolymer_list[key[1]] = [(new_chromosome[i], key[0][0], key[0][1])]
               
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





homopolymer_data = MakeMultiNucFile(a1, a5)
polymer_counts = GetPolymerCounts(homopolymer_data[0],  homopolymer_data[2], homopolymer_data[3], homopolymer_data[4], homopolymer_data[5], homopolymer_data[6], homopolymer_data[7], homopolymer_data[8])[0]
homopolymer1 = GetPolymerCounts(homopolymer_data[0],  homopolymer_data[2], homopolymer_data[3], homopolymer_data[4], homopolymer_data[5], homopolymer_data[6], homopolymer_data[7], homopolymer_data[8])[1]
homopolymer_data2 = MakeMultiNucFile(f1a, b5)
polymer_counts2 = GetPolymerCounts(homopolymer_data2[0],  homopolymer_data2[2], homopolymer_data2[3], homopolymer_data2[4], homopolymer_data2[5], homopolymer_data2[6], homopolymer_data2[7], homopolymer_data2[8])[0]
homopolymer2 = GetPolymerCounts(homopolymer_data2[0],  homopolymer_data2[2], homopolymer_data2[3], homopolymer_data2[4], homopolymer_data2[5], homopolymer_data2[6], homopolymer_data2[7], homopolymer_data2[8])[1]
homopolymer_data3 = MakeMultiNucFile(b1, b6)
polymer_counts3 = GetPolymerCounts(homopolymer_data3[0],  homopolymer_data3[2], homopolymer_data3[3], homopolymer_data3[4], homopolymer_data3[5], homopolymer_data3[6], homopolymer_data3[7], homopolymer_data3[8])[0]
homopolymer3 = GetPolymerCounts(homopolymer_data3[0],  homopolymer_data3[2], homopolymer_data3[3], homopolymer_data3[4], homopolymer_data3[5], homopolymer_data3[6], homopolymer_data3[7], homopolymer_data3[8])[1]
homopolymer_data7 = MakeMultiNucFile(a7, a11)
polymer_counts7 = GetPolymerCounts(homopolymer_data7[0],  homopolymer_data7[2], homopolymer_data7[3], homopolymer_data7[4], homopolymer_data7[5], homopolymer_data7[6], homopolymer_data7[7], homopolymer_data7[8])[0]
homopolymer7 = GetPolymerCounts(homopolymer_data7[0],  homopolymer_data7[2], homopolymer_data7[3], homopolymer_data7[4], homopolymer_data7[5], homopolymer_data7[6], homopolymer_data7[7], homopolymer_data7[8])[1]


rad30_data = homopolymer_data3
rad30_isolate = rad30_data[4]
rad30_chromosome = rad30_data[2]
rad30_position1 = rad30_data[9]
rad30_position = rad30_data[3]
rad30_trinucleotide = rad30_data[5]
rad30_mutation = rad30_data[6]
rad30_strand = rad30_data[1]
rad30_multinucleotide = rad30_data[0] 
rad30_lengths = rad30_data[8]
rad30_mutationtype = rad30_data[10]
rad30_allele = rad30_data[7]


rad30_del_lengths = homopolymer_data[8]


rad30_ins_lengths = homopolymer_data2[8]


rad30_sub_lengths = homopolymer_data7[8]

rad30_isolatelist = {}
for i1 in range(0, len(rad30_isolate)):

    if rad30_isolate[i1] not in rad30_isolatelist.keys():
            rad30_isolatelist[rad30_isolate[i1]] = [(rad30_position[i1],  str(rad30_chromosome[i1]), rad30_trinucleotide[i1], rad30_mutation[i1], rad30_strand[i1], rad30_multinucleotide[i1], rad30_isolate[i1], rad30_position1[i1], rad30_lengths[i1], rad30_mutationtype[i1], rad30_allele[i1])]
    else:
        rad30_isolatelist[rad30_isolate[i1]] .append((rad30_position[i1],  str(rad30_chromosome[i1]), rad30_trinucleotide[i1], rad30_mutation[i1], rad30_strand[i1], rad30_multinucleotide[i1], rad30_isolate[i1], rad30_position1[i1], rad30_lengths[i1], rad30_mutationtype[i1], rad30_allele[i1]))


homopolymer_data4 = MakeMultiNucFile(c1, c6)
polymer_counts4 = GetPolymerCounts(homopolymer_data4[0],  homopolymer_data4[2], homopolymer_data4[3], homopolymer_data4[4], homopolymer_data4[5], homopolymer_data4[6], homopolymer_data4[7], homopolymer_data4[8])[0]
homopolymer4 = GetPolymerCounts(homopolymer_data4[0],  homopolymer_data4[2], homopolymer_data4[3], homopolymer_data4[4], homopolymer_data4[5], homopolymer_data4[6], homopolymer_data4[7], homopolymer_data4[8])[1]
homopolymer_data5 = MakeMultiNucFile(d1, d6)
polymer_counts5 = GetPolymerCounts(homopolymer_data5[0],  homopolymer_data5[2], homopolymer_data5[3], homopolymer_data5[4], homopolymer_data5[5], homopolymer_data5[6], homopolymer_data5[7], homopolymer_data5[8])[0]
homopolymer5 = GetPolymerCounts(homopolymer_data5[0],  homopolymer_data5[2], homopolymer_data5[3], homopolymer_data5[4], homopolymer_data5[5], homopolymer_data5[6], homopolymer_data5[7], homopolymer_data5[8])[1]
homopolymer_data6 = MakeMultiNucFile(e1, e6)
polymer_counts6 = GetPolymerCounts(homopolymer_data6[0],  homopolymer_data6[2], homopolymer_data6[3], homopolymer_data6[4], homopolymer_data6[5], homopolymer_data6[6], homopolymer_data6[7], homopolymer_data6[8])[0]
homopolymer6 = GetPolymerCounts(homopolymer_data6[0],  homopolymer_data6[2], homopolymer_data6[3], homopolymer_data6[4], homopolymer_data6[5], homopolymer_data6[6], homopolymer_data6[7], homopolymer_data6[8])[1]
homopolymer_data8 = MakeMultiNucFile(c7, c11)
polymer_counts8 = GetPolymerCounts(homopolymer_data8[0],  homopolymer_data8[2], homopolymer_data8[3], homopolymer_data8[4], homopolymer_data8[5], homopolymer_data8[6], homopolymer_data8[7], homopolymer_data8[8])[0]
homopolymer8 = GetPolymerCounts(homopolymer_data8[0],  homopolymer_data8[2], homopolymer_data8[3], homopolymer_data8[4], homopolymer_data8[5], homopolymer_data8[6], homopolymer_data8[7], homopolymer_data8[8])[1]


WT_data = homopolymer_data6
WT_isolate2 = WT_data[4]
WT_chromosome2 = WT_data[2]
WT_position1_2 = WT_data[9]
WT_position2 = WT_data[3]
WT_trinucleotide2 = WT_data[5]
WT_mutation2 = WT_data[6]
WT_strand2 = WT_data[1]
WT_multinucleotide = WT_data[0] 
WT_lengths = WT_data[8]
WT_mutationtype2 = WT_data[10]
WT_allele2 = WT_data[7]


WT_del_lengths = homopolymer_data4[8]


WT_ins_lengths = homopolymer_data5[8]


WT_sub_lengths = homopolymer_data8[8]

WT_isolatelist2 = {}
for i1 in range(0, len(WT_isolate2)):

    if WT_isolate2[i1] not in WT_isolatelist2.keys():
            WT_isolatelist2[WT_isolate2[i1]] = [(WT_position2[i1],  str(WT_chromosome2[i1]), WT_trinucleotide2[i1], WT_mutation2[i1], WT_strand2[i1], WT_multinucleotide[i1], WT_isolate2[i1], WT_position1_2[i1], WT_lengths[i1], WT_mutationtype2[i1], WT_allele2[i1])]
    else:
        WT_isolatelist2[WT_isolate2[i1]] .append((WT_position2[i1],  str(WT_chromosome2[i1]), WT_trinucleotide2[i1], WT_mutation2[i1], WT_strand2[i1], WT_multinucleotide[i1], WT_isolate2[i1], WT_position1_2[i1], WT_lengths[i1], WT_mutationtype2[i1], WT_allele2[i1]))




MakeSequenceFile(sequencefile1, sequencefile2, rad30_multinucleotide, mutation, allele, strand)


sequencefile1.close()
sequencefile2.close()

MakeTriNucFile(a1, a2, a3, a4,rad30_del_lengths )
MakeInsertionFile(f1a, f1b, f1c, f1d, rad30_ins_lengths)
MakeTriNucFile(c1, c2, c3, c4, WT_del_lengths )
MakeInsertionFile(f2a, f2b, f2c, f2d, WT_ins_lengths)
MakeTriNucFile(a7, a8, a9, a10,rad30_sub_lengths )
MakeTriNucFile(c7, c8, c9, c10, WT_sub_lengths )

f1b.close()
f1a.close()
f1c.close()
f1d.close()
a1.close()
a2.close()
a3.close()
a4.close()
a5.close()
a6.close()
a7.close()
a8.close()
a9.close()
a10.close()
a11.close()
b1.close()
b2.close()
b3.close()
b4.close()
b6.close()
c1.close()
c6.close()
c7.close()
c8.close()
c9.close()
c10.close()
c11.close()
d1.close()
d6.close()
e1.close()
e6.close()
f2a.close()
f2b.close()
f2c.close()
f2d.close()
c2.close()
c3.close()
c4.close()


homopolymer_file = open('homopolymer_file', 'w+')



#Obtain the total number of bases in the yeast genome in sequences of each homopolymer length
def NormalizeHomopolymer(sequences, lengths, file):
    homopolymer_dict = {}
    C_dict = {}
    T_dict = {}
    homopolymer_counts = {}
    chromosome_sequences = {}
    chromosome_sequences['chr1'] = sequences[0]
    chromosome_sequences['chr2'] = sequences[1]
    chromosome_sequences['chr3'] = sequences[2]
    chromosome_sequences['chr4'] = sequences[3]
    chromosome_sequences['chr5'] = sequences[4]
    chromosome_sequences['chr6'] = sequences[5]
    chromosome_sequences['chr7'] = sequences[6]
    chromosome_sequences['chr8'] = sequences[7]
    chromosome_sequences['chr9'] = sequences[8]
    chromosome_sequences['chr10'] = sequences[9]
    chromosome_sequences['chr11'] = sequences[10]
    chromosome_sequences['chr12'] = sequences[11]
    chromosome_sequences['chr13'] = sequences[12]
    chromosome_sequences['chr14'] = sequences[13]
    chromosome_sequences['chr15'] = sequences[14]
    chromosome_sequences['chr16'] = sequences[15]
  
    


   
    homopolymer_dict['11-15'] = 0
    homopolymer_dict['16-20'] = 0
    homopolymer_dict['21+'] = 0
    C_dict['11-15'] = 0
    C_dict['16-20'] = 0
    C_dict['21+'] = 0
    T_dict['11-15'] = 0
    T_dict['16-20'] = 0
    T_dict['21+'] = 0

    for l in range(0, len(lengths)):
        if lengths[l] not in homopolymer_counts.keys():
            homopolymer_counts[lengths[l]] = 1
            homopolymer_dict[lengths[l]] = 0
        else:
            homopolymer_counts[lengths[l]] = homopolymer_counts[lengths[l]] + 1

    homopolymer_dict[1] = 0
    C_dict[1] = 0
    T_dict[1] = 0
    for key in chromosome_sequences.keys():
        coordinates = []
        sequence = chromosome_sequences[key]
        if sequence[len(sequence) - 1] != sequence[len(sequence) - 2]:
            homopolymer_dict[1] = homopolymer_dict[1] + 1
            if sequence[len(sequence) - 1] == 'C' or sequence[len(sequence) - 1] == 'G':
                C_dict[1] = C_dict[1] + 1
            elif sequence[len(sequence) - 1] == 'T' or sequence[len(sequence) - 1] == 'A':
                T_dict[1] = T_dict[1] + 1
        b = 0
        while b < len(sequence) - 1:
        #for b in range(0, len(sequence) - 1):
            polymer_range = {}
            length = 1
            if sequence[b] != 'N':# and b not in coordinates:
                length = 1
                i = 1
                polymer_range[0] = b
                while sequence[b] == sequence[b + i]:
                    length = length + 1
                    polymer_range[1] = b + i
                    if (b + i) < len(sequence) - 1:
                        i = i + 1
                    else:
                        break
              
                polymer_count = length
                if length >= 4:
                    if (b + 1, b + length) not in coordinates:
                        coordinates.append((b + 1, b + length))

                    #for a in range(b + 1, b + length + 1):
                                #if a not in coordinates:
                                    #coordinates.append(a)
                
                        #polymer_range = number[1]
                #if length > 1:
                if polymer_count in homopolymer_dict.keys():
                            #homopolymer_dict[polymer_count] = homopolymer_dict[polymer_count] + 1
                            #for i in range (polymer_range[0], polymer_range[1] + 1):
                            homopolymer_dict[polymer_count] = homopolymer_dict[polymer_count] + length
                                #coordinates.append(i)
                            if sequence[b] == 'C' or sequence[b] == 'G':
                                if polymer_count in C_dict.keys():
                                    #C_dict[polymer_count] = C_dict[polymer_count] + 1
                                    #for k in range (polymer_range[0], polymer_range[1] + 1):
                                    C_dict[polymer_count] = C_dict[polymer_count] + length
                                else:
                                    C_dict[polymer_count] = length
                                    #for l in range (polymer_range[0], polymer_range[1] + 1):
                                        #C_dict[polymer_count] = C_dict[polymer_count] + 1
                            if sequence[b] == 'T' or sequence[b] == 'A':
                                if polymer_count in T_dict.keys():
                                    #T_dict[polymer_count] = T_dict[polymer_count] + 1
                                    #for k in range (polymer_range[0], polymer_range[1] + 1):
                                    T_dict[polymer_count] = T_dict[polymer_count] + length
                                else:
                                    T_dict[polymer_count] = length
                                    #for l in range (polymer_range[0], polymer_range[1] + 1):
                                        #T_dict[polymer_count] = T_dict[polymer_count] + 1

    

                else:
                            homopolymer_dict[polymer_count] = length
                            #for j in range (polymer_range[0], polymer_range[1] + 1):
                                #homopolymer_dict[polymer_count] = homopolymer_dict[polymer_count] + 1
                                #coordinates.append(j)
                            if sequence[b] == 'C' or sequence[b] == 'G':
                                C_dict[polymer_count] = length
                                #for m in range (polymer_range[0], polymer_range[1] + 1):
                                    #C_dict[polymer_count] = C_dict[polymer_count] + 1
                            if sequence[b] == 'T' or sequence[b] == 'A':
                                T_dict[polymer_count] = length
                                #for n in range (polymer_range[0], polymer_range[1] + 1):
                if polymer_count >= 11:
                    if polymer_count >= 21:
                        homopolymer_dict['21+'] = homopolymer_dict['21+'] + length
                        if sequence[b] == 'C' or sequence[b] == 'G':
                            C_dict['21+'] = C_dict['21+'] + length
                        if sequence[b] == 'T' or sequence[b] == 'A':
                            T_dict['21+'] = T_dict['21+'] + length
                    elif polymer_count >= 16:
                        homopolymer_dict['16-20'] = homopolymer_dict['16-20'] + length
                        if sequence[b] == 'C' or sequence[b] == 'G':
                            C_dict['16-20'] = C_dict['16-20'] + length
                        if sequence[b] == 'T' or sequence[b] == 'A':
                            T_dict['16-20'] = T_dict['16-20'] + length
                    else:
                        homopolymer_dict['11-15'] = homopolymer_dict['11-15'] + length
                        if sequence[b] == 'C' or sequence[b] == 'G':
                            C_dict['11-15'] = C_dict['11-15'] + length
                        if sequence[b] == 'T' or sequence[b] == 'A':
                            T_dict['11-15'] = T_dict['11-15'] + length                    #T_dict[polymer_count] = T_dict[polymer_count] + 1
               
            b = b + length
    

                                

                         
    for key in homopolymer_dict.keys():
        #homopolymer_dict[key] = homopolymer_counts[key]/homopolymer_dict[key]
        file.write(str(key) + ':' + str(homopolymer_dict[key]))
        file.write('\n')
    for key in C_dict.keys():
        #homopolymer_dict[key] = homopolymer_counts[key]/homopolymer_dict[key]
        file.write(str(key) + ':' + str(C_dict[key]))
        file.write('\n')
    for key in T_dict.keys():
        #homopolymer_dict[key] = homopolymer_counts[key]/homopolymer_dict[key]
        file.write(str(key) + ':' + str(T_dict[key]))
        file.write('\n')
    return homopolymer_dict, C_dict, T_dict, coordinates


numbers = NormalizeHomopolymer(sequences, rad30_lengths, homopolymer_file)
C_dict = numbers[1]
T_dict = numbers[2]

rad30_file = open('rad30_del', 'w+')
rad30_file2 = open('rad30_ins', 'w+')
WT_file = open('WT_del', 'w+')
WT_file2 = open('WT_ins', 'w+')

def PrintFrequencies( polymer_counts, C_dict, T_dict, file):
    for a in range(0, len(polymer_counts)):
        if a %2 == 0:
            for key in polymer_counts[a]:
                file.write(str(key) + " :" + str(polymer_counts[a][key]/C_dict[key]))
                file.write("\n") 
        else:
            for key in polymer_counts[a]:
                file.write(str(key) + " :" + str(polymer_counts[a][key]/T_dict[key])) 
                file.write("\n") 


PrintFrequencies( polymer_counts, C_dict, T_dict, rad30_file)

PrintFrequencies( polymer_counts2, C_dict, T_dict, rad30_file2)


PrintFrequencies( polymer_counts4, C_dict, T_dict, WT_file)


PrintFrequencies( polymer_counts5, C_dict, T_dict, WT_file2)


homopolymer_file.close()
rad30_file.close()
rad30_file2.close()
WT_file.close()
WT_file2.close()



def index(tuple):
    return tuple[1]


#Find indels that are <= 10 base pairs from another indel in the same isolate
    lengths = {}
def FindClusters( isolate_list, list2, homopolymers):
    lengths = {}
    clusters = {}
    nonclusters = {}
    total_clusters = 0
    indel_clusters = 0
    list2_key = {}
    for key in isolate_list.keys():
        #if strain in key[1]:
            sorted_isolate_list = sorted(isolate_list[key])
            sorted_list = sorted(sorted_isolate_list, key = index)
            if key in list2.keys():
                list2_key = list2[key]
            new_list2 = {}
            for i in range(0, len(list2_key)):
                for j in range(0, len(sorted_list)):
                    if sorted_list[j][1] == list2_key[i][1] and sorted_list[j][0] == list2_key[i][0]:
                        new_list2[j] = list2_key[i]

            if sorted_list[0][1] == sorted_list[1][1] and (int(sorted_list[1][0]) <= int((sorted_list[ 0][0])) + 10):# and (int(sorted_list[1][0]) != int((sorted_list[ 0][0]))):
                if (sorted_list[0][6], sorted_list[0][1], sorted_list[0][0]) not in clusters.keys() and (sorted_list[0][5] == "DIP" or sorted_list[0][5] == "Deletion" or sorted_list[0][5] == "Insertion"):
                #if (sorted_list[0][6], sorted_list[0][1], sorted_list[0][0]) not in clusters.keys() and (sorted_list[0][5] == "SNV" or sorted_list[0][5] == "SNP"):    
                    clusters[(sorted_list[0][6], sorted_list[0][1], sorted_list[0][0])] = [sorted_list[0]]
                    #total_clusters = total_clusters + 1
                    #indel_clusters = indel_clusters + 1
                #elif isolate_list[key][0] not in clusters.keys():
                    #total_clusters = total_clusters + 1

                if (sorted_list[1][6], sorted_list[1][1], sorted_list[1][0]) not in clusters.keys() and (sorted_list[1][5] == "DIP" or sorted_list[1][5] == "Deletion" or sorted_list[1][5] == "Insertion"):
                #if (sorted_list[1][6], sorted_list[1][1], sorted_list[1][0]) not in clusters.keys() and (sorted_list[1][5] == "SNV" or sorted_list[1][5] == "SNP"):
                    clusters[(sorted_list[1][6], sorted_list[1][1], sorted_list[1][0])] = [sorted_list[1]]
                    #total_clusters = total_clusters + 1
                    #indel_clusters = indel_clusters + 1
                #elif isolate_list[key][1] not in clusters.keys():
                    #total_clusters = total_clusters + 1
            else:#elif (int(sorted_list[1][0]) != int((sorted_list[ 0][0]))):
                if ('-' not in sorted_list[0][3] and sorted_list[0][5] == "DIP") or sorted_list[0][5] == "Insertion":
                    if sorted_list[0][4] == '-':
                        DNA_strand = '-'
                        end_position = GetHomopolymerLengthIns(GetReverseComplement(new_list2[0][5]), int(new_list2[0][7]), int(new_list2[0][0]), int(len(new_list2[0][5])/2) + 1, DNA_strand)[1][1] 
                    else:
                        DNA_strand = '+'
                        end_position = GetHomopolymerLengthIns(new_list2[0][5], int(new_list2[0][7]), int(new_list2[0][0]), int(len(new_list2[0][5])/2) + 1, DNA_strand)[1][1] 
                else:#if (sorted_list[0][5] == "DIP" or sorted_list[0][5] == "Complex DIP" or sorted_list[0][5] == "Deletion" or sorted_list[0][5] == "Insertion") and 0 in new_list2.keys():
                    if sorted_list[0][4] == '-':
                        end_position = GetHomopolymerLength(GetReverseComplement(new_list2[0][5]), int(new_list2[0][7]), int(new_list2[0][0]), int(len(new_list2[0][5])/2) + 1)[1][1] 
                    else:
                        end_position = GetHomopolymerLength(new_list2[0][5], int(new_list2[0][7]), int(new_list2[0][0]), int(len(new_list2[0][5])/2) + 1)[1][1] 
                if int(sorted_list[1][0]) <= (end_position + 10) and sorted_list[0][1] == sorted_list[1][1]:
                    if (sorted_list[0][5] == "DIP" or sorted_list[0][5] == "Deletion" or sorted_list[0][5] == "Insertion"):
                    #if (sorted_list[0][5] == "SNV" or sorted_list[0][5] == "SNP") and 0 in new_list2.keys():
                        if (sorted_list[0][6], sorted_list[0][1], sorted_list[0][0]) not in clusters.keys():
                            
                            clusters[(sorted_list[0][6], sorted_list[0][1], sorted_list[0][0])] = [sorted_list[0]]
                    if (sorted_list[1][5] == "DIP" or sorted_list[1][5] == "Deletion" or sorted_list[1][5] == "Insertion"):
                    #if (sorted_list[1][5] == "SNV" or sorted_list[1][5] == "SNP") and 1 in new_list2.keys():
                        if (sorted_list[1][6], sorted_list[1][1], sorted_list[1][0]) not in clusters.keys():
                            clusters[(sorted_list[1][6], sorted_list[1][1], sorted_list[1][0])] = [sorted_list[1]]
                if ('-' not in sorted_list[1][3] and sorted_list[0][5] == "DIP") or sorted_list[1][5] == "Insertion":
                    if sorted_list[1][4] == '-':
                        DNA_strand = '-'
                        start_position = GetHomopolymerLengthIns(GetReverseComplement(new_list2[1][5]), int(new_list2[1][7]), int(new_list2[1][0]), int(len(new_list2[1][5])/2) + 1, DNA_strand)[1][0]
                    else:
                        DNA_strand = '+'
                        start_position = GetHomopolymerLengthIns(new_list2[1][5], int(new_list2[1][7]), int(new_list2[1][0]), int(len(new_list2[1][5])/2) + 1, DNA_strand)[1][0] 
                else:    
                    if sorted_list[1][4] == '-':
                        start_position = GetHomopolymerLength(GetReverseComplement(new_list2[1][5]), int(new_list2[1][7]), int(new_list2[1][0]), int(len(new_list2[1][5])/2) + 1)[1][0] 
                    else:
                        start_position = GetHomopolymerLength(new_list2[1][5], int(new_list2[1][7]), int(new_list2[1][0]), int(len(new_list2[1][5])/2) + 1)[1][0] 
                if start_position <= (int(sorted_list[0][0]) + 10) and sorted_list[0][1] == sorted_list[1][1] :
                    if (sorted_list[1][5] == "DIP" or sorted_list[1][5] == "Deletion" or sorted_list[1][5] == "Insertion"):
                    #if (sorted_list[1][5] == "SNV" or sorted_list[1][5] == "SNP") and 1 in new_list2.keys():
                        if (sorted_list[1][6], sorted_list[1][1], sorted_list[1][0]) not in clusters.keys():
                            clusters[(sorted_list[1][6], sorted_list[1][1], sorted_list[1][0])] = [sorted_list[1]]
                    if (sorted_list[0][5] == "DIP" or sorted_list[0][5] == "Deletion" or sorted_list[0][5] == "Insertion"):
                    #if (sorted_list[0][5] == "SNV" or sorted_list[0][5] == "SNP") and 0 in new_list2.keys():
                        if (sorted_list[0][6], sorted_list[0][1], sorted_list[0][0]) not in clusters.keys():
                            clusters[(sorted_list[0][6], sorted_list[0][1], sorted_list[0][0])] = [sorted_list[0]]
            #(sorted_list[0][6], sorted_list[0][1], sorted_list[0][0]) not in clusters.keys() and (sorted_list[0][5] == "DIP" or sorted_list[0][5] == "Complex DIP" or sorted_list[0][5] == "Deletion" or sorted_list[0][5] == "Insertion"):
                #if sorted_list[0][1] == sorted_list[1][1] and (int(sorted_list[1][0]) <= int((sorted_list[0 ][0])) + 10 + (homopolymers[(sorted_list[0][6], int(sorted_list[0][0]), sorted_list[0][1])][0][1] - sorted_list[0][0])):
            #elif sorted_list[0][1] == sorted_list[1][1] and (sorted_list[0][6], int(sorted_list[0][0]), sorted_list[0][1]) in homopolymers.keys() and (int(sorted_list[1][0]) <= int((sorted_list[0 ][0])) + 10 + (homopolymers[(sorted_list[key][0][6], int(sorted_list[key][0][0]), sorted_list[key][0][1])][0][1] - sorted_list[key][0][0])):
                
                   
            
            #else:
                #if isolate_list[key][0] not in nonclusters.keys() and (isolate_list[key][0][7] == "DIP" or isolate_list[key][0][7] == "Complex DIP" or isolate_list[key][0][7] == "Deletion" or isolate_list[key][0][7] == "Insertion"):
                    #nonclusters[(isolate_list[key][0][6], isolate_list[key][0][1], isolate_list[key][0][0])] = [isolate_list[key][0]]

            if sorted_list[len (sorted_list) - 2][1] == sorted_list[len (sorted_list) - 1][1] and (int(sorted_list[len (sorted_list) - 1][0]) <= int((sorted_list[ len (sorted_list) - 2][0])) + 10):# and (int(sorted_list[len (sorted_list) - 1][0]) != int((sorted_list[ len (sorted_list) - 2][0]))):
                if (sorted_list[len(sorted_list) - 1][6], sorted_list[len(sorted_list) - 1][1], sorted_list[len(sorted_list) - 1][0]) not in clusters.keys() and (sorted_list[len(sorted_list) - 1][5] == "DIP" or sorted_list[len(sorted_list) - 1][5] == "Deletion" or sorted_list[len(sorted_list) - 1][5] == "Insertion"):
                #if (sorted_list[len(sorted_list) - 1][6], sorted_list[len(sorted_list) - 1][1], sorted_list[len(sorted_list) - 1][0]) not in clusters.keys() and (sorted_list[len(sorted_list) - 1][5] == "SNV" or sorted_list[len(sorted_list) - 1][5] == "SNP"):    
                    clusters[(sorted_list[len(sorted_list) - 1][6], sorted_list[len(sorted_list) - 1][1], sorted_list[len(sorted_list) - 1][0])] = [sorted_list[len(sorted_list) - 1]]
                   
                if (sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][0]) not in clusters.keys() and (sorted_list[len(sorted_list) - 2][5] == "DIP" or sorted_list[len(sorted_list) - 2][5] == "Deletion" or sorted_list[len(sorted_list) - 2][5] == "Insertion"):
                #if (sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][0]) not in clusters.keys() and (sorted_list[len(sorted_list) - 2][5] == "SNV" or sorted_list[len(sorted_list) - 2][5] == "SNP"): 
                    clusters[(sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][0])] = [sorted_list[len(sorted_list) - 2]]
                   
            else:#elif (int(sorted_list[len(sorted_list) - 2][0]) != int((sorted_list[ len(sorted_list) - 1][0]))):#if strand[0] == '+':
                if (sorted_list[len(sorted_list) - 2][5] == "DIP" and '-' not in sorted_list[len(sorted_list) - 2][3]) or sorted_list[len(sorted_list) - 2][5] == "Insertion":
                    if sorted_list[len(sorted_list) - 2][4] == '-':
                        DNA_strand = "-"
                        end_position = GetHomopolymerLengthIns(GetReverseComplement(new_list2[len(sorted_list) - 2][5]), int(new_list2[len(sorted_list) - 2][7]), int(new_list2[len(sorted_list) - 2][0]), int(len(new_list2[len(sorted_list) - 2][5])/2) + 1, DNA_strand)[1][1] 
                    else:
                        DNA_strand = "+"
                        end_position = GetHomopolymerLengthIns(new_list2[len(sorted_list) - 2][5], int(new_list2[len(sorted_list) - 2][7]), int(new_list2[len(sorted_list) - 2][0]), int(len(new_list2[len(sorted_list) - 2][5])/2) + 1, DNA_strand)[1][1] 
                else:    
                    if sorted_list[len(sorted_list) - 2][4] == '-':
                        end_position = GetHomopolymerLength(GetReverseComplement(new_list2[len(sorted_list) - 2][5]), int(new_list2[len(sorted_list) - 2][7]), int(new_list2[len(sorted_list) - 2][0]), int(len(new_list2[len(sorted_list) - 2][5])/2) + 1)[1][1] 
                    else:
                        end_position = GetHomopolymerLength(new_list2[len(sorted_list) - 2][5], int(new_list2[len(sorted_list) - 2][7]), int(new_list2[len(sorted_list) - 2][0]), int(len(new_list2[len(sorted_list) - 2][5])/2) + 1)[1][1] 
                if int(sorted_list[len(sorted_list) - 1][0]) <= (end_position + 10) and sorted_list[len(sorted_list) - 1][1] == sorted_list[len(sorted_list) - 2][1]:
                    if (sorted_list[len(sorted_list) - 2][5] == "DIP" or sorted_list[len(sorted_list) - 2][5] == "Deletion" or sorted_list[len(sorted_list) - 2][5] == "Insertion"):
                    #if (sorted_list[len(sorted_list) - 2][5] == "SNV" or sorted_list[len(sorted_list) - 2][5] == "SNP") and len(sorted_list) - 2 in new_list2.keys():
                        if (sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][6]) not in clusters.keys():
                            clusters[(sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][6])] = [sorted_list[len(sorted_list) - 2]]
                    if (sorted_list[len(sorted_list) - 1][5] == "DIP" or sorted_list[len(sorted_list) - 1][5] == "Deletion" or sorted_list[len(sorted_list) - 1][5] == "Insertion"):  
                    #if (sorted_list[len(sorted_list) - 1][5] == "SNV" or sorted_list[len(sorted_list) - 1][5] == "SNP") and len(sorted_list) - 1 in new_list2.keys(): 
                        if (sorted_list[len(sorted_list) - 1][6], sorted_list[len(sorted_list) - 1][1], sorted_list[len(sorted_list) - 1][6]) not in clusters.keys():
                            clusters[(sorted_list[len(sorted_list) - 1][6], sorted_list[len(sorted_list) - 1][1], sorted_list[len(sorted_list) - 1][0])] = [sorted_list[len(sorted_list) - 1]]
                if (sorted_list[len(sorted_list) - 1][5] == "DIP" and '-' not in sorted_list[len(sorted_list) - 1][3]) or sorted_list[len(sorted_list) - 1][5] == "Insertion": 
                    if sorted_list[len(sorted_list) - 1][4] == '-':
                        DNA_strand = "-"
                        start_position = GetHomopolymerLengthIns(GetReverseComplement(new_list2[len(sorted_list) - 1][5]), int(new_list2[len(sorted_list) - 1][7]), int(new_list2[len(sorted_list) - 1][0]), int(len(new_list2[len(sorted_list) - 1][5])/2) + 1, DNA_strand)[1][0]
                    else:
                        DNA_strand = "+"
                        start_position = GetHomopolymerLengthIns(new_list2[len(sorted_list) - 1][5], int(new_list2[len(sorted_list) - 1][7]), int(new_list2[len(sorted_list) - 1][0]), int(len(new_list2[len(sorted_list) - 1][5])/2) + 1, DNA_strand)[1][0]    
                else:    
                    if sorted_list[len(sorted_list) - 1][4] == '-':
                        start_position = GetHomopolymerLength(GetReverseComplement(new_list2[len(sorted_list) - 1][5]), int(new_list2[len(sorted_list) - 1][7]), int(new_list2[len(sorted_list) - 1][0]), int(len(new_list2[len(sorted_list) - 1][5])/2) + 1)[1][0]
                    else:
                        start_position = GetHomopolymerLength(new_list2[len(sorted_list) - 1][5], int(new_list2[len(sorted_list) - 1][7]), int(new_list2[len(sorted_list) - 1][0]), int(len(new_list2[len(sorted_list) - 1][5])/2) + 1)[1][0]
                if start_position <= int((sorted_list[len(sorted_list) - 2][0] + 10)) and sorted_list[len(sorted_list) - 2][1] == sorted_list[len(sorted_list) - 1][1]:
                    if (sorted_list[len(sorted_list) - 1][5] == "DIP" or sorted_list[len(sorted_list) - 1][5] == "Deletion" or sorted_list[len(sorted_list) - 1][5] == "Insertion"):
                    #if (sorted_list[len(sorted_list) - 1][5] == "SNV" or sorted_list[len(sorted_list) - 1][5] == "SNP") and len(sorted_list) - 1 in new_list2.keys():
                        if (sorted_list[len(sorted_list) - 1][6], sorted_list[len(sorted_list) - 1][1], sorted_list[len(sorted_list) - 1][0]) not in clusters.keys():
                            clusters[(sorted_list[len(sorted_list) - 1][6], sorted_list[len(sorted_list) - 1][1], sorted_list[len(sorted_list) - 1][0])] = [sorted_list[len(sorted_list) - 1]]
                    if (sorted_list[len(sorted_list) - 2][5] == "DIP" or sorted_list[len(sorted_list) - 2][5] == "Deletion" or sorted_list[len(sorted_list) - 2][5] == "Insertion"):
                    #if (sorted_list[len(sorted_list) - 2][5] == "SNV" or sorted_list[len(sorted_list) - 2][5] == "SNP") and len(sorted_list) - 2 in new_list2.keys(): 
                        if (sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][0]) not in clusters.keys():
                            clusters[(sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][0])] = [sorted_list[len(sorted_list) - 2]]
                   
            #else:
                #if isolate_list[key][len(isolate_list[key]) - 1] not in nonclusters.keys() and (isolate_list[key][len(isolate_list[key]) - 1][7] == "DIP" or isolate_list[key][0][7] == "Complex DIP" or isolate_list[key][0][7] == "Deletion" or isolate_list[key][0][7] == "Insertion"):
                    #nonclusters[(isolate_list[key][len(isolate_list[key]) - 1][6], isolate_list[key][len(isolate_list[key]) - 1][1], isolate_list[key][len(isolate_list[key]) - 1][0])] = [isolate_list[key][len(isolate_list[key]) - 1]]

            for j in range (1,len (sorted_list) - 1):
                    
                    if (sorted_list[j][1] == sorted_list[j - 1][1] and (int(sorted_list[j][0]) <= int((sorted_list[j - 1][0])) + 10)) or (sorted_list[j + 1][1] == sorted_list[j][1] and (int(sorted_list[j + 1][0]) <= int((sorted_list[j][0])) + 10)):
                        #or ((isolate_list[key][j + n][6], int(isolate_list[key][j + n][0]), isolate_list[key][j + n][1]) in homopolymers.keys() and (int(isolate_list[key][j + n][0]) <= int((isolate_list[key][j + n - 1][0])) + 10 + homopolymers[key])) : 
                            if (sorted_list[j][6], sorted_list[j][1], sorted_list[j][0]) not in clusters.keys() and (sorted_list[j][5] == "DIP" or sorted_list[j][5] == "Deletion" or sorted_list[j][5] == "Insertion"):
                            #if (sorted_list[j][6], sorted_list[j][1], sorted_list[j][0]) not in clusters.keys() and (sorted_list[j][5] == "SNV" or sorted_list[j][5] == "SNP"):    
                                #if j == 0 or isolate_list[key][j - n] not in clusters.keys():
                                clusters[(sorted_list[j][6], sorted_list[j][1], sorted_list[j][0])] = [sorted_list[j]]
                                #total_clusters = total_clusters + 1
                    else:#elif int(sorted_list[j - 1][0]) != int((sorted_list[j][0])) and int(sorted_list[j + 1][0]) != int((sorted_list[j][0])):
                        if (sorted_list[j][5] == "DIP"  and '-' not in sorted_list[j][3])or sorted_list[j][5] == "Insertion" and int(sorted_list[j + 1][0]) != int((sorted_list[j][0])): 
                            if sorted_list[j][4] == '-':
                                DNA_strand = "-"
                                end_position = GetHomopolymerLengthIns(GetReverseComplement(new_list2[j][5]), int(new_list2[j][7]), int(new_list2[j][0]), int(len(new_list2[j][5])/2) + 1, DNA_strand)[1][1] 
                            else:
                                DNA_strand = "+"
                                end_position = GetHomopolymerLengthIns(new_list2[j][5], int(new_list2[j][7]), int(new_list2[j][0]), int(len(new_list2[j][5])/2) + 1, DNA_strand)[1][1] 
                        else:    
                            if sorted_list[j][4] == '-':
                                end_position = GetHomopolymerLength(GetReverseComplement(new_list2[j][5]), int(new_list2[j][7]), int(new_list2[j][0]), int(len(new_list2[j][5])/2) + 1)[1][1] 
                            else:
                                end_position = GetHomopolymerLength(new_list2[j][5], int(new_list2[j][7]), int(new_list2[j][0]), int(len(new_list2[j][5])/2) + 1)[1][1] 
                        if int(sorted_list[j + 1][0]) <= (end_position + 10) and sorted_list[j][1] == sorted_list[j + 1][1]:
                            if (sorted_list[j][5] == "DIP"or sorted_list[j][5] == "Deletion" or sorted_list[j][5] == "Insertion"):
                            #if (sorted_list[j][5] == "SNV"or sorted_list[j][5] == "SNP") and j in new_list2.keys():
                                    #if sorted_list[j][4] == '-':
                                if (sorted_list[j][6], sorted_list[j][1], sorted_list[j][0]) not in clusters.keys():
                                    clusters[(sorted_list[j][6], sorted_list[j][1], sorted_list[j][0])] = [sorted_list[j]]
                            if (sorted_list[j + 1][5] == "DIP"  or sorted_list[j + 1][5] == "Deletion" or sorted_list[j + 1][5] == "Insertion"): 
                            #if (sorted_list[j + 1][5] == "SNV"  or sorted_list[j + 1][5] == "SNP") and j + 1 in new_list2.keys():
                                if (sorted_list[j + 1][6], sorted_list[j + 1][1], sorted_list[j + 1][0]) not in clusters.keys():
                                    clusters[(sorted_list[j + 1][6], sorted_list[j + 1][1], sorted_list[j + 1][0])] = [sorted_list[j + 1]]
                        if (sorted_list[j][5] == "DIP"  and '-' not in sorted_list[j][3])or sorted_list[j][5] == "Insertion" and int(sorted_list[j - 1][0]) != int((sorted_list[j][0])):
                            if sorted_list[j][4] == '-':
                                DNA_strand = "-"
                                start_position = GetHomopolymerLengthIns(GetReverseComplement(new_list2[j][5]), int(new_list2[j][7]), int(new_list2[j][0]), int(len(new_list2[j][5])/2) + 1, DNA_strand)[1][0] 
                            else:
                                DNA_strand = "+"
                                start_position = GetHomopolymerLengthIns(new_list2[j][5], int(new_list2[j][7]), int(new_list2[j][0]), int(len(new_list2[j][5])/2) + 1, DNA_strand)[1][0] 
                        else:    
                            if sorted_list[j][4] == '-':
                                start_position = GetHomopolymerLength(GetReverseComplement(new_list2[j][5]), int(new_list2[j][7]), int(new_list2[j][0]), int(len(new_list2[j][5])/2) + 1)[1][0] 
                            else:
                                start_position = GetHomopolymerLength(new_list2[j][5], int(new_list2[j][7]), int(new_list2[j][0]), int(len(new_list2[j][5])/2) + 1)[1][0] 
                        if start_position <= int((sorted_list[j - 1][0] + 10)) and sorted_list[j - 1][1] == sorted_list[j][1]:
                            if (sorted_list[j][5] == "DIP" or sorted_list[j][5] == "Deletion" or sorted_list[j][5] == "Insertion"): 
                            #if (sorted_list[j][5] == "SNV" or sorted_list[j][5] == "SNP") and j in new_list2.keys(): 
                                if (sorted_list[j][6], sorted_list[j][1], sorted_list[j][0]) not in clusters.keys():
                                    clusters[(sorted_list[j][6], sorted_list[j][1], sorted_list[j][0])] = [sorted_list[j]]  #indel_clusters = indel_clusters + 1
                            if (sorted_list[j - 1][5] == "DIP" or sorted_list[j - 1][5] == "Deletion" or sorted_list[j - 1][5] == "Insertion"): 
                            #if (sorted_list[j - 1][5] == "SNV" or sorted_list[j - 1][5] == "SNP") and j - 1 in new_list2.keys():    
                                if (sorted_list[j - 1][6], sorted_list[j - 1][1], sorted_list[j - 1][0]) not in clusters.keys():
                                    clusters[(sorted_list[j - 1][6], sorted_list[j - 1][1], sorted_list[j - 1][0])] = [sorted_list[j - 1]]
                    
                        
                        
                    #elif isolate_list[key][j] not in nonclusters.keys() and (isolate_list[key][j][7] == "DIP" or isolate_list[key][0][7] == "Complex DIP" or isolate_list[key][0][7] == "Deletion" or isolate_list[key][0][7] == "Insertion"):
                                #if j == 0 or isolate_list[key][j - n] not in clusters.keys():
                                #nonclusters[(isolate_list[key][j][6], isolate_list[key][j][1], isolate_list[key][j][0])] = [isolate_list[key][j]]
            for k in range(0, len(sorted_list)):
                if (sorted_list[k][5] == "DIP" or sorted_list[k][5] == "Deletion" or sorted_list[k][5] == "Insertion"):
                #if (sorted_list[k][5] == "SNV" or sorted_list[k][5] == "SNP"):
                    if (sorted_list[k][6], sorted_list[k][1], sorted_list[k][0]) not in clusters.keys() and (sorted_list[k][6], sorted_list[k][1], sorted_list[k][0]) not in nonclusters.keys():
                        nonclusters[(sorted_list[k][6], sorted_list[k][1], sorted_list[k][0])] = [sorted_list[k]]


    homopolymer_deletions = {}
    homopolymer_insertions = {}
    homopolymer_subs = {}
    complex_deletions = []
    complex_insertions = []
    complex_subs = []
    for key in clusters.keys():
        if len(clusters[key][0][8]) == 1:
            if homopolymers[key] not in homopolymer_subs.keys():
                homopolymer_subs[homopolymers[key]] = 1
            else:
                homopolymer_subs[homopolymers[key]] = homopolymer_subs[homopolymers[key]] + 1
        if '-' in clusters[key][0][3] and len(clusters[key][0][8]) == 1:
            if homopolymers[key] not in homopolymer_deletions.keys():
                homopolymer_deletions[homopolymers[key]] = 1
            else:
                homopolymer_deletions[homopolymers[key]] = homopolymer_deletions[homopolymers[key]] + 1
        if '-' not in clusters[key][0][3] and  len(clusters[key][0][3]) == 1:
            if homopolymers[key] not in homopolymer_insertions.keys():
                homopolymer_insertions[homopolymers[key]] = 1
            else:
                homopolymer_insertions[homopolymers[key]] = homopolymer_insertions[homopolymers[key]] + 1
    for key in homopolymers.keys():
        if int(homopolymers[key]) >= 4 and key in clusters.keys():

            complex_subs.append(key)
        if int(homopolymers[key]) >= 4 and key in clusters.keys() and '-' in clusters[key][0][3]:

            complex_deletions.append(key)

        if int(homopolymers[key]) >= 4 and key in clusters.keys() and '-' not in clusters[key][0][3]:
            complex_insertions.append(key)
    single_insertion = []
    multi_insertion = []
    single_deletion = []
    multi_deletion = []
    for key in clusters.keys():
        if clusters[key][0][3] == '-' and len(clusters[key][0][8]) == 1:
            single_deletion.append(clusters[key])
        elif '-' in clusters[key][0][3]:
            multi_deletion.append(clusters[key])

        elif len(clusters[key][0][3]) == 1:
            single_insertion.append(clusters[key])
        else:
            multi_insertion.append(clusters[key])
    return clusters, nonclusters


clusters = FindClusters(isolate_list, rad30_isolatelist, homopolymer3)
clusters2 = FindClusters(WT_isolate_list, WT_isolatelist2 ,homopolymer6)


def index2(tuple):
    return tuple[2]

chromosome_sequences = {}
chromosome_sequences['chrI'] = sequence1
chromosome_sequences['chrII'] = sequence2
chromosome_sequences['chrIII'] = sequence3
chromosome_sequences['chrIV'] = sequence4
chromosome_sequences['chrV'] = sequence5
chromosome_sequences['chrVI'] = sequence6
chromosome_sequences['chrVII'] = sequence7
chromosome_sequences['chrVIII'] = sequence8
chromosome_sequences['chrIX'] = sequence9
chromosome_sequences['chrX'] = sequence10
chromosome_sequences['chrXI'] = sequence11
chromosome_sequences['chrXII'] = sequence12
chromosome_sequences['chrXIII'] = sequence13
chromosome_sequences['chrXIV'] = sequence14
chromosome_sequences['chrXV'] = sequence15
chromosome_sequences['chrXVI'] = sequence16
    

d1 = open('rad30_sub_clusters.bed', 'w+')
d2 = open('rad30_sub_sequences', 'w+')
d3 = open('rad30_sub_nonclusters.bed', 'w+')
d4 = open('rad30_sub_sequences2', 'w+')





def MakeClusterFiles(clusters, nonclusters, file1, file2, file3, file4):
    cluster_keys = sorted(clusters.keys(), key = index2)
    sorted_list = sorted(cluster_keys, key = index)
    file2.write("Deletions: ")
    file2.write('\n')
    for c in sorted_list:
        if '-' in clusters[c][0][3] and clusters[c][0][7] == 1:
            file1.write(clusters[c][0][1] + '\t' + str((clusters[c][0][0] - 1)) + '\t' + str(clusters[c][0][0]) + '\t' + clusters[c][0][2] + '\t' + clusters[c][0][3] + '\t' + clusters[c][0][4])
            file1.write('\n')
        for key in chromosome_sequences.keys():
            if clusters[c][0][1] == key:
                sequence = chromosome_sequences[key]
                seq = sequence[int(clusters[c][0][0]) - 11: int(clusters[c][0][0]) + 10]
                if clusters[c][0][0] < 11:
                        seq = sequence[0: int(clusters[c][0][0]) + 10]
                        while len(seq) < 21:
                            seq = 'N' + seq
                else:
                        while len(seq) < 21:
                            seq = seq + 'N'

                #file2.write('\n')
                if '-' in clusters[c][0][3]:
                    if clusters[c][0][4] == '-' and  clusters[c][0][7] == 1:
                        file2.write(GetReverseComplement(seq[::-1]))
                    elif clusters[c][0][7] == 1:
                        file2.write(seq)
                    file2.write('\n')
    file2.write('Insertions: ')
    file2.write('\n')
    for c in sorted_list:
        for key in chromosome_sequences.keys():
            if clusters[c][0][1] == key:
                sequence = chromosome_sequences[key]
                seq = sequence[int(clusters[c][0][0]) - 12: int(clusters[c][0][0]) + 9]
                if clusters[c][0][0] < 10:
                    seq = sequence[0: int(clusters[c][0][0]) + 10]
                    while len(seq) < 21:
                        seq = 'N' + seq
                else:
                     while len(seq) < 21:
                            seq = seq + 'N'
                #file2.write('\n')
                if '-' not in clusters[c][0][3]:

                    if clusters[c][0][4] == '-' and len(clusters[c][0][3]) == 1:
                        file2.write(GetReverseComplement(seq[::-1]))
                    elif len(clusters[c][0][3]) == 1:
                        file2.write(seq)
                    file2.write('\n')
    noncluster_keys = sorted(nonclusters.keys(), key = index2)
    sorted_list2 = sorted(noncluster_keys, key = index)
    file4.write('Deletions: ')
    file4.write('\n')
    for d in sorted_list2:
        if '-' in nonclusters[d][0][3] and nonclusters[d][0][7] == 1:
            file3.write(nonclusters[d][0][1] + '\t' + str((nonclusters[d][0][0] - 1)) + '\t' + str(nonclusters[d][0][0]) + '\t' + nonclusters[d][0][2] + '\t' + nonclusters[d][0][3] + '\t' + nonclusters[d][0][4])
            file3.write('\n')
        for key in chromosome_sequences.keys():
            if nonclusters[d][0][1] == key:
                sequence = chromosome_sequences[key]
                seq = sequence[int(nonclusters[d][0][0]) - 11: int(nonclusters[d][0][0]) + 10]
                if nonclusters[d][0][0] < 11:
                        seq = sequence[0: int(nonclusters[d][0][0]) + 10]
                        while len(seq) < 21:
                            seq = 'N' + seq
                else:
                        while len(seq) < 21:
                            seq = seq + 'N'


                #file4.write('\n')
                if '-' in nonclusters[d][0][3]:

                    if nonclusters[d][0][4] == '-' and nonclusters[d][0][7] == 1:
                        file4.write(GetReverseComplement(seq[::-1]))
                    elif nonclusters[d][0][7] == 1:
                        file4.write(seq)
                    file4.write('\n')
    file4.write('Insertions: ')
    file4.write('\n')
    for d in sorted_list2:
        for key in chromosome_sequences.keys():
            if nonclusters[d][0][1] == key:
                sequence = chromosome_sequences[key]
                seq = sequence[int(nonclusters[d][0][0]) - 12: int(nonclusters[d][0][0]) + 9]
                if nonclusters[d][0][0] < 10:
                        seq = sequence[0: int(nonclusters[d][0][0]) + 10]
                        while len(seq) < 21:
                            seq = 'N' + seq
                else:
                        while len(seq) < 21:
                            seq = seq + 'N'

                #file4.write('\n')
                if '-' not in nonclusters[d][0][3]:

                    if nonclusters[d][0][4] == '-' and len(nonclusters[d][0][3]) == 1:
                        file4.write(GetReverseComplement(seq[::-1]))
                    elif len(nonclusters[d][0][3]) == 1:
                        file4.write(seq)
                    file4.write('\n')

                

MakeClusterFiles(clusters[0], clusters[1], d1, d2, d3, d4)


d1.close()
d2.close()
d3.close()
d4.close()



h1 = open("homopolymer_deletions.txt", 'w+')
h2 = open("homopolymer_insertions.txt", 'w+')
h3 = open("nonhomopolymer_deletions.txt", 'w+')
h4 = open("nonhomopolymer_insertions.txt", 'w+')
h5 = open("homopolymer_indels.txt", 'w+')
h6 = open("nonhomopolymer_indels.txt", 'w+')
h13 = open("all_indels.txt", 'w+')
h15 = open("all_deletions.txt", 'w+')
h16 = open("all_insertions.txt", 'w+')

h7 = open("WT_homopolymer_deletions.txt", 'w+')
h8 = open("WT_homopolymer_insertions.txt", 'w+')
h9 = open("WT_nonhomopolymer_deletions.txt", 'w+')
h10 = open("WT_nonhomopolymer_insertions.txt", 'w+')
h11 = open("WT_homopolymer_indels.txt", 'w+')
h12 = open("WT_nonhomopolymer_indels.txt", 'w+')
h14 = open("WT_all_indels.txt", 'w+')
h17 = open("WT_all_deletions.txt", 'w+')
h18 = open("WT_all_insertions.txt", 'w+')


def HomopolymerMutationsPerIsolate( isolate_list, h1, h2, h3, h4, h5, h6, h13, h15, h16):
    HP_Deletions = []
    HP_Insertions = []
    HP_rad30_deletion = 0
    HP_rad30_insertion = 0
    nonHP_Deletions = []
    nonHP_Insertions = []
    nonHP_rad30_deletion = 0
    nonHP_rad30_insertion = 0
    HP_indels = []
    nonHP_indels = []
    HP_count = 0
    nonHP_count = 0

    for key in isolate_list.keys():
        #if key[0] in isolate_keys:

        d_count = 0
        i_count = 0
        d_count2 = 0
        i_count2 = 0
        h_count = 0
        h_count2 = 0
        for a in range(0, len(isolate_list[key])):
            if isolate_list[key][a][9] == 'DIP'or isolate_list[key][a][9] == 'Deletion' or isolate_list[key][a][9] == 'Insertion':
            #if isolate_list[key][a][9] == 'Complex DIP' or isolate_list[key][a][9] == 'DIP'or isolate_list[key][a][9] == 'Deletion' or isolate_list[key][a][9] == 'Insertion':
                if int(isolate_list[key][a][8]) >=4:
                    h_count = h_count + 1

                    if '-' in isolate_list[key][a][3]:# and  len(isolate_list[key][a][7]) == 1:
                        d_count = d_count + 1
                    elif '-' not in isolate_list[key][a][3]:# and len(isolate_list[key][a][3]) == 1:
                        i_count = i_count + 1
                else:
                    if int(isolate_list[key][a][8]) < 4:
                        h_count2 = h_count2 + 1

                        if '-' in isolate_list[key][a][3]:# and  len(isolate_list[key][a][7]) == 1:
                            d_count2 = d_count2 + 1
                        elif '-' not in isolate_list[key][a][3]:# and len(isolate_list[key][a][3]) == 1:
                            i_count2 = i_count2 + 1

      
        HP_Deletions.append(d_count)
        h1.write(str(d_count))
        h1.write("\n")
        HP_Insertions.append(i_count)
        h2.write(str(i_count))
        h2.write("\n")
        
        nonHP_Deletions.append(d_count2)
        h3.write(str(d_count2))
        h3.write("\n")
        nonHP_Insertions.append(i_count2)
        h4.write(str(i_count2))
        h4.write("\n")

        
        HP_indels.append(h_count)
        h5.write(str(h_count))
        h5.write("\n")
        nonHP_indels.append(h_count2)
        h6.write(str(h_count2))
        h6.write("\n")

        h13.write(str(h_count + h_count2))
        h13.write("\n")

        h15.write(str(d_count + d_count2))
        h15.write("\n")

        h16.write(str(i_count + i_count2))
        h16.write("\n")
    
    for e in HP_Deletions:
        HP_rad30_deletion = HP_rad30_deletion+ e
    HP_rad30_deletion = HP_rad30_deletion/len(HP_Deletions)
    print(HP_Deletions)

    for i in HP_Insertions:
        HP_rad30_insertion = HP_rad30_insertion + i
    HP_rad30_insertion = HP_rad30_insertion/len(HP_Insertions)
    print(HP_Insertions)

    for f in nonHP_Deletions:
        nonHP_rad30_deletion = nonHP_rad30_deletion+ f
    nonHP_rad30_deletion = nonHP_rad30_deletion/len(nonHP_Deletions)
    print(nonHP_Deletions)

    for j in nonHP_Insertions:
        nonHP_rad30_insertion = nonHP_rad30_insertion + j
    nonHP_rad30_insertion = nonHP_rad30_insertion/len(nonHP_Insertions)
    print(nonHP_Insertions)

    for g in HP_indels:
        HP_count = HP_count + g
    HP_count = HP_count/len(HP_indels)
    print(HP_indels)

    for k in nonHP_indels:
        nonHP_count = nonHP_count + k
    nonHP_count = nonHP_count/len(nonHP_indels)
    print(nonHP_indels)

    
    return HP_Deletions, HP_Insertions, nonHP_Deletions, nonHP_Insertions, HP_indels, nonHP_indels
       

HomopolymerMutationsPerIsolate(rad30_isolatelist, h1, h2, h3, h4, h5, h6, h13, h15, h16)
HomopolymerMutationsPerIsolate(WT_isolatelist2, h7, h8, h9, h10, h11, h12, h14, h17, h18)

h1.close()
h2.close()
h3.close()
h4.close()
h5.close()
h6.close()
h7.close()
h8.close()
h9.close()
h10.close()
h11.close()
h12.close()
h13.close()
h14.close()
h15.close()
h16.close()
h17.close()
h18.close()

x1 = open("complex_deletions.txt", 'w+')
x2 = open("complex_insertions.txt", 'w+')
x3 = open("noncomplex_deletions.txt", 'w+')
x4 = open("noncomplex_insertions.txt", 'w+')
x5 = open("complex_indels.txt", 'w+')
x6 = open("noncomplex_indels.txt", 'w+')

x7 = open("WT_complex_deletions.txt", 'w+')
x8 = open("WT_complex_insertions.txt", 'w+')
x9 = open("WT_noncomplex_deletions.txt", 'w+')
x10 = open("WT_noncomplex_insertions.txt", 'w+')
x11 = open("WT_complex_indels.txt", 'w+')
x12 = open("WT_noncomplex_indels.txt", 'w+')


def ComplexMutationsPerIsolate( isolate_list, clusters, x1, x2, x3, x4, x5, x6):
    Complex_Deletions = []
    Complex_Insertions = []
    Complex_rad30_deletion = 0
    Complex_rad30_insertion = 0
    nonComplex_Deletions = []
    nonComplex_Insertions = []
    nonComplex_rad30_deletion = 0
    nonComplex_rad30_insertion = 0
    Complex_indels = []
    nonComplex_indels = []
    complex_count = 0
    noncomplex_count = 0

    for key in isolate_list.keys():
        #if key[0] in isolate_keys:

        d_count = 0
        i_count = 0
        d_count2 = 0
        i_count2 = 0
        c_count = 0
        c_count2 = 0
        for a in range(0, len(isolate_list[key])):
            if isolate_list[key][a][9] == 'DIP'or isolate_list[key][a][9] == 'Deletion' or isolate_list[key][a][9] == 'Insertion':
            #if isolate_list[key][a][9] == 'Complex DIP' or isolate_list[key][a][9] == 'DIP'or isolate_list[key][a][9] == 'Deletion' or isolate_list[key][a][9] == 'Insertion':
                if (isolate_list[key][a][6], isolate_list[key][a][1], isolate_list[key][a][0]) in clusters.keys():
                    c_count = c_count + 1

                    if '-' in isolate_list[key][a][3]:# and  len(isolate_list[key][a][7]) == 1:
                        d_count = d_count + 1
                    elif '-' not in isolate_list[key][a][3]:# and len(isolate_list[key][a][3]) == 1:
                        i_count = i_count + 1
                else:
                    if (isolate_list[key][a][6], isolate_list[key][a][1], isolate_list[key][a][0]) not in clusters.keys():
                        c_count2 = c_count2 + 1

                        if '-' in isolate_list[key][a][3]:# and  len(isolate_list[key][a][7]) == 1:
                            d_count2 = d_count2 + 1
                        elif '-' not in isolate_list[key][a][3]:# and len(isolate_list[key][a][3]) == 1:
                            i_count2 = i_count2 + 1

      
        Complex_Deletions.append(d_count)
        x1.write(str(d_count))
        x1.write("\n")
        Complex_Insertions.append(i_count)
        x2.write(str(i_count))
        x2.write("\n")
        
        nonComplex_Deletions.append(d_count2)
        x3.write(str(d_count2))
        x3.write("\n")
        nonComplex_Insertions.append(i_count2)
        x4.write(str(i_count2))
        x4.write("\n")

        
        Complex_indels.append(c_count)
        x5.write(str(c_count))       
        x5.write("\n")
        nonComplex_indels.append(c_count2)
        x6.write(str(c_count2))
        x6.write("\n")
    
    for e in Complex_Deletions:
        Complex_rad30_deletion = Complex_rad30_deletion+ e
    Complex_rad30_deletion = Complex_rad30_deletion/len(Complex_Deletions)
    print(Complex_Deletions)

    for i in Complex_Insertions:
        Complex_rad30_insertion = Complex_rad30_insertion + i
    Complex_rad30_insertion = Complex_rad30_insertion/len(Complex_Insertions)
    print(Complex_Insertions)

    for f in nonComplex_Deletions:
        nonComplex_rad30_deletion = nonComplex_rad30_deletion+ f
    nonComplex_rad30_deletion = nonComplex_rad30_deletion/len(nonComplex_Deletions)
    print(nonComplex_Deletions)

    for j in nonComplex_Insertions:
        nonComplex_rad30_insertion = nonComplex_rad30_insertion + j
    nonComplex_rad30_insertion = nonComplex_rad30_insertion/len(nonComplex_Insertions)
    print(nonComplex_Insertions)

    for g in Complex_indels:
        complex_count = complex_count + g
    complex_count = complex_count/len(Complex_indels)
    print(Complex_indels)

    for k in nonComplex_indels:
        noncomplex_count = noncomplex_count + k
    noncomplex_count = noncomplex_count/len(nonComplex_indels)
    print(nonComplex_indels)
    
    return Complex_Deletions, Complex_Insertions, nonComplex_Deletions, nonComplex_Insertions, Complex_indels, nonComplex_indels
       

ComplexMutationsPerIsolate(rad30_isolatelist, clusters[0], x1, x2, x3, x4, x5, x6)
ComplexMutationsPerIsolate(WT_isolatelist2, clusters2[0], x7, x8, x9, x10, x11, x12)

x1.close()
x2.close()
x3.close()
x4.close()
x5.close()
x6.close()
x7.close()
x8.close()
x9.close()
x10.close()
x11.close()
x12.close()

t1 = open("C_deletions_nonclusters.txt", 'w+')
t2 = open("c_insertions_nonclusters.txt", 'w+')
t3 = open("t_deletions_nonclusters.txt", 'w+')
t4 = open("t_insertions_nonclusters.txt", 'w+')
t5 = open("c_indels_nonclusters.txt", 'w+')
t6 = open("t_indels_nonclusters.txt", 'w+')

t7 = open("WT_c_deletions_nonclusters.txt", 'w+')
t8 = open("WT_c_insertions_nonclusters.txt", 'w+')
t9 = open("WT_t_deletions_nonclusters.txt", 'w+')
t10 = open("WT_t_insertions_nonclusters.txt", 'w+')
t11 = open("WT_c_indels_nonclusters.txt", 'w+')
t12 = open("WT_t_indels_nonclusters.txt", 'w+')



def CTMutationsPerIsolate( isolate_list, t1, t2, t3, t4, t5, t6, clusters):
    C_Deletions = []
    C_Insertions = []
    C_rad30_deletion = 0
    C_rad30_insertion = 0
    T_Deletions = []
    T_Insertions = []
    T_rad30_deletion = 0
    T_rad30_insertion = 0
    C_indels = []
    T_indels = []
    C_count = 0
    T_count = 0

    for key in isolate_list.keys():
        #if key[0] in isolate_keys:

        d_count = 0
        i_count = 0
        d_count2 = 0
        i_count2 = 0
        h_count = 0
        h_count2 = 0
        for a in range(0, len(isolate_list[key])):
            if isolate_list[key][a][9] == 'DIP'or isolate_list[key][a][9] == 'Deletion' or isolate_list[key][a][9] == 'Insertion':
            #if isolate_list[key][a][9] == 'Complex DIP' or isolate_list[key][a][9] == 'DIP'or isolate_list[key][a][9] == 'Deletion' or isolate_list[key][a][9] == 'Insertion':
                if (isolate_list[key][a][6], isolate_list[key][a][1], isolate_list[key][a][0]) in clusters.keys():
                    if isolate_list[key][a][10] == 'C':
                        d_count = d_count + 1
                        h_count = h_count + 1
                    if isolate_list[key][a][3] == 'C':
                        i_count = i_count + 1
                        h_count = h_count + 1
                    if isolate_list[key][a][10] == 'T':
                        d_count2 = d_count2 + 1
                        h_count2 = h_count2 + 1
                    if isolate_list[key][a][3] == 'T':
                        i_count2 = i_count2 + 1
                        h_count2 = h_count2 + 1

                   

      
        C_Deletions.append(d_count)
        t1.write(str(d_count))
        t1.write("\n")
        C_Insertions.append(i_count)
        t2.write(str(i_count))
        t2.write("\n")
        
        T_Deletions.append(d_count2)
        t3.write(str(d_count2))
        t3.write("\n")
        T_Insertions.append(i_count2)
        t4.write(str(i_count2))
        t4.write("\n")

        
        C_indels.append(h_count)
        t5.write(str(h_count))
        t5.write("\n")
        T_indels.append(h_count2)
        t6.write(str(h_count2))
        t6.write("\n")

    
    return C_Deletions, C_Insertions, T_Deletions, T_Insertions, C_indels, T_indels


CTMutationsPerIsolate(rad30_isolatelist, t1, t2, t3, t4, t5, t6, clusters[1])
CTMutationsPerIsolate(WT_isolatelist2, t7, t8, t9, t10, t11, t12, clusters2[1])

t1.close()
t2.close()
t3.close()
t4.close()
t5.close()
t6.close()
t7.close()
t8.close()
t9.close()
t10.close()
t11.close()
t12.close()






