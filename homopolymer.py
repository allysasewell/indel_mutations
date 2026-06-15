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
    str2 = ''
    for i in range(0, len(str)):
        str2 = str2 + dict[str[i]]
    return str2[::-1]

def GetComplement(str):
    if '-' in str:
        return str
    if 'N' in str:
        return 'N'
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


#Return reference chromosome sequence from the SacCer3 genome
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
        


        


#with open("UV_yeast_Merged_SNP_DIPs_sorted_anz2_tandem.txt") as f:
with open("UVC_filter_mutations.txt") as f:
    data = f.read()
lines = [s.strip().split() for s in data.splitlines()]
genotypes = {'WT': 0, 'rad16':0, 'rad26':0, 'rad30': 0}
WT_chromosomes = {}
rad16_chromosomes = {}
rad26_chromosomes = {}
rad30_chromosomes = {}
WT_indels = 0
rad16_indels = 0
rad26_indels = 0
rad30_indels = 0
genotype = []
background = []
chromosome = []
position1 = []
position2 = []
mutation_type = []
length = []
allele = []
mutation = []
frequency = []
left_reference = []
right_reference = []
left_consensus = []
right_consensus = []
trinucleotide = []
strand = []
isolate = []
isolate_list = {}
lines.pop(0)
for line in lines:
   
    if 'Mito' not in line[6] and 'ySR128' not in line[3] and line[0] != 'dRad16_25':#and line[10] == "DIP" or line[10] == "Complex DIP":
        isolate.append(line[0])
    
        genotype.append(line[4])
    
        mutation_type.append(line[10])
        if 'WT' in line[4]:# and line[10] == "SNP" or line[10] == 'SNV':
            genotypes['WT'] = genotypes['WT'] + 1
        if 'rad16' in line[4]:# and  line[10] == "SNP" or line[10] == 'SNV':
            genotypes['rad16'] = genotypes['rad16'] + 1
        if 'rad26' in line[4]:# and  line[10] == "SNP" or line[10] == 'SNV':
            genotypes['rad26'] = genotypes['rad26'] + 1
        if 'rad30' in line[4]:# and  line[10] == "SNP" or line[10] == 'SNV':
            genotypes['rad30'] = genotypes['rad30'] + 1
        chromosome.append('chr' + str(line[6]))
        if 'WT' in line[4]:
            if str(line[6]) in WT_chromosomes.keys():
                WT_chromosomes[str(line[6])] = WT_chromosomes[str(line[6])] + 1

            else:
                WT_chromosomes[str(line[6])] = 1
        if 'rad16' in line[4]:
            if str(line[6]) in rad16_chromosomes.keys():
                rad16_chromosomes[str(line[6])] = rad16_chromosomes[str(line[6])] + 1
            else:
                rad16_chromosomes[str(line[6])] = 1
        if 'rad26' in line[4]:
            if str(line[6]) in rad26_chromosomes.keys():
                rad26_chromosomes[str(line[6])] = rad26_chromosomes[str(line[6])] + 1
            else:
                rad26_chromosomes[str(line[6])] = 1
        if 'rad30' in line[4]:
            if str(line[6]) in rad30_chromosomes.keys():
                rad30_chromosomes[str(line[6])] = rad30_chromosomes[str(line[6])] + 1
            else:
                rad30_chromosomes[str(line[6])] = 1
        position1.append(int(line[9]) - len(line[13]))
       
        position2.append(int(line[9]))
        length.append(line[11])
        allele.append(line[12])
        if len(line) >= 22: 
            if '-' not in  line[18] and line[18].isalpha() == False:
                    if float(line[18]) >= 90:
                        frequency.append('Homozygous')
                    else:
                        frequency.append('Heterozygous')
            elif float(line[19]) >= 90:
                    frequency.append('Homozygous')
            else:
                    frequency.append('Heterozygous')
            left_reference.append(line[23])
            right_reference.append(line[24])
            left_consensus.append(line[25])
            right_consensus.append(line[26])
            if line[27].isalpha() == False and len(line[24]) <= 10:
                    if line[24][1] == "G" or line[24][1] == "A":
                        #trinucleotide.append(GetReverseComplement(line[24]))
                        #strand.append('-')
                        if '-' not in line[13]:
                            mutation.append(GetReverseComplement(line[13]))
                        else:
                            mutation.append(line[13])
                    else:
                        #trinucleotide.append(line[24])
                        mutation.append(line[13])
                        #strand.append('+')
            elif len(line[27]) <= 10:
                    if line[27][1] == "A" or line[27][1] == "G":
                        #trinucleotide.append(GetReverseComplement(line[27]))
                        #strand.append('-')
                        if '-' not in line[13]:
                            mutation.append(GetReverseComplement(line[13]))
                        else:
                            mutation.append(line[13])
                    else:
                        #trinucleotide.append(line[27])
                        mutation.append(line[13])
                        #strand.append('+')
            else:
                    if line[28][1] == "A" or line[28][1] == "G":
                        #trinucleotide.append(GetReverseComplement(line[28]))
                        #strand.append('-')
                        if '-' not in line[13]:
                            mutation.append(GetReverseComplement(line[13]))
                        else:
                            mutation.append(line[13])
                    else:
                        #trinucleotide.append(line[28])
                        mutation.append(line[13])
                        #strand.append('+')
        


    
        else:
            frequency.append('NA')
            left_reference.append('NA')
            right_reference.append('NA')
            left_consensus.append('NA')
            right_consensus.append('NA')
            #trinucleotide.append('NA')
            mutation.append(line[13])
       


for pos in range(0, len(position2)):
    #if mutation_type_2[pos] == 'Insertion' or mutation_type_2[pos] == 'Deletion':# or mutation_type_2[pos] == 'Deletion' and len(allele2[pos]) == 1:
        if chromosome[pos] == 'chrI':
            sequence = GetSequence('chr1.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrII':
            sequence = GetSequence('chr2.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrIII':
            sequence = GetSequence('chr3.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrIV':
            sequence = GetSequence('chr4.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrV':
            sequence = GetSequence('chr5.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrVI':
            sequence = GetSequence('chr6.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrVII':
            sequence = GetSequence('chr7.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrVIII':
            sequence = GetSequence('chr8.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrIX':
            sequence = GetSequence('chr9.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrX':
            sequence = GetSequence('chr10.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrXI':
            sequence = GetSequence('chr11.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrXII':
            sequence = GetSequence('chr12.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrXIII':
            sequence = GetSequence('chr13.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrXIV':
            sequence = GetSequence('chr14.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrXV':
            sequence = GetSequence('chr15.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
                strand.append('-')
        elif chromosome[pos] == 'chrXVI':
            sequence = GetSequence('chr16.txt')
            if sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'C' or sequence[int(position1[pos]) - 1: int(position2[pos]) + 1][1] == 'T':
                trinucleotide.append(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1])
                strand.append('+')
            else:
                trinucleotide.append(GetReverseComplement(sequence[int(position1[pos]) - 1: int(position2[pos]) + 1]))
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

    if (isolate[i1], genotype[i1]) not in isolate_list.keys():
            isolate_list[(isolate[i1], genotype[i1])] = [(position2[i1], 'chr' + str(chromosome[i1]), trinucleotide[i1], mutation[i1], strand[i1], mutation_type[i1], isolate[i1])]
    else:
        isolate_list[(isolate[i1], genotype[i1])] .append((position2[i1], 'chr' + str(chromosome[i1]), trinucleotide[i1], mutation[i1], strand[i1], mutation_type[i1], isolate[i1]))
    






def getStrforInt(string):
    new_str = ""
    for s in string:
        if s.isnumeric():
            new_str = new_str + s
        else:
            return new_str
    return new_str

with open("UVB_filter_mutations.txt") as f2:
    data2 = f2.read()
lines2 = [s.strip().split() for s in data2.splitlines()]
lines2.pop(0)
genotypes2 = {'WT': 0, 'Rad16':0, 'Rad26':0, 'Rad30': 0}
WT_indels_2 = 0
rad16_indels_2 = 0
rad26_indels_2 = 0
rad30_indels_2 = 0
isolates2 = []
genotype2 = []
chromosome2 = []
position1_2 = []
position2_2 = []
mutation_type_2 = []
length2 = []
allele2 = []
mutation2 = []
frequency2 = []
trinucleotide2 = []
strand2 = []
isolate_list2 = {}
WT_chromosomes2 = {}
rad16_chromosomes2 = {}
rad26_chromosomes2 = {}
rad30_chromosomes2 = {}
for line in lines2:
    if 'Mito' not in line[5] and line[2] != 'RP' and line[0] != 'B01' and line[0] != 'C13':
        isolates2.append(line[0])
        genotype2.append(line[1])
        mutation_type_2.append(line[13])
        if 'WT' in line[1]:
            genotypes2['WT'] = genotypes2['WT'] + 1
        if 'Rad16' in line[1]:
            genotypes2['Rad16'] = genotypes2['Rad16'] + 1
        if 'Rad26' in line[1]:
            genotypes2['Rad26'] = genotypes2['Rad26'] + 1
        if 'Rad30' in line[1]:
            genotypes2['Rad30'] = genotypes2['Rad30'] + 1
        if line[2] != 'RP':
            chromosome2.append(line[4])
        else:
            chromosome2.append(line[5])
        if 'WT' in line[1]:
            if chromosome2[len(chromosome2) - 1] in WT_chromosomes2.keys():
                WT_chromosomes2[chromosome2[len(chromosome2) - 1]] = WT_chromosomes2[chromosome2[len(chromosome2) - 1]] + 1

            else:
                WT_chromosomes2[chromosome2[len(chromosome2) - 1]] = 1
            if line[13] == "Deletion" or line[13] == "Insertion" or line[14] == "Deletion" or line[14] == "Insertion":
                WT_indels_2 = WT_indels_2 + 1
        if 'Rad16' in line[1]:
            if chromosome2[len(chromosome2) - 1] in rad16_chromosomes2.keys():
                rad16_chromosomes2[chromosome2[len(chromosome2) - 1]] = rad16_chromosomes2[chromosome2[len(chromosome2) - 1]] + 1
            else:
                rad16_chromosomes2[chromosome2[len(chromosome2) - 1]] = 1
            if line[13] == "Deletion" or line[13] == "Insertion" or line[14] == "Deletion" or line[14] == "Insertion":
                rad16_indels_2 = rad16_indels_2 + 1
        if 'Rad26' in line[1]:
            if chromosome2[len(chromosome2) - 1] in rad26_chromosomes2.keys():
                rad26_chromosomes2[chromosome2[len(chromosome2) - 1]] = rad26_chromosomes2[chromosome2[len(chromosome2) - 1]] + 1
            else:
                rad26_chromosomes2[chromosome2[len(chromosome2) - 1]] = 1
            if line[13] == "Deletion" or line[13] == "Insertion" or line[14] == "Deletion" or line[14] == "Insertion":
                rad26_indels_2 = rad26_indels_2 + 1
        if 'Rad30' in line[1]:
            if chromosome2[len(chromosome2) - 1] in rad30_chromosomes2.keys():
                rad30_chromosomes2[chromosome2[len(chromosome2) - 1]] = rad30_chromosomes2[chromosome2[len(chromosome2) - 1]] + 1
            else:
                rad30_chromosomes2[chromosome2[len(chromosome2) - 1]] = 1
            if line[13] == "Deletion" or line[13] == "Insertion" or line[14] == "Deletion" or line[14] == "Insertion":
                rad30_indels_2 = rad30_indels_2 + 1
        if line[6].isnumeric():
            if line[14] == '-':
                position1_2.append(int(line[6]) - len(line[15]))
            else:
                position1_2.append(int(line[6]) - len(line[14]))
            position2_2.append(int(line[6]))
        else:
           
            new_str = getStrforInt(line[6])

            position2_2.append(int(new_str))
            position1_2.append(int(new_str) - len(line[14]))
        allele2.append(line[14])
        #if line[14] == 'C' or line[14] == 'T':
        mutation2.append(line[15])
         


for pos in range(0, len(position2_2)):
    #if mutation_type_2[pos] == 'Insertion' or mutation_type_2[pos] == 'Deletion':# or mutation_type_2[pos] == 'Deletion' and len(allele2[pos]) == 1:
        if chromosome2[pos] == 'chrI':
            sequence = GetSequence('chr1.txt')
            if sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'C' or sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'T':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
                strand2.append('+')
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
                strand2.append('-')

        elif chromosome2[pos] == 'chrII':
            sequence = GetSequence('chr2.txt')
            if sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'C' or sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'T':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
                strand2.append('+')
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
                strand2.append('-')
        elif chromosome2[pos] == 'chrIII':
            sequence = GetSequence('chr3.txt')
            if sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'C' or sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'T':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
                strand2.append('+')
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
                strand2.append('-')
        elif chromosome2[pos] == 'chrIV':
            sequence = GetSequence('chr4.txt')
            if sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'C' or sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'T':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
                strand2.append('+')
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
                strand2.append('-')
        elif chromosome2[pos] == 'chrV':
            sequence = GetSequence('chr5.txt')
            if sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'C' or sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'T':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
                strand2.append('+')
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
                strand2.append('-')
        elif chromosome2[pos] == 'chrVI':
            sequence = GetSequence('chr6.txt')
            if sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'C' or sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'T':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
                strand2.append('+')
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
                strand2.append('-')
        elif chromosome2[pos] == 'chrVII':
            sequence = GetSequence('chr7.txt')
            if sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'C' or sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'T':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
                strand2.append('+')
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
                strand2.append('-')
        elif chromosome2[pos] == 'chrVIII':
            sequence = GetSequence('chr8.txt')
            if sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'C' or sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'T':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
                strand2.append('+')
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
                strand2.append('-')
        elif chromosome2[pos] == 'chrIX':
            sequence = GetSequence('chr9.txt')
            if sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'C' or sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'T':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
                strand2.append('+')
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
                strand2.append('-')
        elif chromosome2[pos] == 'chrX':
            sequence = GetSequence('chr10.txt')
            if sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'C' or sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'T':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
                strand2.append('+')
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
                strand2.append('-')
        elif chromosome2[pos] == 'chrXI':
            sequence = GetSequence('chr11.txt')
            if sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'C' or sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'T':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
                strand2.append('+')
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
                strand2.append('-')
        elif chromosome2[pos] == 'chrXII':
            sequence = GetSequence('chr12.txt')
            if sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'C' or sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'T':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
                strand2.append('+')
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
                strand2.append('-')
        elif chromosome2[pos] == 'chrXIII':
            sequence = GetSequence('chr13.txt')
            if sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'C' or sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'T':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
                strand2.append('+')
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
                strand2.append('-')
        elif chromosome2[pos] == 'chrXIV':
            sequence = GetSequence('chr14.txt')
            if sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'C' or sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'T':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
                strand2.append('+')
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
                strand2.append('-')
        elif chromosome2[pos] == 'chrXV':
            sequence = GetSequence('chr15.txt')
            if sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'C' or sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'T':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
                strand2.append('+')
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
                strand2.append('-')
        elif chromosome2[pos] == 'chrXVI':
            sequence = GetSequence('chr16.txt')
            if sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'C' or sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1][1] == 'T':
                trinucleotide2.append(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1])
                strand2.append('+')
            else:
                trinucleotide2.append(GetReverseComplement(sequence[int(position1_2[pos]) - 1: int(position2_2[pos]) + 1]))
                strand2.append('-')
        else:
            trinucleotide2.append('NA')
            strand2.append('+')
    #else:
        #trinucleotide2.append('NA')
for m1 in  range(0, len(mutation2)):
    if mutation_type_2[m1] == 'Insertion':
        if mutation2[m1] == 'A' or mutation2[m1] == 'G':
            strand2[m1] = '-'
        else:
            strand2[m1] = '+'
    elif strand2[m1] == '-':
        mutation2[m1] = GetReverseComplement(mutation2[m1])
for i1 in range(0, len(isolates2)):

    if (isolates2[i1], genotype2[i1]) not in isolate_list2.keys():
            isolate_list2[(isolates2[i1], genotype2[i1])] = [(position2_2[i1], 'chr' + str(chromosome2[i1]), trinucleotide2[i1], mutation2[i1], strand2[i1], mutation_type_2[i1], isolates2[i1])]
    else:
        isolate_list2[(isolates2[i1], genotype2[i1])] .append((position2_2[i1], 'chr' + str(chromosome2[i1]), trinucleotide2[i1], mutation2[i1], strand2[i1], mutation_type_2[i1], isolates2[i1]))

#f1a = open("WT_substitutions_sorted_tandem.bed", 'w+')
#f1b = open("rad16_substitutions_sorted_tandem.bed", 'w+')
#f1c = open("rad26_substitutions_sorted_tandem.bed", 'w+')
#f1d = open("rad30_substitutions_sorted_tandem.bed", 'w+')
f1a  = open("WT_mutations_sorted_tandem.bed", 'w+')
f1b = open("rad16_mutations_sorted_tandem.bed", 'w+')
f1c = open("rad26_mutations_sorted_tandem.bed", 'w+')
f1d = open("rad30_mutations_sorted_tandem.bed", 'w+')

f1ua  = open("WT_UVB_mutations_sorted_tandem.bed", 'w+')
f1ub = open("rad16_UVB_mutations_sorted_tandem.bed", 'w+')
f1uc = open("rad26_UVB_mutations_sorted_tandem.bed", 'w+')
f1ud = open("rad30_UVB_mutations_sorted_tandem.bed", 'w+')

#Makes bed files 
def WriteChr(chr_name, file1a, file1b, file1c, file1d, chromosome, position1, position2,
             trinucleotide, allele, mutation, genotype, mutation_type ,isolate, strand):
    file_data = []
    chromosome1_dict = {}
    chromosome1 =[]
    for x in range(0, len(mutation_type)):
        #if mutation_type[x] == 'SNV' or mutation_type[x] == 'SNP':
                if trinucleotide[x].isalpha():# len(trinucleotide[x]) == 3 and trinucleotide[x].isalpha():
                    if chromosome[x] == chr_name: 
                        chromosome1_dict[int(position2[x]), isolate[x]] = x
                


    for key1 in chromosome1_dict.keys():
        chromosome1.append((int(key1[0]), key1[1]))
    sorted_chromosome1= sorted(chromosome1)
    #sorted_chromosome1.sort()
    for y in range(0, len(sorted_chromosome1)):
        file_data.append(chromosome[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + str(int(position1[chromosome1_dict[sorted_chromosome1[y]]])) + '\t' + str(position2[chromosome1_dict[sorted_chromosome1[y]]]) + '\t' + trinucleotide[chromosome1_dict[sorted_chromosome1[y]]]  + '\t' + mutation[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + strand[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + allele[chromosome1_dict[sorted_chromosome1[y]]] + '\t' )
        #print(chromosome[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + str(int(position1[chromosome1_dict[sorted_chromosome1[y]]]) - 1) + '\t' + position2[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + trinucleotide[chromosome1_dict[sorted_chromosome1[y]]]  + '\t' + mutation[chromosome1_dict[sorted_chromosome1[y]]])
    for fd in range(0, len(file_data)):
        
            if genotype[chromosome1_dict[sorted_chromosome1[fd]]] == 'WT':
                file1a.write(file_data[fd])
               
                file1a.write(isolate[chromosome1_dict[sorted_chromosome1[fd]]] + '\t' + mutation_type[chromosome1_dict[sorted_chromosome1[fd]]])
                file1a.write('\n')
                
            if 'rad16' in genotype[chromosome1_dict[sorted_chromosome1[fd]]] or 'Rad16' in  genotype[chromosome1_dict[sorted_chromosome1[fd]]]:
                file1b.write(file_data[fd])
                
                file1b.write(isolate[chromosome1_dict[sorted_chromosome1[fd]]] + '\t' + mutation_type[chromosome1_dict[sorted_chromosome1[fd]]])
                file1b.write('\n')
              
            if 'rad26' in genotype[chromosome1_dict[sorted_chromosome1[fd]]] or 'Rad26' in genotype[chromosome1_dict[sorted_chromosome1[fd]]]:
    
                file1c.write(file_data[fd])
               
                file1c.write(isolate[chromosome1_dict[sorted_chromosome1[fd]]] + '\t' + mutation_type[chromosome1_dict[sorted_chromosome1[fd]]])
                file1c.write('\n')
              
            if 'rad30' in genotype[chromosome1_dict[sorted_chromosome1[fd]]] or 'Rad30' in  genotype[chromosome1_dict[sorted_chromosome1[fd]]]:
                file1d.write(file_data[fd])
              
                file1d.write(isolate[chromosome1_dict[sorted_chromosome1[fd]]] + '\t' + mutation_type[chromosome1_dict[sorted_chromosome1[fd]]])
                file1d.write('\n')
               

       

WriteChr('chrI',f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrII', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrIII',  f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrIV', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrIX', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrV', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrVI', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrVII', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrVIII', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrX', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXI', f1a, f1b, f1c, f1d,  chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXII', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXIII', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXIV', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand) 
WriteChr('chrXV', f1a, f1b, f1c, f1d,  chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXVI', f1a, f1b, f1c, f1d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)


WriteChr('chrI',f1ua, f1ub, f1uc, f1ud, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrII', f1ua, f1ub, f1uc, f1ud, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrIII',  f1ua, f1ub, f1uc, f1ud, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrIV', f1ua, f1ub, f1uc, f1ud, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrIX', f1ua, f1ub, f1uc, f1ud, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrV', f1ua, f1ub, f1uc, f1ud, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrVI', f1ua, f1ub, f1uc, f1ud, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrVII', f1ua, f1ub, f1uc, f1ud, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrVIII', f1ua, f1ub, f1uc, f1ud, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrX', f1ua, f1ub, f1uc, f1ud, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrXI', f1ua, f1ub, f1uc, f1ud, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrXII',f1ua, f1ub, f1uc, f1ud, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrXIII', f1ua, f1ub, f1uc, f1ud, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrXIV', f1ua, f1ub, f1uc, f1ud, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2) 
WriteChr('chrXV', f1ua, f1ub, f1uc, f1ud, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrXVI',f1ua, f1ub, f1uc, f1ud, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)


f1a.close()
f1b.close()
f1c.close()
f1d.close()

f1ua.close()
f1ub.close()
f1uc.close()
f1ud.close()

#f11a.close()
#f11b.close()
#f11c.close()
#f11d.close()

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

    
    
    return multinucleotide, new_strand, new_chromosome, new_position, new_trinucleotide, new_mutation, new_allele, lengths, new_position1, new_isolate, new_mutation_type

#Counts the number of mutations in homopolymers of different lengths
def GetPolymerCounts(pentanucleotide,  new_chromosome, new_position, new_trinucleotide, new_mutation, new_allele, lengths, new_isolate, new_mutation_type):
    homopolymers = {}
    homopolymer_list = {}
    polymer_counts = []
    homopolymer_number = 0
    for j in range(0, 4):
        polymer_counts.append({})
    for i in range (0, len(lengths)):
        if len(new_mutation[i]) == 1:
            print(pentanucleotide[i])
            print(new_trinucleotide[i])
            key = lengths[i]
            
            #elif strand[i] == '-':
                #key = GetHomopolymerLength(GetReversComplement(sequence), )[1]
            
            if key >= 4 and (new_mutation_type[i] == 'SNP' or new_mutation_type[i] == 'SNV'):
                homopolymer_number = homopolymer_number + 1
                homopolymers[(new_isolate[i], new_chromosome[i], new_position[i])] = key
                #homopolymer_list[key[1]] = [(new_chromosome[i], key[0][0], key[0][1])]
            print(new_trinucleotide[i])
            print(pentanucleotide[i])
            print(key)

            if (new_allele[i] == 'C' or new_allele[i] == 'G') and (new_mutation_type[i] == 'SNP' or new_mutation_type[i] == 'SNV'):

                if key in polymer_counts[0].keys():
                    polymer_counts[0][key] = polymer_counts[0][key] + 1
                else:
                    polymer_counts[0][key] = 1
            if (new_allele[i] == 'A' or new_allele[i] == 'T') and (new_mutation_type[i] == 'SNP' or new_mutation_type[i] == 'SNV'):

                if key in polymer_counts[1].keys():
                    polymer_counts[1][key] = polymer_counts[1][key] + 1
                else:
                    polymer_counts[1][key] = 1
            
    return polymer_counts, homopolymers, homopolymer_list

#f2 = open("WT_UV_bothruns_Muts_allSNVs_sorted.bed")
f2 = open("WT_mutations_sorted_tandem.bed")
f2a = open("WT_mutations_sorted_multinucleotide.bed", 'w+')
#f2a = open("WT_all_sorted_multinucleotide.bed", 'w+')
#f3 = open("rad16_UV_bothruns_Muts_allSNVs_sorted.bed")
f3 = open("rad16_mutations_sorted_tandem.bed")
f3a = open("rad16_mutations_multinucleotide.bed", 'w+')
#f3a = open("rad16_all_sorted_multinucleotide.bed", 'w+')
#f4 = open("rad26_UV_bothruns_Muts_allSNVs_sorted.bed")
f4 = open("rad26_mutations_sorted_tandem.bed")
f4a = open("rad26_mutations_multinucleotide.bed", 'w+')
#f4a = open("rad26_all_multinucleotide.bed", 'w+')
#f10 = open("rad30_UV_bothruns_Muts_allSNVs_sorted.bed")
f10a = open("rad30_substitutions_multinucleotide.bed", 'w+')
f10 = open("rad30_mutations_sorted_tandem.bed")
#f10a = open("rad30_all_multinucleotide.bed", 'w+')
f2b = open('WT_multinucleotide_subs', 'w+')
f3b = open('rad16_multinucleotide_subs', 'w+')
f4b = open('rad26_multinucleotide_subs', 'w+')
f10b = open('rad30_multinucleotide_subs', 'w+')
#f11 = open("WT_UVB_Muts_allSNVs_sorted.bed")
f11 = open("WT_UVB_mutations_sorted_tandem.bed")
#f12 = open("rad16_UVB_Muts_allSNVs_sorted.bed")
f12 = open("rad16_UVB_mutations_sorted_tandem.bed")
#f13 = open("rad26_UVB_Muts_allSNVs_sorted.bed")
f13 = open("rad26_UVB_mutations_sorted_tandem.bed")
#f14 = open("rad30_UVB_Muts_allSNVs_sorted.bed")
f14 = open("rad30_UVB_mutations_sorted_tandem.bed")
#f11a = open("WT_UVB_substitutions_sorted_multinucleotide.bed", 'w+')
#f12a = open("rad16_UVB_substitutions_sorted_multinucleotide.bed", 'w+')
#f13a = open("rad26_UVB_substitutions_sorted_multinucleotide.bed", 'w+')
#f14a = open("rad30_UVB_substitutions_sorted_multinucleotide.bed", 'w+')
f11a = open("WT_UVB_all_sorted_multinucleotide.bed", 'w+')
f12a = open("rad16_UVB_all_sorted_multinucleotide.bed", 'w+')
f13a = open("rad26_UVB_all_sorted_multinucleotide.bed", 'w+')
f14a = open("rad30_UVB_all_sorted_multinucleotide.bed", 'w+')
f11b = open('WT_UVB_multinucleotide_subs', 'w+')
f12b = open('rad16_UVB_multinucleotide_subs', 'w+')
f13b = open('rad26_UVB_multinucleotide_subs', 'w+')
f14b = open('rad30_UVB_multinucleotide_subs', 'w+')

WT_data = MakeMultiNucFile(f2, f2a)
WT_isolate = WT_data[9]
WT_chromosome = WT_data[2]
WT_position1 = WT_data[8]
WT_position = WT_data[3]
WT_trinucleotide = WT_data[4]
WT_mutation = WT_data[5]
WT_strand = WT_data[1]
WT_multinucleotide = WT_data[0] 
WT_lengths = WT_data[7]
WT_mutation_type = WT_data[10]


rad16_data = MakeMultiNucFile(f3, f3a)
rad16_isolate = rad16_data[9]
rad16_chromosome = rad16_data[2]
rad16_position1 = rad16_data[8]
rad16_position = rad16_data[3]
rad16_trinucleotide = rad16_data[4]
rad16_mutation = rad16_data[5]
rad16_strand = rad16_data[1]
rad16_multinucleotide = rad16_data[0] 
rad16_lengths = rad16_data[7]
rad16_mutation_type = rad16_data[10]

rad26_data = MakeMultiNucFile(f4, f4a)

rad26_isolate = rad26_data[9]
rad26_chromosome = rad26_data[2]
rad26_position1 = rad26_data[8]
rad26_position = rad26_data[3]
rad26_trinucleotide = rad26_data[4]
rad26_mutation = rad26_data[5]
rad26_strand = rad26_data[1]
rad26_multinucleotide = rad26_data[0] 
rad26_lengths = rad26_data[7]

rad30_data = MakeMultiNucFile(f10, f10a)
rad30_isolate = rad30_data[9]
rad30_chromosome = rad30_data[2]
rad30_position1 = rad30_data[8]
rad30_position = rad30_data[3]
rad30_trinucleotide = rad30_data[4]
rad30_mutation = rad30_data[5]
rad30_strand = rad30_data[1]
rad30_multinucleotide = rad30_data[0] 
rad30_lengths = rad30_data[7]
rad26_mutation_type = rad26_data[10]


WT_dataB = MakeMultiNucFile(f11, f11a)
WT_isolateB = WT_dataB[9]
WT_chromosomeB = WT_dataB[2]
WT_position1B = WT_dataB[8]
WT_positionB = WT_dataB[3]
WT_trinucleotideB = WT_dataB[4]
WT_mutationB = WT_dataB[5]
WT_strandB = WT_dataB[1]
WT_multinucleotideB = WT_dataB[0] 
WT_lengthsB = WT_dataB[7]
WT_mutation_typeB = WT_dataB[10]


rad16_dataB = MakeMultiNucFile(f12, f12a)
rad16_isolateB = rad16_dataB[9]
rad16_chromosomeB = rad16_dataB[2]
rad16_position1B = rad16_dataB[8]
rad16_positionB = rad16_dataB[3]
rad16_trinucleotideB = rad16_dataB[4]
rad16_mutationB = rad16_dataB[5]
rad16_strandB = rad16_dataB[1]
rad16_multinucleotideB = rad16_dataB[0] 
rad16_lengthsB = rad16_dataB[7]
rad16_mutation_typeB = rad16_dataB[10]

rad26_dataB = MakeMultiNucFile(f13, f13a)
rad26_isolateB = rad26_dataB[9]

rad26_chromosomeB = rad26_dataB[2]
rad26_position1B = rad26_dataB[8]
rad26_positionB = rad26_dataB[3]
rad26_trinucleotideB = rad26_dataB[4]
rad26_mutationB = rad26_dataB[5]
rad26_strandB = rad26_dataB[1]
rad26_multinucleotideB = rad26_dataB[0] 
rad26_lengthsB = rad26_dataB[7]
rad26_mutation_typeB = rad26_dataB[10]

rad30_dataB = MakeMultiNucFile(f14, f14a)
rad30_isolateB = rad30_dataB[9]
rad30_chromosomeB = rad30_dataB[2]
rad30_position1B = rad30_dataB[8]
rad30_positionB = rad30_dataB[3]
rad30_trinucleotideB = rad30_dataB[4]
rad30_mutationB = rad30_dataB[5]
rad30_strandB = rad30_dataB[1]
rad30_multinucleotideB = rad30_dataB[0] 
rad30_lengthsB = rad30_dataB[7]



WT_polymer_counts = GetPolymerCounts(WT_data[0],  WT_data[2], WT_data[3],  WT_data[4], WT_data[5], WT_data[6], WT_data[7], WT_data[9], WT_data[10])[0]
WT_Ccount = WT_polymer_counts[0]
WT_Tcount = WT_polymer_counts[1]
rad16_polymer_counts = GetPolymerCounts(rad16_data[0],  rad16_data[2], rad16_data[3], rad16_data[4] , rad16_data[5], rad16_data[6], rad16_data[7], rad16_data[9], rad16_data[10])[0]
rad16_Ccount = rad16_polymer_counts[0]
rad16_Tcount = rad16_polymer_counts[1]
WT_homopolymers = GetPolymerCounts(WT_data[0],  WT_data[2], WT_data[3],  WT_data[4], WT_data[5], WT_data[6], WT_data[7], WT_data[9], WT_data[10])[1]
rad16_homopolymers = GetPolymerCounts(rad16_data[0],  rad16_data[2], rad16_data[3],  rad16_data[4], rad16_data[5], rad16_data[6], rad16_data[7],  rad16_data[9], rad16_data[10])[1]
rad26_homopolymers = GetPolymerCounts(rad26_data[0],  rad26_data[2], rad26_data[3],  rad26_data[4], rad26_data[5], rad26_data[6], rad26_data[7],  rad26_data[9], rad26_data[10])[1]
rad26_polymer_counts = GetPolymerCounts(rad26_data[0],  rad26_data[2], rad26_data[3], rad26_data[4], rad26_data[5], rad26_data[6], rad26_data[7],  rad26_data[9], rad26_data[10])[0]
rad26_Ccount = rad26_polymer_counts[0]
rad26_Tcount = rad26_polymer_counts[1]
rad30_polymer_counts = GetPolymerCounts(rad30_data[0],  rad30_data[2], rad30_data[3], rad30_data[4], rad30_data[5], rad30_data[6], rad30_data[7], rad30_data[9], rad30_data[10])[0]
rad30_Ccount = rad30_polymer_counts[0]
rad30_Tcount = rad30_polymer_counts[1]

WT_polymer_countsB = GetPolymerCounts(WT_dataB[0],  WT_dataB[2], WT_dataB[3],  WT_dataB[4], WT_dataB[5], WT_dataB[6], WT_dataB[7], WT_dataB[9], WT_dataB[10])[0]
WT_CcountB = WT_polymer_countsB[0]
WT_TcountB = WT_polymer_countsB[1]
rad16_polymer_countsB = GetPolymerCounts(rad16_dataB[0],  rad16_dataB[2], rad16_dataB[3], rad16_dataB[4] , rad16_dataB[5], rad16_dataB[6], rad16_dataB[7], rad16_dataB[9], rad16_dataB[10])[0]
rad16_CcountB = rad16_polymer_countsB[0]
rad16_TcountB = rad16_polymer_countsB[1]
WT_homopolymersB = GetPolymerCounts(WT_dataB[0],  WT_dataB[2], WT_dataB[3],  WT_dataB[4], WT_dataB[5], WT_dataB[6], WT_dataB[7], WT_dataB[9], WT_dataB[10])[1]
rad16_homopolymersB = GetPolymerCounts(rad16_dataB[0],  rad16_dataB[2], rad16_dataB[3],  rad16_dataB[4], rad16_dataB[5], rad16_dataB[6], rad16_dataB[7], rad16_dataB[9], rad16_dataB[10])[1]
rad26_homopolymersB = GetPolymerCounts(rad26_dataB[0],  rad26_dataB[2], rad26_dataB[3],  rad26_dataB[4], rad26_dataB[5], rad26_dataB[6], rad26_dataB[7], rad26_dataB[9], rad26_dataB[10])[1]
rad26_polymer_countsB = GetPolymerCounts(rad26_dataB[0],  rad26_dataB[2], rad26_dataB[3], rad26_dataB[4], rad26_dataB[5], rad26_dataB[6], rad26_dataB[7], rad26_dataB[9], rad26_dataB[10])[0]
rad26_CcountB = rad26_polymer_countsB[0]
rad26_TcountB = rad26_polymer_countsB[1]
rad30_polymer_countsB = GetPolymerCounts(rad30_dataB[0],  rad30_dataB[2], rad30_dataB[3], rad30_dataB[4], rad30_dataB[5], rad30_dataB[6], rad30_dataB[7], rad30_dataB[9], rad30_dataB[10])[0]
rad30_CcountB = rad30_polymer_countsB[0]
rad30_TcountB = rad30_polymer_countsB[1]

f2.close()
f2a.close()
f2b.close()
f3.close()
f3a.close()
f3b.close()
f4.close()
f4a.close()
f4b.close()
f10.close()
f10a.close()
f10b.close()



f11.close()
f11a.close()
f11b.close()
f12.close()
f12a.close()
f12b.close()
f13.close()
f13a.close()
f13b.close()
f14.close()
f14a.close()
f14b.close()


WT_isolate_list = {}
for i1 in range(0, len(WT_isolate)):

    if WT_isolate[i1] not in WT_isolate_list.keys():
            WT_isolate_list[WT_isolate[i1]] = [(WT_position[i1],  str(WT_chromosome[i1]), WT_trinucleotide[i1], WT_mutation[i1], WT_strand[i1], WT_multinucleotide[i1], WT_isolate[i1], WT_position1[i1], WT_mutation_type[i1])]
    else:
        WT_isolate_list[WT_isolate[i1]] .append((WT_position[i1],  str(WT_chromosome[i1]), WT_trinucleotide[i1], WT_mutation[i1], WT_strand[i1], WT_multinucleotide[i1], WT_isolate[i1], WT_position1[i1], WT_mutation_type[i1]))


rad16_isolate_list = {}
for i1 in range(0, len(rad16_isolate)):

    if rad16_isolate[i1] not in rad16_isolate_list.keys():
            rad16_isolate_list[rad16_isolate[i1]] = [(rad16_position[i1],  str(rad16_chromosome[i1]), rad16_trinucleotide[i1], rad16_mutation[i1], rad16_strand[i1], rad16_multinucleotide[i1], rad16_isolate[i1], rad16_position1[i1], rad16_mutation_type[i1])]
    else:
        rad16_isolate_list[rad16_isolate[i1]] .append((rad16_position[i1],  str(rad16_chromosome[i1]), rad16_trinucleotide[i1], rad16_mutation[i1], rad16_strand[i1], rad16_multinucleotide[i1], rad16_isolate[i1], rad16_position1[i1], rad16_mutation_type[i1]))


rad26_isolate_list = {}
for i1 in range(0, len(rad26_isolate)):

    if rad26_isolate[i1] not in rad26_isolate_list.keys():
            rad26_isolate_list[rad26_isolate[i1]] = [(rad26_position[i1],  str(rad26_chromosome[i1]), rad26_trinucleotide[i1], rad26_mutation[i1], rad26_strand[i1], rad26_multinucleotide[i1], rad26_isolate[i1], rad26_position1[i1], rad26_mutation_type[i1])]
    else:
        rad26_isolate_list[rad26_isolate[i1]] .append((rad26_position[i1],  str(rad26_chromosome[i1]), rad26_trinucleotide[i1], rad26_mutation[i1], rad26_strand[i1], rad26_multinucleotide[i1], rad26_isolate[i1], rad26_position1[i1], rad26_mutation_type[i1]))



WT_isolate_listB = {}
for i1 in range(0, len(WT_isolateB)):

    if WT_isolateB[i1] not in WT_isolate_listB.keys():
            WT_isolate_listB[WT_isolateB[i1]] = [(WT_positionB[i1], str(WT_chromosomeB[i1]), WT_trinucleotideB[i1], WT_mutationB[i1], WT_strandB[i1], WT_multinucleotideB[i1], WT_isolateB[i1], WT_position1B[i1],WT_mutation_typeB[i1] )]
    else:
        WT_isolate_listB[WT_isolateB[i1]].append((WT_positionB[i1],  str(WT_chromosomeB[i1]), WT_trinucleotideB[i1], WT_mutationB[i1], WT_strandB[i1], WT_multinucleotideB[i1], WT_isolateB[i1], WT_position1B[i1], WT_mutation_typeB[i1]))

rad16_isolate_listB = {}
for i1 in range(0, len(rad16_isolateB)):

    if rad16_isolateB[i1] not in rad16_isolate_listB.keys():
            rad16_isolate_listB[rad16_isolateB[i1]] = [(rad16_positionB[i1],  str(rad16_chromosomeB[i1]), rad16_trinucleotideB[i1], rad16_mutationB[i1], rad16_strandB[i1], rad16_multinucleotideB[i1], rad16_isolateB[i1], rad16_position1B[i1], rad16_mutation_typeB[i1])]
    else:
        rad16_isolate_listB[rad16_isolateB[i1]] .append((rad16_positionB[i1],  str(rad16_chromosomeB[i1]), rad16_trinucleotideB[i1], rad16_mutationB[i1], rad16_strandB[i1], rad16_multinucleotideB[i1], rad16_isolateB[i1], rad16_position1B[i1], rad16_mutation_typeB[i1]))

rad26_isolate_listB = {}
for i1 in range(0, len(rad26_isolateB)):

    if rad26_isolateB[i1] not in rad26_isolate_listB.keys():
            rad26_isolate_listB[rad26_isolateB[i1]] = [(rad26_positionB[i1],  str(rad26_chromosomeB[i1]), rad26_trinucleotideB[i1], rad26_mutationB[i1], rad26_strandB[i1], rad26_multinucleotideB[i1], rad26_isolateB[i1], rad26_position1B[i1], rad26_mutation_typeB[i1])]
    else:
        rad26_isolate_listB[rad26_isolateB[i1]] .append((rad26_positionB[i1],  str(rad26_chromosomeB[i1]), rad26_trinucleotideB[i1], rad26_mutationB[i1], rad26_strandB[i1], rad26_multinucleotideB[i1], rad26_isolateB[i1], rad26_position1B[i1], rad26_mutation_typeB[i1]))

#Obtain the total number of bases in the yeast genome in sequences of each homopolymer length
def NormalizeHomopolymer(sequences,  file):
    homopolymer_dict = {}
    C_dict = {}
    T_dict = {}
    
    homopolymer_dict['11-15'] = 0
    homopolymer_dict['16-20'] = 0
    homopolymer_dict['21+'] = 0
    C_dict['11-15'] = 0
    C_dict['16-20'] = 0
    C_dict['21+'] = 0
    T_dict['11-15'] = 0
    T_dict['16-20'] = 0
    T_dict['21+'] = 0
  
    homopolymer_dict[1] = 0
    C_dict[1] = 0
    T_dict[1] = 0
    for seq in sequences:
        #coordinates = []
        sequence = seq
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
                                   #T_dict[polymer_count] = T_dict[polymer_count] + 1
             
                b = b + length
    for key in homopolymer_dict.keys():
        file.write(str(key) + '\t' + str(homopolymer_dict[key]))
        file.write('\n')
    
    return homopolymer_dict, C_dict, T_dict
f5 = open('Homopolymer_numbers', 'w+')

new_homopolymer_dict = NormalizeHomopolymer(sequences, f5 )
f5.close()

f6 = open('WT_CHomopolymer_counts', 'w+')
f7 = open('rad16_CHomopolymer_counts', 'w+')
f8 = open('rad26_CHomopolymer_counts', 'w+')
f9 = open('rad30_CHomopolymer_counts', 'w+')
f6a = open('WT_THomopolymer_counts', 'w+')
f7a = open('rad16_THomopolymer_counts', 'w+')
f8a = open('rad26_THomopolymer_counts', 'w+')
f9a = open('rad30_THomopolymer_counts', 'w+')

f16 = open('WT_UVB_Homopolymer_counts', 'w+')
f17 = open('rad16_UVB_Homopolymer_counts', 'w+')
f18 = open('rad26_UVB_Homopolymer_counts', 'w+')
f19 = open('rad30_UVB_Homopolymer_counts', 'w+')
f16a = open('WT_UVB_CHomopolymer_counts', 'w+')
f17a = open('rad16_UVB_CHomopolymer_counts', 'w+')
f18a = open('rad26_UVB_CHomopolymer_counts', 'w+')
f19a = open('rad30_UVB_CHomopolymer_counts', 'w+')
f16b = open('WT_UVB_THomopolymer_counts', 'w+')
f17b = open('rad16_UVB_THomopolymer_counts', 'w+')
f18b = open('rad26_UVB_THomopolymer_counts', 'w+')
f19b = open('rad30_UVB_THomopolymer_counts', 'w+')

def PrintFrequencies(homopolymer_dict, lengths, file):
    homopolymer_counts = {}
    homopolymer_counts['21+'] = 0
    homopolymer_counts['16-20'] = 0
    homopolymer_counts['11-15'] = 0
    for l in range(0, len(lengths)):
        if lengths[l] != 'NA':
            #if len(allele[l]) == 1: #and mutation_type[l] == "deletion") or (len(mutation[l]) == 1 and mutation_type[l] == "insertion") :
                if lengths[l] not in homopolymer_counts.keys():
                    homopolymer_counts[lengths[l]] = 1
                    #homopolymer_dict[lengths[l]] = 0
                else:
                    homopolymer_counts[lengths[l]] = homopolymer_counts[lengths[l]] + 1
                if int(lengths[l]) >= 11:
                    if int(lengths[l]) >= 21:
                            homopolymer_counts['21+'] = homopolymer_counts['21+'] + 1
                    elif int(lengths[l]) >= 16:
                            homopolymer_counts['16-20'] = homopolymer_counts['16-20'] + 1
                    else:
                            homopolymer_counts['11-15'] = homopolymer_counts['11-15'] + 1
                
    for key in homopolymer_counts.keys():
        if homopolymer_dict[key] != 0:
            homopolymer_counts[key] = homopolymer_counts[key]/homopolymer_dict[key]
            file.write(str(key) + ':' + str(homopolymer_counts[key]))
            file.write('\n')
   
#PrintFrequencies(homopolymer_dict, WT_lengths, f6 )
#PrintFrequencies(homopolymer_dict, rad16_lengths, f7 )
#PrintFrequencies(homopolymer_dict, rad26_lengths, f8 )
#PrintFrequencies(homopolymer_dict, rad30_lengths, f9 )
#f6.close()
#f7.close()
#f8.close()
#f9.close()

def PrintCTFrequencies(homopolymer_dict, counts, file):
    homopolymer_counts = counts
    homopolymer_counts['21+'] = 0
    homopolymer_counts['16-20'] = 0
    homopolymer_counts['11-15'] = 0
    for key in counts.keys():
        if key != '11-15' and key != '16-20' and key != '21+':
            #if len(allele[l]) == 1: #and mutation_type[l] == "deletion") or (len(mutation[l]) == 1 and mutation_type[l] == "insertion") :
                #if key not in homopolymer_counts.keys():
                    #homopolymer_counts[key] = 1
                    #homopolymer_dict[lengths[l]] = 0
                #else:
                    #homopolymer_counts[key] = homopolymer_counts[key] + 1
                if int(key) >= 11:
                    if int(key) >= 21:
                            homopolymer_counts['21+'] = homopolymer_counts['21+'] + homopolymer_counts[key]
                    elif int(key) >= 16:
                            homopolymer_counts['16-20'] = homopolymer_counts['16-20'] + homopolymer_counts[key]
                    else:
                            homopolymer_counts['11-15'] = homopolymer_counts['11-15'] + homopolymer_counts[key]
                
    for key in homopolymer_counts.keys():
        if key in homopolymer_dict.keys():
            if homopolymer_dict[key] != 0:
                homopolymer_counts[key] = homopolymer_counts[key]/homopolymer_dict[key]
                file.write(str(key) + ':' + str(homopolymer_counts[key]))
                file.write('\n')

PrintFrequencies(new_homopolymer_dict[0], WT_lengthsB, f16 )
PrintFrequencies(new_homopolymer_dict[0], rad16_lengthsB, f17 )
PrintFrequencies(new_homopolymer_dict[0], rad26_lengthsB, f18 )
PrintFrequencies(new_homopolymer_dict[0], rad30_lengthsB, f19 )
f16.close()
f17.close()
f18.close()
f19.close()

PrintCTFrequencies(new_homopolymer_dict[1], WT_Ccount, f6)
PrintCTFrequencies(new_homopolymer_dict[1], rad16_Ccount, f7)
PrintCTFrequencies(new_homopolymer_dict[1], rad26_Ccount, f8)
PrintCTFrequencies(new_homopolymer_dict[1], rad30_Ccount, f9)
PrintCTFrequencies(new_homopolymer_dict[2], WT_Tcount, f6a)
PrintCTFrequencies(new_homopolymer_dict[2], rad16_Tcount, f7a)
PrintCTFrequencies(new_homopolymer_dict[2], rad26_Tcount, f8a)
PrintCTFrequencies(new_homopolymer_dict[2], rad30_Tcount, f9a)
PrintCTFrequencies(new_homopolymer_dict[1], WT_CcountB, f16a)
PrintCTFrequencies(new_homopolymer_dict[1], rad16_CcountB, f17a)
PrintCTFrequencies(new_homopolymer_dict[1], rad26_CcountB, f18a)
PrintCTFrequencies(new_homopolymer_dict[1], rad30_CcountB, f19a)
PrintCTFrequencies(new_homopolymer_dict[2], WT_TcountB, f16b)
PrintCTFrequencies(new_homopolymer_dict[2], rad16_TcountB, f17b)
PrintCTFrequencies(new_homopolymer_dict[2], rad26_TcountB, f18b)
PrintCTFrequencies(new_homopolymer_dict[2], rad30_TcountB, f19b)

f6.close()
f6a.close()
f7.close()
f7a.close()
f8.close()
f8a.close()
f9.close()
f9a.close()

f16a.close()
f16b.close()
f17a.close()
f17b.close()
f18a.close()
f18b.close()
f19a.close()
f19b.close()


def index(tuple):
    return tuple[1]   

def FindClusters( isolate_list, homopolymers):
    lengths = {}
    clusters = {}
    nonclusters = {}
    total_clusters = 0
    for key in isolate_list.keys():
            sorted_isolate_list = sorted(isolate_list[key])
            sorted_list = sorted(sorted_isolate_list, key = index)
            if sorted_list[0][1] == sorted_list[1][1] and (int(sorted_list[1][0]) <= int((sorted_list[ 0][0])) + 10):
                if (sorted_list[0][6], sorted_list[0][1], sorted_list[0][6]) not in clusters.keys() and (sorted_list[0][8] == 'SNV' or sorted_list[0][8] == 'SNP'):
                    clusters[(sorted_list[0][6], sorted_list[0][1], sorted_list[0][6])] = [sorted_list[0]]
                   
                if (sorted_list[1][6], sorted_list[1][1], sorted_list[1][6]) not in clusters.keys() and (sorted_list[1][8] == 'SNV' or sorted_list[1][8] == 'SNP'):
                    clusters[(sorted_list[1][6], sorted_list[1][1], sorted_list[1][6])] = [sorted_list[1]]
                    
            else:#if strand[0] == '+':
                if (sorted_list[0][8] == "DIP"  and '-' not in sorted_list[0][3])or sorted_list[0][8] == "Insertion": 
                            if sorted_list[0][4] == '-':
                                DNA_strand = "-"
                                end_position = GetHomopolymerLengthIns(GetReverseComplement(sorted_list[0][5]), int(sorted_list[j][7]), int(sorted_list[0][0]), int(len(sorted_list[0][5])/2) + 1, DNA_strand)[1][1] 
                            else:
                                DNA_strand = "+"
                                end_position = GetHomopolymerLengthIns(sorted_list[0][5], int(sorted_list[0][7]), int(sorted_list[0][0]), int(len(sorted_list[0][5])/2) + 1, DNA_strand)[1][1] 
                else:    
                            if sorted_list[0][4] == '-':
                                end_position = GetHomopolymerLength(GetReverseComplement(sorted_list[0][5]), int(sorted_list[0][7]), int(sorted_list[0][0]), int(len(sorted_list[0][5])/2) + 1)[1][1] 
                            else:
                                end_position = GetHomopolymerLength(sorted_list[0][5], int(sorted_list[0][7]), int(sorted_list[0][0]), int(len(sorted_list[0][5])/2) + 1)[1][1] 

                if int(sorted_list[1][0]) <= (end_position + 10) and sorted_list[0][1] == sorted_list[1][1]:
                    if (sorted_list[0][6], sorted_list[0][1], sorted_list[0][6]) not in clusters.keys() and (sorted_list[0][8] == 'SNV' or sorted_list[0][8] == 'SNP'):
                        clusters[(sorted_list[0][6], sorted_list[0][1], sorted_list[0][6])] = [sorted_list[0]]
                    if (sorted_list[1][6], sorted_list[1][1], sorted_list[1][6]) not in clusters.keys() and (sorted_list[1][8] == 'SNV' or sorted_list[1][8] == 'SNP'):
                        clusters[(sorted_list[1][6], sorted_list[1][1], sorted_list[1][6])] = [sorted_list[1]]
                
                if ('-' not in sorted_list[1][3] and sorted_list[0][8] == "DIP") or sorted_list[1][8] == "Insertion":
                    if sorted_list[1][4] == '-':
                        DNA_strand = '-'
                        start_position = GetHomopolymerLengthIns(GetReverseComplement(sorted_list[1][5]), int(sorted_list[1][7]), int(sorted_list[1][0]), int(len(sorted_list[1][5])/2) + 1, DNA_strand)[1][0]
                    else:
                        DNA_strand = '+'
                        start_position = GetHomopolymerLengthIns(sorted_list[1][5], int(sorted_list[1][7]), int(sorted_list[1][0]), int(len(sorted_list[1][5])/2) + 1, DNA_strand)[1][0] 
                else:    
                    if sorted_list[1][4] == '-':
                        start_position = GetHomopolymerLength(GetReverseComplement(sorted_list[1][5]), int(sorted_list[1][7]), int(sorted_list[1][0]), int(len(sorted_list[1][5])/2) + 1)[1][0] 
                    else:
                        start_position = GetHomopolymerLength(sorted_list[1][5], int(sorted_list[1][7]), int(sorted_list[1][0]), int(len(sorted_list[1][5])/2) + 1)[1][0] 
                
                
                #start_position = GetHomopolymerLength(sorted_list[1][5], int(sorted_list[1][7]), int(sorted_list[1][0]), int(len(sorted_list[1][5])/2) + 1)[1][0]
                if start_position <= (sorted_list[0][0] + 10) and sorted_list[0][1] == sorted_list[1][1]:
                    if (sorted_list[0][6], sorted_list[0][1], sorted_list[0][6]) not in clusters.keys() and (sorted_list[0][8] == 'SNV' or sorted_list[0][8] == 'SNP'):
                        clusters[(sorted_list[0][6], sorted_list[0][1], sorted_list[0][6])] = [sorted_list[0]]
                    if (sorted_list[1][6], sorted_list[1][1], sorted_list[1][6]) not in clusters.keys() and (sorted_list[1][8] == 'SNV' or sorted_list[1][8] == 'SNP'):
                        clusters[(sorted_list[1][6], sorted_list[1][1], sorted_list[1][6])] = [sorted_list[1]]
            if sorted_list[len (sorted_list) - 2][1] == sorted_list[len (sorted_list) - 1][1] and (int(sorted_list[len (sorted_list) - 1][0]) <= int((sorted_list[ len (sorted_list) - 2][0])) + 10):
                if (sorted_list[len(sorted_list) - 1][6], sorted_list[len(sorted_list) - 1][1], sorted_list[len(sorted_list) - 1][6]) not in clusters.keys() and (sorted_list[len(sorted_list) - 1][8] == 'SNV' or sorted_list[len(sorted_list) - 1][8] == 'SNP'):
                    clusters[(sorted_list[len(sorted_list) - 1][6], sorted_list[len(sorted_list) - 1][1], sorted_list[len(sorted_list) - 1][6])] = [sorted_list[len(sorted_list) - 1]]
                  
                if (sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][6]) not in clusters.keys() and (sorted_list[len(sorted_list) - 2][8] == 'SNV' or sorted_list[len(sorted_list) - 2][8] == 'SNP'):
                    clusters[(sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][6])] = [sorted_list[len(sorted_list) - 2]]
                    
            else:#if strand[0] == '+':
                #end_position = GetHomopolymerLength(sorted_list[len(sorted_list) - 2][5], int(sorted_list[len(sorted_list) - 2][7]), int(sorted_list[len(sorted_list) - 2][0]), int(len(sorted_list[len(sorted_list) - 2][5])/2) + 1)[1][1]
                if (sorted_list[len(sorted_list) - 2][8] == "DIP" and '-' not in sorted_list[len(sorted_list) - 2][3]) or sorted_list[len(sorted_list) - 2][8] == "Insertion":
                    if sorted_list[len(sorted_list) - 2][4] == '-':
                        DNA_strand = "-"
                        end_position = GetHomopolymerLengthIns(GetReverseComplement(sorted_list[len(sorted_list) - 2][5]), int(sorted_list[len(sorted_list) - 2][7]), int(sorted_list[len(sorted_list) - 2][0]), int(len(sorted_list[len(sorted_list) - 2][5])/2) + 1, DNA_strand)[1][1] 
                    else:
                        DNA_strand = "+"
                        end_position = GetHomopolymerLengthIns(sorted_list[len(sorted_list) - 2][5], int(sorted_list[len(sorted_list) - 2][7]), int(sorted_list[len(sorted_list) - 2][0]), int(len(sorted_list[len(sorted_list) - 2][5])/2) + 1, DNA_strand)[1][1] 
                else:    
                    if sorted_list[len(sorted_list) - 2][4] == '-':
                        end_position = GetHomopolymerLength(GetReverseComplement(sorted_list[len(sorted_list) - 2][5]), int(sorted_list[len(sorted_list) - 2][7]), int(sorted_list[len(sorted_list) - 2][0]), int(len(sorted_list[len(sorted_list) - 2][5])/2) + 1)[1][1] 
                    else:
                        end_position = GetHomopolymerLength(sorted_list[len(sorted_list) - 2][5], int(sorted_list[len(sorted_list) - 2][7]), int(sorted_list[len(sorted_list) - 2][0]), int(len(sorted_list[len(sorted_list) - 2][5])/2) + 1)[1][1] 
                
                
                if int(sorted_list[len(sorted_list) - 1][0]) <= (end_position + 10) and sorted_list[len(sorted_list) - 1][1] == sorted_list[len(sorted_list) - 2][1]:
                    if (sorted_list[len(sorted_list) - 1][6], sorted_list[len(sorted_list) - 1][1], sorted_list[len(sorted_list) - 1][6]) not in clusters.keys() and (sorted_list[len(sorted_list) - 1][8] == 'SNV' or sorted_list[len(sorted_list) - 1][8] == 'SNP'):
                        clusters[(sorted_list[len(sorted_list) - 1][6], sorted_list[len(sorted_list) - 1][1], sorted_list[len(sorted_list) - 1][6])] = [sorted_list[len(sorted_list) - 1]]
                    if (sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][6]) not in clusters.keys() and (sorted_list[len(sorted_list) - 2][8] == 'SNV' or sorted_list[len(sorted_list) - 2][8] == 'SNP'):
                        clusters[(sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][6])] = [sorted_list[len(sorted_list) - 2]]
                #start_position = GetHomopolymerLength(sorted_list[len(sorted_list) - 1][5], int(sorted_list[len(sorted_list) - 1][7]), int(sorted_list[len(sorted_list) - 1][0]), int(len(sorted_list[len(sorted_list) - 1][5])/2) + 1)[1][0]
                if (sorted_list[len(sorted_list) - 1][8] == "DIP" and '-' not in sorted_list[len(sorted_list) - 1][3]) or sorted_list[len(sorted_list) - 1][8] == "Insertion": 
                    if sorted_list[len(sorted_list) - 1][4] == '-':
                        DNA_strand = "-"
                        start_position = GetHomopolymerLengthIns(GetReverseComplement(sorted_list[len(sorted_list) - 1][5]), int(sorted_list[len(sorted_list) - 1][7]), int(sorted_list[len(sorted_list) - 1][0]), int(len(sorted_list[len(sorted_list) - 1][5])/2) + 1, DNA_strand)[1][0]
                    else:
                        DNA_strand = "+"
                        start_position = GetHomopolymerLengthIns(sorted_list[len(sorted_list) - 1][5], int(sorted_list[len(sorted_list) - 1][7]), int(sorted_list[len(sorted_list) - 1][0]), int(len(sorted_list[len(sorted_list) - 1][5])/2) + 1, DNA_strand)[1][0]    
                else:    
                    if sorted_list[len(sorted_list) - 1][4] == '-':
                        start_position = GetHomopolymerLength(GetReverseComplement(sorted_list[len(sorted_list) - 1][5]), int(sorted_list[len(sorted_list) - 1][7]), int(sorted_list[len(sorted_list) - 1][0]), int(len(sorted_list[len(sorted_list) - 1][5])/2) + 1)[1][0]
                    else:
                        start_position = GetHomopolymerLength(sorted_list[len(sorted_list) - 1][5], int(sorted_list[len(sorted_list) - 1][7]), int(sorted_list[len(sorted_list) - 1][0]), int(len(sorted_list[len(sorted_list) - 1][5])/2) + 1)[1][0]
                if start_position <= (sorted_list[len(sorted_list) - 2][0] + 10) and sorted_list[len(sorted_list) - 2][1] == sorted_list[len(sorted_list) - 1][1]:
                    if (sorted_list[len(sorted_list) - 1][6], sorted_list[len(sorted_list) - 1][1], sorted_list[len(sorted_list) - 1][6]) not in clusters.keys() and (sorted_list[len(sorted_list) - 1][8] == 'SNV' or sorted_list[len(sorted_list) - 1][8] == 'SNP'):
                        clusters[(sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][6])] = [sorted_list[len(sorted_list) - 1]]
                    if (sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][6]) not in clusters.keys() and (sorted_list[len(sorted_list) - 2][8] == 'SNV' or sorted_list[len(sorted_list) - 2][8] == 'SNP'):
                        clusters[(sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][6])] = [sorted_list[len(sorted_list) - 2]]
                   
            for j in range (1,len (sorted_list) - 1):
                    if (sorted_list[j][1] == sorted_list[j - 1][1] and (int(sorted_list[j][0]) <= int((sorted_list[j - 1][0])) + 10)) or (sorted_list[j + 1][1] == sorted_list[j][1] and (int(sorted_list[j + 1][0]) <= int((sorted_list[j][0])) + 10)): #or ((isolate_list[key][j + n][6], int(isolate_list[key][j + n][0]), isolate_list[key][j + n][1]) in homopolymers.keys() and (int(isolate_list[key][j + n][0]) <= int((isolate_list[key][j + n - 1][0])) + 10 + homopolymers[key])) : 
                            if (sorted_list[j][6], sorted_list[j][1], sorted_list[j][0]) not in clusters.keys() and (sorted_list[j][8] == 'SNV' or sorted_list[j][8] == 'SNP'):
                                #if j == 0 or isolate_list[key][j - n] not in clusters.keys():
                                clusters[(sorted_list[j][6], sorted_list[j][1], sorted_list[j][0])] = [sorted_list[j]]
                    else:
                        #end_position = GetHomopolymerLength(sorted_list[j][5], int(sorted_list[j][7]), int(sorted_list[j][0]), int(len(sorted_list[j][5])/2) + 1)[1][1]
                        if (sorted_list[j][8] == "DIP"  and '-' not in sorted_list[j][3])or sorted_list[j][8] == "Insertion": 
                            if sorted_list[j][4] == '-':
                                DNA_strand = "-"
                                end_position = GetHomopolymerLengthIns(GetReverseComplement(sorted_list[j][5]), int(sorted_list[j][7]), int(sorted_list[j][0]), int(len(sorted_list[j][5])/2) + 1, DNA_strand)[1][1] 
                            else:
                                DNA_strand = "+"
                                end_position = GetHomopolymerLengthIns(sorted_list[j][5], int(sorted_list[j][7]), int(sorted_list[j][0]), int(len(sorted_list[j][5])/2) + 1, DNA_strand)[1][1] 
                        else:    
                            if sorted_list[j][4] == '-':
                                end_position = GetHomopolymerLength(GetReverseComplement(sorted_list[j][5]), int(sorted_list[j][7]), int(sorted_list[j][0]), int(len(sorted_list[j][5])/2) + 1)[1][1] 
                            else:
                                end_position = GetHomopolymerLength(sorted_list[j][5], int(sorted_list[j][7]), int(sorted_list[j][0]), int(len(sorted_list[j][5])/2) + 1)[1][1] 
                        if int(sorted_list[j + 1][0]) <= (end_position + 10) and sorted_list[j][1] == sorted_list[j + 1][1]:
                            if (sorted_list[j][6], sorted_list[j][1], sorted_list[j][0]) not in clusters.keys() and (sorted_list[j][8] == 'SNV' or sorted_list[j][8] == 'SNP'):
                                clusters[(sorted_list[j][6], sorted_list[j][1], sorted_list[j][0])] = [sorted_list[j]]
                            if (sorted_list[j + 1][6], sorted_list[j + 1][1], sorted_list[j + 1][0]) not in clusters.keys() and (sorted_list[j + 1][8] == 'SNV' or sorted_list[j + 1][8] == 'SNP'):
                                clusters[(sorted_list[j + 1][6], sorted_list[j + 1][1], sorted_list[j + 1][0])] = [sorted_list[j + 1]]
                        #start_position = GetHomopolymerLength(sorted_list[j][5], int(sorted_list[j][7]), int(sorted_list[j][0]), int(len(sorted_list[j][5])/2) + 1)[1][0]
                        if (sorted_list[j][8] == "DIP"  and '-' not in sorted_list[j][3])or sorted_list[j][8] == "Insertion":
                            if sorted_list[j][4] == '-':
                                DNA_strand = "-"
                                start_position = GetHomopolymerLengthIns(GetReverseComplement(sorted_list[j][5]), int(sorted_list[j][7]), int(sorted_list[j][0]), int(len(sorted_list[j][5])/2) + 1, DNA_strand)[1][0] 
                            else:
                                DNA_strand = "+"
                                start_position = GetHomopolymerLengthIns(sorted_list[j][5], int(sorted_list[j][7]), int(sorted_list[j][0]), int(len(sorted_list[j][5])/2) + 1, DNA_strand)[1][0] 
                        elif int(sorted_list[j][0]) != int((sorted_list[j + 1][0])) and int(sorted_list[j][0]) != int((sorted_list[j - 1][0])):    
                            if sorted_list[j][4] == '-':
                                start_position = GetHomopolymerLength(GetReverseComplement(sorted_list[j][5]), int(sorted_list[j][7]), int(sorted_list[j][0]), int(len(sorted_list[j][5])/2) + 1)[1][0] 
                            else:
                                start_position = GetHomopolymerLength(sorted_list[j][5], int(sorted_list[j][7]), int(sorted_list[j][0]), int(len(sorted_list[j][5])/2) + 1)[1][0] 
                        if start_position <= (sorted_list[j - 1][0] + 10) and sorted_list[j - 1][1] == sorted_list[j][1]:
                            if (sorted_list[j - 1][6], sorted_list[j - 1][1], sorted_list[j - 1][0]) not in clusters.keys() and (sorted_list[j - 1][8] == 'SNV' or sorted_list[j - 1][8] == 'SNP'):
                                clusters[(sorted_list[j - 1][6], sorted_list[j - 1][1], sorted_list[j - 1][0])] = [sorted_list[j - 1]]
                            if (sorted_list[j][6], sorted_list[j][1], sorted_list[j][0]) not in clusters.keys() and (sorted_list[j][8] == 'SNV' or sorted_list[j][8] == 'SNP'):
                                clusters[(sorted_list[j][6], sorted_list[j][1], sorted_list[j][0])] = [sorted_list[j]]           
                    #elif (sorted_list[j][1] == sorted_list[j + 1][1]) and(sorted_list[j][6], int(sorted_list[j][0]), sorted_list[j][1]) in homopolymers.keys() and (int(sorted_list[j + 1][0]) <= int((sorted_list[j][0])) + 10 + (homopolymers[(sorted_list[j][6], int(sorted_list[j][0]), sorted_list[j][1])][0][1] - sorted_list[j][0])):
                        #if sorted_list[j] not in clusters.keys() and (sorted_list[j][7] == "DIP" or sorted_list[j][7] == "Complex DIP" or sorted_list[key][j][7] == "Deletion" or sorted_list[j][7] == "Insertion"):
                            #clusters[(sorted_list[j][6], sorted_list[j][1], sorted_list[j][0])] = [sorted_list[j]]
                        #if sorted_list[j + 1] not in clusters.keys() and (sorted_list[j + 1][7] == "DIP" or sorted_list[j + 1][7] == "Complex DIP" or sorted_list[j + 1][7] == "Deletion" or sorted_list[j + 1][7] == "Insertion"):
                            #clusters[(sorted_list[j + 1][6], sorted_list[j + 1][1], sorted_list[j + 1][0])] = [sorted_list[j + 1]]
                        
                    #elif (sorted_list[j][1] == sorted_list[j - 1][1]) and (sorted_list[j][6], sorted_list[j][0], sorted_list[j][1]) in homopolymers.keys() and (int(sorted_list[j][0]) <= int((sorted_list[j-1][0])) + 10 + (sorted_list[j][0] - homopolymers[sorted_list[j][6], int(sorted_list[j][0]), sorted_list[j][1]][0][0] )):
                        #if sorted_list[j] not in clusters.keys() and (sorted_list[j][7] == "DIP" or sorted_list[j][7] == "Complex DIP" or sorted_list[j][7] == "Deletion" or sorted_list[j][7] == "Insertion"):
                            #clusters[(sorted_list[j][6], sorted_list[j][1], sorted_list[j][0])] = [sorted_list[j]]
                        #if sorted_list[j - 1] not in clusters.keys() and (sorted_list[j - 1][7] == "DIP" or sorted_list[j - 1][7] == "Complex DIP" or sorted_list[j - 1][7] == "Deletion" or sorted_list[j - 1][7] == "Insertion"):
                            #clusters[(sorted_list[j - 1][6], sorted_list[j - 1][1], sorted_list[j - 1][0])] = [sorted_list[j - 1]]
        
   

    complex_count = []
    for key in homopolymers.keys():
        if key in clusters.keys():
            complex_count.append(key)
    
    return clusters

clusters = FindClusters(WT_isolate_list, WT_homopolymers)
clusters16 = FindClusters(rad16_isolate_list, rad16_homopolymers)
clusters26 = FindClusters(rad26_isolate_list, rad26_homopolymers)
clustersB = FindClusters(WT_isolate_listB, WT_homopolymersB)
clusters16B = FindClusters(rad16_isolate_listB, rad16_homopolymersB)
clusters26B = FindClusters(rad26_isolate_listB, rad26_homopolymersB)


a = open('WT_sub_clustersv2.bed', 'w+')
b = open('rad16_sub_clustersv2.bed', 'w+')
d = open('rad26_sub_clustersv2.bed', 'w+')


e = open('WT_sub_clustersUVB.bed', 'w+')
f = open('rad16_sub_clustersUVB.bed', 'w+')
g = open('rad26_sub_clustersUVB.bed', 'w+')


def index2(tuple):
    return tuple[2]

def MakeClusterFiles(clusters, file1):
    for c in clusters:
        file1.write(clusters[c][0][1] + '\t' + str((clusters[c][0][0] - 1)) + '\t' + str(clusters[c][0][0]) + '\t' + clusters[c][0][2] + '\t' + clusters[c][0][3] + '\t' + clusters[c][0][4])
        file1.write('\n')
    

                
MakeClusterFiles(clusters,  a)
MakeClusterFiles(clusters16, b)
MakeClusterFiles(clusters26, d)
MakeClusterFiles(clustersB, e)
MakeClusterFiles(clusters16B, f)
MakeClusterFiles(clusters26B, g)


a.close()
b.close()
d.close()
e.close()
f.close()
g.close()
