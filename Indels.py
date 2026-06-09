




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




def GetSequence(file): #Return reference chromosome sequence from the SacCer3 genome
    f1 = open(file)
    sequence1 = f1.read()
    sequence1 = str(sequence1.strip().upper())
    sequence1 = ''.join(sequence1.split('\n'))
    return sequence1


sequences = []
sequence1 = GetSequence('chr1.txt')
sequences.append(sequence1)
print(len(sequence1))
sequence2 = GetSequence('chr2.txt')
sequences.append(sequence2)
print(len(sequence2))
sequence3 = GetSequence('chr3.txt')
sequences.append(sequence3)
print(len(sequence3))
sequence4 = GetSequence('chr4.txt')
sequences.append(sequence4)
print(len(sequence4))
sequence5 = GetSequence('chr5.txt')
sequences.append(sequence5)
print(len(sequence5))
sequence6 = GetSequence('chr6.txt')
print(len(sequence6))
sequences.append(sequence6)
sequence7 = GetSequence('chr7.txt')
sequences.append(sequence7)
print(len(sequence7))
sequence8 = GetSequence('chr8.txt')
sequences.append(sequence8)
print(len(sequence8))
sequence9 = GetSequence('chr9.txt')
sequences.append(sequence9)
print(len(sequence9))
sequence10 = GetSequence('chr10.txt')
sequences.append(sequence10)
print(len(sequence10))
sequence11 = GetSequence('chr11.txt')
print(len(sequence11))
sequences.append(sequence11)
sequence12 = GetSequence('chr12.txt')
print(len(sequence12))
sequences.append(sequence12)
sequence13 = GetSequence('chr13.txt')
sequences.append(sequence13)
print(len(sequence13))
sequence14 = GetSequence('chr14.txt')
sequences.append(sequence14)
print(len(sequence14))
sequence15 = GetSequence('chr15.txt')
sequences.append(sequence15)
print(len(sequence15))
sequence16 = GetSequence('chr16.txt')
sequences.append(sequence16)
print(len(sequence16))


#Returns the start and end positions and the length of the homopolymer containing a deleted base
def GetHomopolymerLength(sequence, mut_pos1, mut_pos, mut_position):
    #mut_position = int(len(sequence)/2) + 1
    mut_pos1 = int(mut_pos1) #0 index position on chromosome
    mut_pos = int(mut_pos) #1 index position on chromosome
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
               
                if mut_position <= end_position  and mut_position >= (end_position - (length - 1 )):
                     # or [(mut_position2 <= end_position)  and (mut_position2 >= end_position - (length - 1 ))]:
                     lengths[(end_position  + (mut_pos - mut_position) - (length - 1 ) , end_position  + (mut_pos - mut_position))] = length
                elif mut_position2 <= end_position  and mut_position2 >= (end_position - (length - 1 )):
                    lengths[(end_position  + (mut_pos - mut_position) - (length - 1 ) , end_position  + (mut_pos - mut_position))] = length
    
        if mut_position <= end_position  and mut_position >= (end_position - (length - 1 )):
        
    
            lengths[(end_position  + (mut_pos - mut_position) - (length - 1 ) , end_position  + (mut_pos - mut_position) )] = length
        elif mut_position2 <= end_position  and mut_position2 >= (end_position - (length - 1 )):
                    lengths[(end_position  + (mut_pos - mut_position) - (length - 1 ) , end_position  + (mut_pos - mut_position))] = length   
    if len(lengths) == 0:
        return 1, (mut_pos, mut_pos)
    for key in lengths:
        if lengths[key] == max(lengths.values()):
            return lengths[key], key

#Returns the start and end positions and the length of the homopolymer containing an inserted base
def GetHomopolymerLengthIns(sequence, mut_pos1, mut_pos, mut_position, strand):
    #mut_position = int(len(sequence)/2) + 1
    mut_position = mut_position 
    mut_pos1 = int(mut_pos1) #0 index position on chromosome
    mut_pos = int(mut_pos) #1 index position on chromosome
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
               
                if mut_position <= end_position + 1  and (mut_position) >= (end_position - (length - 1)):
                        lengths[(end_position  + (mut_pos - mut_position) - (length - 1 ) , end_position  + (mut_pos - mut_position))] = length
              
    
            if mut_position <= end_position + 1   and mut_position >= (end_position - (length - 1)):
                lengths[(end_position  + (mut_pos - mut_position) - (length - 1 ) , end_position  + (mut_pos - mut_position))] = length
         
    
    if len(lengths) == 0:
        return 1, (mut_pos, mut_pos)
    for key in lengths:
        if lengths[key] == max(lengths.values()):
            return lengths[key], key



with open("UVC_filter_mutations.txt") as f: #Input UVC whole-genome sequencing data
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
trinucleotide = []
strand = []
isolate = []
isolate_list = {}
lines.pop(0)
isolate_keys = []
WT_double = 0
rad16_double = 0
rad26_double = 0
rad30_double = 0
WT_subs = 0
rad16_subs = 0
rad26_subs = 0
rad30_subs = 0

for line in lines:

    if 'Mito' not in line[6] and 'ySR128' not in line[3] and line[0] != 'dRad16_25':
        isolate.append(line[0])
    
        genotype.append(line[4])
    
        mutation_type.append(line[10])
       
       
        position2.append(int(line[9]))
        length.append(line[11])
        position1.append(int(line[9]) - len(line[13]))
        allele.append(line[12])
        chromosome.append('chr' + str(line[6]))

        mutation.append(line[13])
                    
        
#Count mutations of different types
for g in range(0, len(genotype)):
    if 'WT' in genotype[g]:
        genotypes['WT'] = genotypes['WT'] + 1
    if 'rad16' in genotype[g]:
            genotypes['rad16'] = genotypes['rad16'] + 1
    if 'rad26' in genotype[g]:
            genotypes['rad26'] = genotypes['rad26'] + 1
    if 'rad30' in genotype[g]:
            genotypes['rad30'] = genotypes['rad30'] + 1
    if 'WT' in genotype[g]:
        if str(chromosome[g]) in WT_chromosomes.keys():
                WT_chromosomes[str(chromosome[g])] = WT_chromosomes[str(chromosome[g])] + 1

        else:
            WT_chromosomes[str(chromosome[g])] = 1
        if mutation_type[g] == "DIP":
            WT_indels = WT_indels + 1
        if mutation_type[g] == "DNP" or mutation_type[g] == "DNV":
            WT_double = WT_double + 1
        if mutation_type[g] == "SNP" or mutation_type[g] == "SNV":
            WT_subs = WT_subs + 1
           
    if 'rad16' in genotype[g]:
            if str(chromosome[g]) in rad16_chromosomes.keys():
                rad16_chromosomes[chromosome[g]] = rad16_chromosomes[chromosome[g]] + 1
            else:
                rad16_chromosomes[chromosome[g]] = 1
            if mutation_type[g] == "DIP":
                rad16_indels = rad16_indels + 1
            if mutation_type[g] == "DNP" or mutation_type[g] == "DNV":
                rad16_double = rad16_double + 1
            if mutation_type[g] == "SNP" or mutation_type[g] == "SNV":
                rad16_subs = rad16_subs + 1
    if 'rad26' in genotype[g]:
            if str(chromosome[g]) in rad26_chromosomes.keys():
                rad26_chromosomes[str(chromosome[g])] = rad26_chromosomes[str(chromosome[g])] + 1
            else:
                rad26_chromosomes[str(chromosome[g])] = 1
            if mutation_type[g] == "DIP":
                rad26_indels = rad26_indels + 1
            if mutation_type[g] == "DNV" or mutation_type[g] == "DNP":
                rad26_double = rad26_double + 1
            if mutation_type[g] == "SNV" or mutation_type[g] == "SNP":
                rad26_subs = rad26_subs + 1
    if 'rad30' in genotype[g]:
            if str(chromosome[g]) in rad30_chromosomes.keys():
                rad30_chromosomes[str(chromosome[g])] = rad30_chromosomes[str(chromosome[g])] + 1
            else:
                rad30_chromosomes[str(chromosome[g])] = 1
            if mutation_type[g] == "DIP":
                rad30_indels = rad30_indels + 1
            if mutation_type[g] == "SNV" or mutation_type[g] == "SNP":
                rad30_subs = rad30_subs + 1
            if mutation_type[g] == "DNV" or mutation_type[g] == "DNP":
                rad30_double = rad30_double + 1
        

for pos in range(0, len(position2)):  #Find trinucleotide sequence contexts
    
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
            isolate_list[(isolate[i1], genotype[i1])] = [(position2[i1], str(chromosome[i1]), trinucleotide[i1], mutation[i1], strand[i1], mutation_type[i1], isolate[i1], allele[i1])]
    else:
        isolate_list[(isolate[i1], genotype[i1])] .append((position2[i1], str(chromosome[i1]), trinucleotide[i1], mutation[i1], strand[i1], mutation_type[i1], isolate[i1], allele[i1]))
    
def MutationsPerIsolate( isolate_list, file1, file2, file3):  #Finds the number of deletions and insertions for each isolate and calculates the average 
    #mutations per isolate for WT, rad16, rad26, and rad30
    Deletions = {'WT': [], 'rad16': [], 'rad26': [], 'rad30': []}
    Insertions = {'WT': [], 'rad16': [], 'rad26': [], 'rad30': []}
    All = {'WT': [], 'rad16': [], 'rad26': [], 'rad30': []}
    average_WT_deletion = 0
    average_rad16_deletion = 0
    average_rad26_deletion = 0
    average_rad30_deletion = 0
    average_WT_insertion = 0
    average_rad16_insertion = 0
    average_rad26_insertion = 0
    average_rad30_insertion = 0

    for key in isolate_list.keys():

        d_count = 0
        i_count = 0
        count = 0

        for a in range(0, len(isolate_list[key])):
            if isolate_list[key][a][5] == 'DIP' or isolate_list[key][a][5] == 'Deletion' or isolate_list[key][a][5] == 'Insertion':
            #if isolate_list[key][a][5] == 'SNP' or isolate_list[key][a][5] == 'SNV':

                count = count + 1
                if '-' in isolate_list[key][a][3]:# and  len(isolate_list[key][a][7]) == 1:
                    d_count = d_count + 1
                elif '-' not in isolate_list[key][a][3]:# and len(isolate_list[key][a][3]) == 1:
                    i_count = i_count + 1
        if key[1] == 'WT':
            Deletions['WT'].append(d_count)
            Insertions['WT'].append(i_count)
            All['WT'].append(count)
        if 'rad16' in key[1] or 'Rad16' in key[1]:
            Deletions['rad16'].append(d_count)
            Insertions['rad16'].append(i_count)
            All['rad16'].append(count)
        if 'rad26' in key[1] or 'Rad26' in key[1]:
            Deletions['rad26'].append(d_count)
            Insertions['rad26'].append(i_count)
            All['rad26'].append(count)
        if 'rad30' in key[1] or 'Rad30' in key[1]:
            Deletions['rad30'].append(d_count)
            Insertions['rad30'].append(i_count) 
            All['rad30'].append(count) 
    file1.write('Deletions: ')
    file1.write('\n')
    file2.write('Deletions: ')
    file2.write('\n')
    file3.write('Deletions: ')
    file3.write('\n')
    for b in Deletions['WT']:
        average_WT_deletion = average_WT_deletion+ b
        file1.write(str(b))
        file1.write('\n')
    average_WT_deletion = average_WT_deletion/len(Deletions['WT'])
    for c in Deletions['rad16']:
        file2.write(str(c))
        file2.write('\n')
        average_rad16_deletion = average_rad16_deletion+ c
    average_rad16_deletion = average_rad16_deletion/(len(Deletions['rad16']))
    for d in Deletions['rad26']:
        file3.write(str(d))
        file3.write('\n')
        average_rad26_deletion = average_rad26_deletion+ d
    average_rad26_deletion = average_rad26_deletion/len(Deletions['rad26'])
    for e in Deletions['rad30']:
        average_rad30_deletion = average_rad30_deletion+ e
    average_rad30_deletion = average_rad30_deletion/len(Deletions['rad30'])

    file1.write('Insertions: ')
    file1.write('\n')
    file2.write('Insertions: ')
    file2.write('\n')
    file3.write('Insertions: ')
    file3.write('\n')

    for f in Insertions['WT']:
        file1.write(str(f))
        file1.write('\n')
        average_WT_insertion = average_WT_insertion + f
    average_WT_insertion = average_WT_insertion/len(Insertions['WT'])
    for g in Insertions['rad16']:
        file2.write(str(g))
        file2.write('\n')
        average_rad16_insertion = average_rad16_insertion + g
    average_rad16_insertion = average_rad16_insertion/len(Insertions['rad16'])
    for h in Insertions['rad26']:
        file3.write(str(h))
        file3.write('\n')
        average_rad26_insertion = average_rad26_insertion + h
    average_rad26_insertion = average_rad26_insertion/len(Insertions['rad26'])
    for i in Insertions['rad30']:
        average_rad30_insertion = average_rad30_insertion + i
    average_rad30_insertion = average_rad30_insertion/len(Insertions['rad30'])
    print(Deletions)
    print(Insertions)
    print(All)
    
    return Deletions, Insertions, All
       
file1 = open('WT_mutationsperisolate', 'w+')
file2 = open('rad16_mutationsperisolate', 'w+')
file3 = open('rad26_mutationsperisolate', 'w+')
MutationsPerIsolate(isolate_list, file1, file2, file3)


WT_deletions = {}
WT_insertions = {} 
rad16_deletions = {}
rad16_insertions = {} 
rad26_deletions = {}
rad26_insertions = {}
rad30_deletions = {}  
rad30_insertions = {}   

#Makes a dictionary with the number of deletions and insertions of each length
for a in range(0, len(allele)):
    if mutation_type[a] == 'DIP':
        if genotype[a] == 'WT':      
                if '-' in allele[a]:
                    if len(mutation[a]) in WT_insertions.keys():
                        WT_insertions[len(mutation[a])] = WT_insertions[len(mutation[a])] + 1
                    else: 
                        WT_insertions[len(mutation[a])] = 1
                elif '-' in mutation[a]:
                    if len(allele[a]) in WT_deletions.keys():
                        WT_deletions[len(allele[a])] = WT_deletions[len(allele[a])] + 1
                    else: 
                        WT_deletions[len(allele[a])] = 1
        if 'rad16' in genotype[a]:      
                if '-' in allele[a]:
                    if len(mutation[a]) in rad16_insertions.keys():
                        rad16_insertions[len(mutation[a])] = rad16_insertions[len(mutation[a])] + 1
                    else: 
                        rad16_insertions[len(mutation[a])] = 1
                elif '-' in mutation[a]:
                    if len(allele[a]) in rad16_deletions.keys():
                        rad16_deletions[len(allele[a])] = rad16_deletions[len(allele[a])] + 1
                    else: 
                        rad16_deletions[len(allele[a])] = 1
        if 'rad26' in genotype[a]:      
                if '-' in allele[a]:
                    if len(mutation[a]) in rad26_insertions.keys():
                        rad26_insertions[len(mutation[a])] = rad26_insertions[len(mutation[a])] + 1
                    else: 
                        rad26_insertions[len(mutation[a])] = 1
                elif '-' in mutation[a]:
                    if len(allele[a]) in rad26_deletions.keys():
                        rad26_deletions[len(allele[a])] = rad26_deletions[len(allele[a])] + 1
                    else: 
                        rad26_deletions[len(allele[a])] = 1
        if 'rad30' in genotype[a]:      
                if '-' in allele[a]:
                    if len(mutation[a]) in rad30_insertions.keys():
                        rad30_insertions[len(mutation[a])] = rad30_insertions[len(mutation[a])] + 1
                    else: 
                        rad30_insertions[len(mutation[a])] = 1
                elif '-' in mutation[a]:
                    if len(allele[a]) in rad30_deletions.keys():
                        rad30_deletions[len(allele[a])] = rad30_deletions[len(allele[a])] + 1
                    else: 
                        rad30_deletions[len(allele[a])] = 1


f1a = open("WT_insertionsv3_sorted_tandem.bed", 'w+')
f1b = open("rad16_insertionsv3_sorted_tandem.bed", 'w+')
f1c = open("rad26_insertionsv3_sorted_tandem.bed", 'w+')
f1d = open("rad30_insertionsv3_sorted_tandem.bed", 'w+')
f1 = open("LOH_DIPsv2_sorted_tandem.bed", 'w+')
f2 = open("UV_yeast_Merged_SNP_DIPsv2_sorted_tandem_insertions.bed", 'w+')
f2a = open("WT_deletionsv3_sorted_tandem.bed", 'w+')
f2b = open("rad16_deletionsv3_sorted_tandem.bed", 'w+')
f2c = open("rad26_deletionsv3_sorted_tandem.bed", 'w+')
f2d = open("rad30_deletionsv3_sorted_tandem.bed", 'w+')
f3a = open("WT_allv3_sorted_tandem.bed", 'w+')
f3b = open("rad16_allv3_sorted_tandem.bed", 'w+')
f3c = open("rad26_allv3_sorted_tandem.bed", 'w+')
f3d = open("rad30_allv3_sorted_tandem.bed", 'w+')

f4a = open("WT_allcombined_sorted_tandem.bed", 'w+')
f4b = open("rad16_allcombined_sorted_tandem.bed", 'w+')
f4c = open("rad26_allcombined_sorted_tandem.bed", 'w+')
f4d = open("rad30_allcombined_sorted_tandem.bed", 'w+')


#Makes bed files for the indels, deletions, and insertion for each strain
def WriteChr(chr_name, file1a, file1b, file1c, file1d,file2a, file2b, file2c, file2d, file3a, file3b, file3c, file3d, chromosome, position1, position2,
             trinucleotide, allele, mutation, genotype, mutation_type ,isolate, strand):
    file_data = []
    chromosome1_dict = {}
    chromosome1 =[]
    for x in range(0, len(mutation_type)):
        if mutation_type[x] == 'DIP' or mutation_type[x] == 'Insertion' or mutation_type[x] == 'Deletion':
                
                if trinucleotide[x].isalpha():
                    if chromosome[x] == chr_name: 
                        chromosome1_dict[int(position2[x]), isolate[x]] = x
                


    for key1 in chromosome1_dict.keys():
        chromosome1.append((int(key1[0]), key1[1]))
    sorted_chromosome1= sorted(chromosome1)
    for y in range(0, len(sorted_chromosome1)):
        file_data.append(chromosome[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + str(int(position1[chromosome1_dict[sorted_chromosome1[y]]])) + '\t' + str(position2[chromosome1_dict[sorted_chromosome1[y]]]) + '\t' + trinucleotide[chromosome1_dict[sorted_chromosome1[y]]]  + '\t' + mutation[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + strand[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + allele[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + isolate[chromosome1_dict[sorted_chromosome1[y]]])
    for fd in range(0, len(file_data)):
        if '-' in allele[chromosome1_dict[sorted_chromosome1[fd]]]: #Insertions
            if genotype[chromosome1_dict[sorted_chromosome1[fd]]] == 'WT':
                file1a.write(file_data[fd])
                file3a.write(file_data[fd])
                file1a.write('\n')
                file3a.write('\n')
            if 'rad16' in genotype[chromosome1_dict[sorted_chromosome1[fd]]] or 'Rad16' in  genotype[chromosome1_dict[sorted_chromosome1[fd]]]:
                file1b.write(file_data[fd])
                file3b.write(file_data[fd])
                file1b.write('\n')
                file3b.write('\n')
            if 'rad26' in genotype[chromosome1_dict[sorted_chromosome1[fd]]] or 'Rad26' in genotype[chromosome1_dict[sorted_chromosome1[fd]]]:
                file1c.write(file_data[fd])
                file3c.write(file_data[fd])
                file1c.write('\n')
                file3c.write('\n')
            if 'rad30' in genotype[chromosome1_dict[sorted_chromosome1[fd]]] or 'Rad30' in  genotype[chromosome1_dict[sorted_chromosome1[fd]]]:
                file1d.write(file_data[fd])
                file3d.write(file_data[fd])
                file1d.write('\n')
                file3d.write('\n')

        elif '-' in mutation[chromosome1_dict[sorted_chromosome1[fd]]]:#Deletions
            if genotype[chromosome1_dict[sorted_chromosome1[fd]]] == 'WT':
                file2a.write(file_data[fd])
                file3a.write(file_data[fd])
                file2a.write('\n')
                file3a.write('\n')
                
            if 'rad16' in  genotype[chromosome1_dict[sorted_chromosome1[fd]]] or 'Rad16' in  genotype[chromosome1_dict[sorted_chromosome1[fd]]]:
                file2b.write(file_data[fd])
                file3b.write(file_data[fd])
                file2b.write('\n')
                file3b.write('\n')
            if 'rad26' in  genotype[chromosome1_dict[sorted_chromosome1[fd]]] or 'Rad26' in genotype[chromosome1_dict[sorted_chromosome1[fd]]]:
                file2c.write(file_data[fd])
                file3c.write(file_data[fd])
                file2c.write('\n')
                file3c.write('\n')
            if 'rad30' in  genotype[chromosome1_dict[sorted_chromosome1[fd]]] or 'Rad30' in  genotype[chromosome1_dict[sorted_chromosome1[fd]]]:
                file2d.write(file_data[fd])
                file3d.write(file_data[fd])
                file2d.write('\n')
                file3d.write('\n')
        


def WriteAllChr(chr_name, file1a, file1b, file1c, file1d, chromosome, position1, position2,
             trinucleotide, allele, mutation, genotype, mutation_type ,isolate, strand):
    #Makes bed files with all the mutations
    file_data = []
    chromosome1_dict = {}
    chromosome1 =[]
    for x in range(0, len(mutation_type)):
                
                    if chromosome[x] == chr_name: 
                        chromosome1_dict[int(position2[x]), isolate[x]] = x
                


    for key1 in chromosome1_dict.keys():
        chromosome1.append((int(key1[0]), key1[1]))
    sorted_chromosome1= sorted(chromosome1)
    for y in range(0, len(sorted_chromosome1)):
        file_data.append(chromosome[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + str(int(position1[chromosome1_dict[sorted_chromosome1[y]]])) + '\t' + str(position2[chromosome1_dict[sorted_chromosome1[y]]]) + '\t' + trinucleotide[chromosome1_dict[sorted_chromosome1[y]]]  + '\t' + mutation[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + strand[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + allele[chromosome1_dict[sorted_chromosome1[y]]] + '\t' + isolate[chromosome1_dict[sorted_chromosome1[y]]])
        
    for fd in range(0, len(file_data)):
            if genotype[chromosome1_dict[sorted_chromosome1[fd]]] == 'WT':
                file1a.write(file_data[fd])
               
                file1a.write('\n')
    
            if 'rad16' in genotype[chromosome1_dict[sorted_chromosome1[fd]]] or 'Rad16' in  genotype[chromosome1_dict[sorted_chromosome1[fd]]]:
                file1b.write(file_data[fd])
                file1b.write('\n')
            if 'rad26' in genotype[chromosome1_dict[sorted_chromosome1[fd]]] or 'Rad26' in genotype[chromosome1_dict[sorted_chromosome1[fd]]]:
                file1c.write(file_data[fd])
                file1c.write('\n')
            if 'rad30' in genotype[chromosome1_dict[sorted_chromosome1[fd]]] or 'Rad30' in  genotype[chromosome1_dict[sorted_chromosome1[fd]]]:
                file1d.write(file_data[fd])
                file1d.write('\n')
            
            



WriteChr('chrI',f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, f3a, f3b, f3c, f3d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrII', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, f3a, f3b, f3c, f3d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrIII',  f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, f3a, f3b, f3c, f3d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrIV', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, f3a, f3b, f3c, f3d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrIX', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, f3a, f3b, f3c, f3d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrV', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, f3a, f3b, f3c, f3d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrVI', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, f3a, f3b, f3c, f3d,chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrVII', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, f3a, f3b, f3c, f3d,  chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrVIII', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, f3a, f3b, f3c, f3d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrX', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, f3a, f3b, f3c, f3d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXI', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, f3a, f3b, f3c, f3d,  chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXII', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, f3a, f3b, f3c, f3d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXIII', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, f3a, f3b, f3c, f3d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXIV', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, f3a, f3b, f3c, f3d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand) 
WriteChr('chrXV', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, f3a, f3b, f3c, f3d,  chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteChr('chrXVI', f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, f3a, f3b, f3c, f3d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)


WriteAllChr('chrI',f4a, f4b, f4c, f4d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrII', f4a, f4b, f4c, f4d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrIII',  f4a, f4b, f4c, f4d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrIV', f4a, f4b, f4c, f4d,  chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrIX', f4a, f4b, f4c, f4d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrV', f4a, f4b, f4c, f4d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrVI', f4a, f4b, f4c, f4d ,chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrVII', f4a, f4b, f4c, f4d,  chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrVIII', f4a, f4b, f4c, f4d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrX', f4a, f4b, f4c, f4d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrXI', f4a, f4b, f4c, f4d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrXII', f4a, f4b, f4c, f4d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrXIII', f4a, f4b, f4c, f4d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrXIV', f4a, f4b, f4c, f4d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand) 
WriteAllChr('chrXV', f4a, f4b, f4c, f4d,  chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)
WriteAllChr('chrXVI', f4a, f4b, f4c, f4d, chromosome, position1, position2,
trinucleotide, allele, mutation, genotype, mutation_type, isolate, strand)

f1.close()
f1a.close()
f1b.close()
f1c.close()
f1d.close()
f2.close()
f2a.close()
f2b.close()
f2c.close()
f2d.close()
f3a.close()
f3b.close()
f3c.close()
f3d.close()
f4a.close()
f4b.close()
f4c.close()
f4d.close()


def MakeTriNucFile(file1, file2, file3, file4, lengths): #Makes bed files with single-base deletions
    file1.seek(0)
    new_data = file1.read()
    new_lines = [s.strip().split() for s in new_data.splitlines()]
    n = -1
    for line in new_lines:
            n = n + 1
            if len(line[6]) == 1:
               
                    file2.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' + line[5]  + '\t' + line[6])
                    file2.write('\n')
                
                    if lengths[n] >= 4: #Deletions in homopolymers
                        file3.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] +  '\t' + line[4] + '\t' + line[5] + '\t' + line[6])
                        file3.write('\n')
                    else: #Deletions not in homopolymers
                        file4.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4]+ '\t' + line[5] + '\t' + line[6])
                        file4.write('\n')
              

def MakeInsertionFile(file1, file2, file3, file4, lengths):  #Makes bed files for single-base insertions
    file1.seek(0)
    new_data = file1.read()
    new_lines = [s.strip().split() for s in new_data.splitlines()]
    n = -1
    for line in new_lines:
            n = n + 1
            if len(line[4]) == 1:
                
                file2.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' + line[5]  + '\t' + line[6])
                file2.write('\n')
               
                if lengths[n] >= 4: #Insertions in homopolymers
                        file3.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' + line[5]  + '\t' + line[6])
                        file3.write('\n')
                else: #Insertions not in homopolymers
                        file4.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' + line[5]  + '\t' + line[6])
                        file4.write('\n')

               

f1a = open("WT_insertionsv3_sorted_tandem.bed")
f1b = open("rad16_insertionsv3_sorted_tandem.bed")
f1c = open("rad26_insertionsv3_sorted_tandem.bed")
f1e = open("WT_singlebase_insertions.bed", 'w+')
f1f = open("rad16_singlebase_insertions.bed", 'w+')
f1g = open("rad26_singlebase_insertions.bed", 'w+')





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

    
    
    return multinucleotide, new_strand, new_chromosome, new_position, new_isolate, new_trinucleotide, new_mutation, new_allele, lengths, new_position1

sequencefile1 = open("WT_single_deletion_sequences_T.txt", 'w+')
sequencefile2 = open("WT_single_insertion_sequences_T.txt", 'w+')
sequencefile3 = open("rad16_single_deletion_sequences_T.txt", 'w+')
sequencefile4 = open("rad16_single_insertion_sequences_T.txt", 'w+') 
sequencefile5 = open("rad26_single_deletion_sequences_T.txt", 'w+')
sequencefile6 = open("rad26_single_insertion_sequences_T.txt", 'w+') 



def MakeSequenceFile(file1, file2, pentanucleotide, mutation, allele, strand):
    for p in range(0, len(pentanucleotide)):
        #if allele[p] == 'C' or allele[p] == 'G':
        if allele[p] == 'T' or allele[p] == 'A':
            if '-' in mutation[p] and len(allele[p]) == 1:
                if strand[p] == '-':
                    file1.write(GetComplement(GetReverseComplement(pentanucleotide[p])))
                    file1.write('\n')
                else:
                    file1.write(pentanucleotide[p])
                    file1.write('\n')
        #elif mutation[p]  == 'C' or mutation[p] == 'G':
        if mutation[p] == 'T' or mutation[p] == 'A':
            if '-' in allele[p]  and len(mutation[p]) == 1:
                if strand[p] == '-':
                    file2.write(GetComplement(GetReverseComplement(pentanucleotide[p])))
                    file2.write('\n')
                else:
                    file2.write(pentanucleotide[p])
                    file2.write('\n')

            
#Counts the number of mutations in homopolymers of different lengths
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
                if '-' in new_mutation[i] and '-' not in new_allele[i]:
                #if new_mutation[i] == '-':
                    homopolymer_deletions = homopolymer_deletions + 1
                    if key in deletion_numbers.keys():
                        deletion_numbers[key] = deletion_numbers[key] + 1
                    else:
                        deletion_numbers[key] = 1
                if '-' in new_allele[i] and '-' not in new_mutation[i]:
                #if new_allele[i] == '-':
                    homopolymer_insertions = homopolymer_insertions + 1
                    if key in insertion_numbers.keys():
                        insertion_numbers[key] = insertion_numbers[key] + 1
                    else:
                        insertion_numbers[key] = 1
                

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




a1 = open("WT_deletionsv3_sorted_tandem.bed")
a2 = open('WT_singlebase_deletions_Trinuc.bed', 'w+')
#a3 = open('WT_singlebase_deletions_Trinuc_dipyrimidine.bed', 'w+')
a3 = open('WT_singlebase_deletions_Trinuc_homopolymer.bed', 'w+')
#a4 = open('WT_singlebase_deletions_Trinuc_nondipyrimidine.bed', 'w+')
a4 = open('WT_singlebase_deletions_Trinuc_nonhomopolymer.bed', 'w+')
a5 = open('WT_sorted_deletions_Pentanuc.bed', 'w+')
a6 = open('WT_sorted_singledeletions_25nuc.bed', 'w+')
#a7 = open("WT_allv3_sorted_tandem.bed")
a7 = open("WT_allcombined_sorted_tandem.bed")
a8 = open('WT_sorted_combined_25nuc.bed', 'w+')

b1 = open("rad16_deletionsv3_sorted_tandem.bed")
b2 = open('rad16_singlebase_deletions_Trinuc.bed', 'w+')
#b3 = open('rad16_singlebase_deletions_Trinuc_dipyrimidine.bed', 'w+')
b3 = open('rad16_singlebase_deletions_Trinuc_homopolymer.bed', 'w+')
#b4 = open('rad16_singlebase_deletions_Trinuc_nondipyrimidine.bed', 'w+')
b4 = open('rad16_singlebase_deletions_Trinuc_nonhomopolymer.bed', 'w+')
b5 = open('rad16_sorted_deletions_Pentanuc.bed', 'w+')
b6 = open('rad16_sorted_single_deletions_25nuc.bed', 'w+')
#b7 = open("rad16_allv3_sorted_tandem.bed")
b7 = open("rad16_allcombined_sorted_tandem.bed")
b8 = open('rad16_sorted_combined_25nuc.bed', 'w+')

d1 = open("rad26_deletionsv3_sorted_tandem.bed")
d2 = open('rad26_singlebase_deletions_Trinuc.bed', 'w+')
#d3 = open('rad26_singlebase_deletions_Trinuc_dipyrimidine.bed', 'w+')
d3 = open('rad26_singlebase_deletions_Trinuc_homopolymer.bed', 'w+')
#d4 = open('rad26_singlebase_deletions_Trinuc_nondipyrimidine.bed', 'w+')
d4 = open('rad26_singlebase_deletions_Trinuc_nonhomopolymer.bed', 'w+')
d5 = open('rad26_sorted_deletions_Pentanuc.bed', 'w+')
d6 = open('rad26_sorted_single_deletions_25nuc.bed', 'w+')
#d7 = open("rad26_allv3_sorted_tandem.bed")
d7 = open("rad26_allcombined_sorted_tandem.bed")
d8 = open('rad26_sorted_combined_25nuc.bed', 'w+')

e1 = open("rad30_deletions_sorted_tandem.bed")
e2 = open('rad30_singlebase_deletions_Trinuc.bed', 'w+')
#e3 = open('rad30_singlebase_deletions_Trinuc_dipyrimidine.bed', 'w+')
#e4 = open('rad30_singlebase_deletions_Trinuc_nondipyrimidine.bed', 'w+')
e3 = open('rad30_singlebase_deletions_Trinuc_homopolymer.bed', 'w+')
e4 = open('rad30_singlebase_deletions_Trinuc_nonhomopolymer.bed', 'w+')
e5 = open('rad30_sorted_deletions_Pentanuc.bed', 'w+')
e6 = open('rad30_sorted_single_deletions_25nuc.bed', 'w+')

g1 = open("WT_insertionsv3_sorted_tandem.bed")
g2 = open('WT_singlebase_insertions_Trinuc.bed', 'w+')
#g3 = open('WT_singlebase_combined_Trinuc_dipyrimidine.bed', 'w+')
#g4 = open('WT_singlebase_combined_Trinuc_nondipyrimidine.bed', 'w+')
g3 = open('WT_singlebase_insertions_Trinuc_homopolymer.bed', 'w+')
g4 = open('WT_singlebase_insertions_Trinuc_nonhomopolymer.bed', 'w+')
#g5 = open('WT_sorted_insertions_Pentanuc.bed', 'w+')
g6 = open('WT_sorted_nonhomopolymer_singleinsertions_25nuc.bed', 'w+')

h1 = open("rad16_insertionsv3_sorted_tandem.bed")
h2 = open('rad16_singlebase_insertions_Trinuc.bed', 'w+')
#h3 = open('rad16_singlebase_combined_Trinuc_dipyrimidine.bed', 'w+')
#h4 = open('rad16_singlebase_combined_Trinuc_nondipyrimidine.bed', 'w+')
h3 = open('rad16_singlebase_insertions_Trinuc_homopolymer.bed', 'w+')
h4 = open('rad16_singlebase_insertions_Trinuc_nonhomopolymer.bed', 'w+')
#h5 = open('rad16_sorted_insertions_Pentanuc.bed', 'w+')
h6 = open('rad16_sorted_nonhomopolymer_singleinsertions_25nuc.bed', 'w+')

j1 = open("rad26_insertionsv3_sorted_tandem.bed")
j2 = open('rad26_singlebase_insertions_Trinuc.bed', 'w+')
#j3 = open('rad26_singlebase_combined_Trinuc_dipyrimidine.bed', 'w+')
#j4 = open('rad26_singlebase_combined_Trinuc_nondipyrimidine.bed', 'w+')
j3 = open('rad26_singlebase_insertions_Trinuc_homopolymer.bed', 'w+')
j4 = open('rad26_singlebase_insertions_Trinuc_nonhomopolymer.bed', 'w+')
#j5 = open('rad26_sorted_insertions_Pentanuc.bed', 'w+')
j6 = open('rad26_sorted_nonhomopolymer_singleinsertions_25nuc.bed', 'w+')

k1 = open("rad30_insertions_sorted_tandem.bed")
k2 = open('rad30_singlebase_insertions_Trinuc.bed', 'w+')
#k3 = open('rad30_singlebase_combined_Trinuc_dipyrimidine.bed', 'w+')
#k4 = open('rad30_singlebase_combined_Trinuc_nondipyrimidine.bed', 'w+')
k3 = open('rad30_singlebase_insertions_Trinuc_homopolymer.bed', 'w+')
k4 = open('rad30_singlebase_insertions_Trinuc_nonhomopolymer.bed', 'w+')
#k5 = open('rad30_sorted_insertions_Pentanuc.bed', 'w+')
k6 = open('rad30_sorted_nonhomopolymer_singleinsertions_25nuc.bed', 'w+')




homopolymer_data = MakeMultiNucFile(g1, g6)
polymer_counts = GetPolymerCounts(homopolymer_data[0],  homopolymer_data[2], homopolymer_data[3], homopolymer_data[4], homopolymer_data[5], homopolymer_data[6], homopolymer_data[7], homopolymer_data[8])[0]
homopolymer = GetPolymerCounts(homopolymer_data[0],  homopolymer_data[2], homopolymer_data[3], homopolymer_data[4], homopolymer_data[5], homopolymer_data[6], homopolymer_data[7], homopolymer_data[8])[1]
homopolymer_data16 = MakeMultiNucFile(h1, h6)
polymer_counts16 = GetPolymerCounts(homopolymer_data16[0],  homopolymer_data16[2], homopolymer_data16[3], homopolymer_data16[4], homopolymer_data16[5], homopolymer_data16[6], homopolymer_data16[7], homopolymer_data16[8])[0]
homopolymer16 = GetPolymerCounts(homopolymer_data16[0],  homopolymer_data16[2], homopolymer_data16[3], homopolymer_data16[4], homopolymer_data16[5], homopolymer_data16[6], homopolymer_data16[7], homopolymer_data16[8])[1]
homopolymer_data26 = MakeMultiNucFile(j1, j6)
polymer_counts26 = GetPolymerCounts(homopolymer_data26[0],  homopolymer_data26[2], homopolymer_data26[3], homopolymer_data26[4], homopolymer_data26[5], homopolymer_data26[6], homopolymer_data26[7], homopolymer_data26[8])[0]
homopolymer26 = GetPolymerCounts(homopolymer_data26[0],  homopolymer_data26[2], homopolymer_data26[3], homopolymer_data26[4], homopolymer_data26[5], homopolymer_data26[6], homopolymer_data26[7], homopolymer_data26[8])[1]
homopolymer_data30 = MakeMultiNucFile(k1, k6)
polymer_counts30 = GetPolymerCounts(homopolymer_data30[0],  homopolymer_data30[2], homopolymer_data30[3], homopolymer_data30[4], homopolymer_data30[5], homopolymer_data30[6], homopolymer_data30[7], homopolymer_data30[8])[0]
homopolymer30 = GetPolymerCounts(homopolymer_data30[0],  homopolymer_data30[2], homopolymer_data30[3], homopolymer_data30[4], homopolymer_data30[5], homopolymer_data30[6], homopolymer_data30[7], homopolymer_data30[8])[1]
homopolymer_data2 = MakeMultiNucFile(a1, a6)
polymer_counts2 = GetPolymerCounts(homopolymer_data2[0],  homopolymer_data2[2], homopolymer_data2[3], homopolymer_data2[4], homopolymer_data2[5], homopolymer_data2[6], homopolymer_data2[7], homopolymer_data2[8])[0]
homopolymer_data2_16 = MakeMultiNucFile(b1, b6)
polymer_counts2_16 = GetPolymerCounts(homopolymer_data2_16[0],  homopolymer_data2_16[2], homopolymer_data2_16[3], homopolymer_data2_16[4], homopolymer_data2_16[5], homopolymer_data2_16[6], homopolymer_data2_16[7], homopolymer_data2_16[8])[0]
homopolymer_data2_30 = MakeMultiNucFile(e1, e6)
polymer_counts2_30 = GetPolymerCounts(homopolymer_data2_30[0],  homopolymer_data2_30[2], homopolymer_data2_30[3], homopolymer_data2_30[4], homopolymer_data2_30[5], homopolymer_data2_30[6], homopolymer_data2_30[7], homopolymer_data2_30[8])[0]
homopolymer_data2_26 = MakeMultiNucFile(d1, d6)
polymer_counts2_26 = GetPolymerCounts(homopolymer_data2_26[0],  homopolymer_data2_26[2], homopolymer_data2_26[3], homopolymer_data2_26[4], homopolymer_data2_26[5], homopolymer_data2_26[6], homopolymer_data2_26[7], homopolymer_data2_26[8])[0]

homopolymer_combined_data = MakeMultiNucFile(a7, a8)
homopolymer_combined = GetPolymerCounts(homopolymer_combined_data[0],  homopolymer_combined_data[2], homopolymer_combined_data[3], homopolymer_combined_data[4], homopolymer_combined_data[5], homopolymer_combined_data[6], homopolymer_combined_data[7], homopolymer_combined_data[8])[1]
homopolymer_combined_data16 = MakeMultiNucFile(b7, b8)
homopolymer16_combined = GetPolymerCounts(homopolymer_combined_data16[0],  homopolymer_combined_data16[2], homopolymer_combined_data16[3], homopolymer_combined_data16[4], homopolymer_combined_data16[5], homopolymer_combined_data16[6], homopolymer_combined_data16[7], homopolymer_combined_data16[8])[1]
homopolymer_combined_data26 = MakeMultiNucFile(d7, d8)
homopolymer26_combined = GetPolymerCounts(homopolymer_combined_data26[0],  homopolymer_combined_data26[2], homopolymer_combined_data26[3], homopolymer_combined_data26[4], homopolymer_combined_data26[5], homopolymer_combined_data26[6], homopolymer_combined_data26[7], homopolymer_combined_data26[8])[1]


WT_data = homopolymer_combined_data
WT_isolate = WT_data[4]
WT_chromosome = WT_data[2]
WT_position1 = WT_data[9]
WT_position = WT_data[3]
WT_trinucleotide = WT_data[5]
WT_mutation = WT_data[6]
WT_strand = WT_data[1]
WT_multinucleotide = WT_data[0] 
WT_lengths = WT_data[8]
WT_allele = WT_data[7]
WT_lengths2 = homopolymer_data2[8]
WT_lengths3 = homopolymer_data[8]


rad16_data = homopolymer_combined_data16
rad16_isolate = rad16_data[4]
rad16_chromosome = rad16_data[2]
rad16_position1 = rad16_data[9]
rad16_position = rad16_data[3]
rad16_trinucleotide = rad16_data[5]
rad16_mutation = rad16_data[6]
rad16_strand = rad16_data[1]
rad16_multinucleotide = rad16_data[0] 
rad16_lengths = rad16_data[8]
rad16_allele = rad16_data[7]
rad16_lengths2 = homopolymer_data2_16[8]
rad16_lengths3 = homopolymer_data16[8]

rad26_data = homopolymer_combined_data26

rad26_isolate = rad26_data[4]
rad26_chromosome = rad26_data[2]
rad26_position1 = rad26_data[9]
rad26_position = rad26_data[3]
rad26_trinucleotide = rad26_data[5]
rad26_mutation = rad26_data[6]
rad26_strand = rad26_data[1]
rad26_multinucleotide = rad26_data[0] 
rad26_lengths = rad26_data[8]
rad26_allele = rad26_data[7]
rad26_lengths2 = homopolymer_data2_26[8]
rad26_lengths3 = homopolymer_data26[8]

rad30_data = homopolymer_data30
rad30_isolate = rad30_data[4]
rad30_chromosome = rad30_data[2]
rad30_position1 = rad30_data[9]
rad30_position = rad30_data[3]
rad30_trinucleotide = rad30_data[5]
rad30_mutation = rad30_data[6]
rad30_strand = rad30_data[1]
rad30_multinucleotide = rad30_data[0] 
rad30_lengths = rad30_data[8]
rad30_lengths2 = homopolymer_data2_30[8]




MakeSequenceFile(sequencefile1, sequencefile2, WT_multinucleotide, WT_mutation, WT_allele, WT_strand)
MakeSequenceFile(sequencefile3, sequencefile4, rad16_multinucleotide, rad16_mutation, rad16_allele, rad16_strand)
MakeSequenceFile(sequencefile5, sequencefile6, rad26_multinucleotide, rad26_mutation, rad26_allele, rad26_strand)

sequencefile1.close()
sequencefile2.close()
sequencefile3.close()
sequencefile4.close()
sequencefile5.close()
sequencefile6.close()

MakeTriNucFile(a1, a2, a3, a4,WT_lengths2 )
MakeTriNucFile(b1, b2, b3, b4, rad16_lengths2)   
MakeTriNucFile(d1, d2, d3, d4, rad26_lengths2)
MakeTriNucFile(e1, e2, e3, e4, rad30_lengths2) 
MakeInsertionFile(g1, g2, g3, g4, WT_lengths3)
MakeInsertionFile(h1, h2, h3, h4, rad16_lengths3)   
MakeInsertionFile(j1, j2, j3, j4, rad26_lengths3)
MakeInsertionFile(k1, k2, k3, k4, rad30_lengths) 
#MakeTriNucFile(ca, cb, cc, cd)
#MakeTriNucFile(c2a, c2b, c2c, c2d)
#MakeTriNucFile(c3a, c3b, c3c, c3d)

g1.close()
h1.close()
j1.close()
k1.close()
g2.close()
g3.close()
g4.close()
h2.close()
h3.close()
h4.close()
j2.close()
j3.close()
j4.close()
k2.close()
k3.close()
k4.close()
a1.close()
b1.close()
e1.close()

a1.close()
a2.close()
a3.close()
a4.close()
a5.close()
a6.close()
b1.close()
b2.close()
b3.close()
b4.close()
b5.close()
b6.close()
d1.close()
d2.close()
d3.close()
d4.close()
d5.close()
d6.close()
e1.close()
e2.close()
e3.close()
e4.close()
e5.close()
e6.close()
g1.close()
g6.close()
h1.close()
h6.close()
j1.close()
j6.close()
k1.close()
k6.close()

WT_isolate_list = {}
for i1 in range(0, len(WT_isolate)):

    if WT_isolate[i1] not in WT_isolate_list.keys():
            WT_isolate_list[WT_isolate[i1]] = [((WT_position[i1]),  str(WT_chromosome[i1]), WT_trinucleotide[i1], WT_mutation[i1], WT_strand[i1], WT_multinucleotide[i1], WT_isolate[i1], int(WT_position1[i1]))]
    else:
        WT_isolate_list[WT_isolate[i1]].append(((WT_position[i1]),  str(WT_chromosome[i1]), WT_trinucleotide[i1], WT_mutation[i1], WT_strand[i1], WT_multinucleotide[i1], WT_isolate[i1], int(WT_position1[i1])))
rad16_isolate_list = {}
for i1 in range(0, len(rad16_isolate)):

    if rad16_isolate[i1] not in rad16_isolate_list.keys():
            rad16_isolate_list[rad16_isolate[i1]] = [(rad16_position[i1],  str(rad16_chromosome[i1]), rad16_trinucleotide[i1], rad16_mutation[i1], rad16_strand[i1], rad16_multinucleotide[i1], rad16_isolate[i1], rad16_position1[i1])]
    else:
        rad16_isolate_list[rad16_isolate[i1]] .append((rad16_position[i1],  str(rad16_chromosome[i1]), rad16_trinucleotide[i1], rad16_mutation[i1], rad16_strand[i1], rad16_multinucleotide[i1], rad16_isolate[i1], rad16_position1[i1]))

rad26_isolate_list = {}
for i1 in range(0, len(rad26_isolate)):

    if rad26_isolate[i1] not in rad26_isolate_list.keys():
            rad26_isolate_list[rad26_isolate[i1]] = [(rad26_position[i1],  str(rad26_chromosome[i1]), rad26_trinucleotide[i1], rad26_mutation[i1], rad26_strand[i1], rad26_multinucleotide[i1], rad26_isolate[i1], rad26_position1[i1])]
    else:
        rad26_isolate_list[rad26_isolate[i1]] .append((rad26_position[i1],  str(rad26_chromosome[i1]), rad26_trinucleotide[i1], rad26_mutation[i1], rad26_strand[i1], rad26_multinucleotide[i1], rad26_isolate[i1], rad26_position1[i1]))

rad30_isolate_list = {}
for i1 in range(0, len(rad30_isolate)):

    if rad30_isolate[i1] not in rad30_isolate_list.keys():
            rad30_isolate_list[rad30_isolate[i1]] = [(rad30_position[i1],  str(rad30_chromosome[i1]), rad30_trinucleotide[i1], rad30_mutation[i1], rad30_strand[i1], rad30_multinucleotide[i1], rad30_isolate[i1], rad30_position1[i1])]
    else:
        rad30_isolate_list[rad30_isolate[i1]] .append((rad30_position[i1],  str(rad30_chromosome[i1]), rad30_trinucleotide[i1], rad30_mutation[i1], rad30_strand[i1], rad30_multinucleotide[i1], rad30_isolate[i1], rad30_position1[i1]))






def getStrforInt(string):
    new_str = ""
    for s in string:
        if s.isnumeric():
            new_str = new_str + s
        else:
            return new_str
    return new_str



with open("UVB_filter_mutations.txt") as f2:  #Input UVB whole-genome sequencing data
    data2 = f2.read()
lines2 = [s.strip().split() for s in data2.splitlines()]
lines2.pop(0)
genotypes2 = {'WT': 0, 'Rad16':0, 'Rad26':0, 'Rad30': 0}
WT_indels2 = 0
rad16_indels2 = 0
rad26_indels2 = 0
rad30_indels2 = 0
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
WT_double2 = 0
rad16_double2 = 0
rad26_double2 = 0
rad30_double2 = 0
WT_subs2 = 0
rad16_subs2 = 0
rad26_subs2 = 0
rad30_subs2 = 0
for line in lines2:
    
    if 'Mito' not in line[5] and line[2] != 'RP' and line[0] != 'B01' and line[0] != 'C13':

        isolates2.append(line[0])
        genotype2.append(line[1])
        mutation_type_2.append(line[13])
        
        if line[2] != 'RP':
            chromosome2.append(line[4])
        else:
            chromosome2.append(line[5])
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
        
        mutation2.append(line[15])
        
#Count mutations of different types
for g in range(0, len(genotype2)):
    if 'WT' in genotype2[g]:
        genotypes2['WT'] = genotypes2['WT'] + 1
    if 'Rad16' in genotype2[g]:
            genotypes2['Rad16'] = genotypes2['Rad16'] + 1
    if 'Rad26' in genotype2[g]:
            genotypes2['Rad26'] = genotypes2['Rad26'] + 1
    if 'Rad30' in genotype2[g]:
            genotypes2['Rad30'] = genotypes2['Rad30'] + 1
    if 'WT' in genotype2[g]:
    
        if mutation_type_2[g] == "Deletion" or mutation_type_2[g] == "Insertion" or mutation_type_2[g] == "DIP":
            WT_indels2 = WT_indels2 + 1
        if mutation_type_2[g] == "MNV" or mutation_type_2[g] == "DNP":
            WT_double2 = WT_double2 + 1
        if mutation_type_2[g] == "SNP" or mutation_type_2[g] == "SNV":
            WT_subs2 = WT_subs2 + 1
           
    if 'Rad16' in genotype2[g]:
            if mutation_type_2[g] == "DIP" or  mutation_type_2[g] == "Deletion" or mutation_type_2[g] == "Insertion":
                rad16_indels2 = rad16_indels2 + 1
            if mutation_type_2[g] == "DNP" or mutation_type_2[g] == "MNV":
                rad16_double2 = rad16_double2 + 1
            if mutation_type_2[g] == "SNP" or mutation_type_2[g] == "SNV":
                rad16_subs2 = rad16_subs2 + 1
    if 'Rad26' in genotype2[g]:
            if mutation_type_2[g] == "DIP" or mutation_type_2[g] == "Deletion" or mutation_type_2[g] == "Insertion":
                rad26_indels2 = rad26_indels2 + 1
            if mutation_type_2[g] == "MNV" or mutation_type_2[g] == "DNP":
                rad26_double2 = rad26_double2 + 1
            if mutation_type_2[g] == "SNV" or mutation_type_2[g] == "SNP":
                rad26_subs2 = rad26_subs2 + 1
    if 'Rad30' in genotype2[g]:
            if mutation_type_2[g] == "DIP" or mutation_type_2[g] == "Deletion" or mutation_type_2[g] == "Insertion":
                rad30_indels2 = rad30_indels2 + 1
            if mutation_type_2[g] == "SNV" or mutation_type_2[g] == "SNP":
                rad30_subs2 = rad30_subs2 + 1
            if mutation_type_2[g] == "MNV" or mutation_type_2[g] == "DNP":
                rad30_double2 = rad30_double2 + 1

#Input trinucleotide sequence context
for pos in range(0, len(position2_2)):
    
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
            isolate_list2[(isolates2[i1], genotype2[i1])] = [(position2_2[i1],  str(chromosome2[i1]), trinucleotide2[i1], mutation2[i1], strand2[i1], mutation_type_2[i1], isolates2[i1], allele2[i1])]
    else:
        isolate_list2[(isolates2[i1], genotype2[i1])] .append((position2_2[i1],  str(chromosome2[i1]), trinucleotide2[i1], mutation2[i1], strand2[i1], mutation_type_2[i1], isolates2[i1], allele2[i1]))


file4 = open('UVB_WT_mutationsperisolate', 'w+')
file5 = open('UVB_rad16_mutationsperisolate', 'w+')
file6 = open('UVB_rad26_mutationsperisolate', 'w+')
MutationsPerIsolate(isolate_list2, file4, file5, file6)

file1.close()
file2.close()
file3.close()
file4.close()
file5.close()
file6.close()

WT_insertions2 = {}
WT_deletions2 = {}
rad16_insertions2 = {}
rad16_deletions2 = {}
rad26_insertions2 = {}
rad26_deletions2 = {}
rad30_insertions2 = {}
rad30_deletions2 = {}


#Count the deletions and insertions of each length

for m in range(0, len(mutation_type_2)):
    if genotype2[m] == 'WT':
        
        if mutation_type_2[m] == 'Insertion':
            if len(mutation2[m]) in WT_insertions2.keys():
                    WT_insertions2[len(mutation2[m])] = WT_insertions2[len(mutation2[m])] + 1
            else: 
                    WT_insertions2[len(mutation2[m])] = 1
        elif mutation_type_2[m] == 'Deletion':
       
            if len(allele2[m]) in WT_deletions2.keys():
                WT_deletions2[len(allele2[m])] = WT_deletions2[len(allele2[m])] + 1
            else: 
                WT_deletions2[len(allele2[m])] = 1
    if 'Rad16' in genotype2[m]:
        if mutation_type_2[m] == 'Insertion':
      
            if len(mutation2[m]) in rad16_insertions2.keys():
                    rad16_insertions2[len(mutation2[m])] = rad16_insertions2[len(mutation2[m])] + 1
            else: 
                    rad16_insertions2[len(mutation2[m])] = 1
        elif mutation_type_2[m] == 'Deletion':
        #elif '-' in mutation2[m]:
            if len(allele2[m]) in rad16_deletions2.keys():
                rad16_deletions2[len(allele2[m])] = rad16_deletions2[len(allele2[m])] + 1
            else: 
                rad16_deletions2[len(allele2[m])] = 1
    if 'Rad26' in genotype2[m]:
        if mutation_type_2[m] == 'Insertion':
        #if '-' in allele2[m]:
            if len(mutation2[m]) in rad26_insertions2.keys():
                    rad26_insertions2[len(mutation2[m])] = rad26_insertions2[len(mutation2[m])] + 1
            else: 
                    rad26_insertions2[len(mutation2[m])] = 1
        elif mutation_type_2[m] == 'Deletion':
        #elif '-' in mutation2[m]:
            if len(allele2[m]) in rad26_deletions2.keys():
                rad26_deletions2[len(allele2[m])] = rad26_deletions2[len(allele2[m])] + 1
            else: 
                rad26_deletions2[len(allele2[m])] = 1
    if 'Rad30' in genotype2[m]:
        #if mutation_type_2[m] == 'Insertion':
        if '-' in allele2[m]:
            if len(mutation2[m]) in rad30_insertions2.keys():
                    rad30_insertions2[len(mutation2[m])] = rad30_insertions2[len(mutation2[m])] + 1
            else: 
                    rad30_insertions2[len(mutation2[m])] = 1
        #elif mutation_type_2[m] == 'Deletion':
        elif '-' in mutation2[m]:
            if len(allele2[m]) in rad30_deletions2.keys():
                rad30_deletions2[len(allele2[m])] = rad30_deletions2[len(allele2[m])] + 1
            else: 
                rad30_deletions2[len(allele2[m])] = 1

g1a = open("UVB_WT_insertions_sorted_tandem.bed", 'w+')
g1b = open("UVB_rad16_insertions_sorted_tandem.bed", 'w+')
g1c = open("UVB_rad26_insertions_sorted_tandem.bed", 'w+')
g1d = open("UVB_rad30_insertions_sorted_tandem.bed", 'w+')
g2a = open("UVB_WT_deletions_sorted_tandem.bed", 'w+')
g2b = open("UVB_rad16_deletions_sorted_tandem.bed", 'w+')
g2c = open("UVB_rad26_deletions_sorted_tandem.bed", 'w+')
g2d = open("UVB_rad30_deletions_sorted_tandem.bed", 'w+')
g3a = open("UVB_WT_all_sorted_tandem.bed", 'w+')
g3b = open("UVB_rad16_all_sorted_tandem.bed", 'w+')
g3c = open("UVB_rad26_all_sorted_tandem.bed", 'w+')
g3d = open("UVB_rad30_all_sorted_tandem.bed", 'w+')

g4a = open("UVB_WT_allcombined_sorted_tandem.bed", 'w+')
g4b = open("UVB_rad16_allcombined_sorted_tandem.bed", 'w+')
g4c = open("UVB_rad26_allcombined_sorted_tandem.bed", 'w+')
g4d = open("UVB_rad30_allcombined_sorted_tandem.bed", 'w+')



WriteChr('chrI',g1a, g1b, g1c, g1d, g2a, g2b, g2c, g2d, g3a, g3b, g3c, g3d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrII', g1a, g1b, g1c, g1d, g2a, g2b, g2c, g2d, g3a, g3b, g3c, g3d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrIII',  g1a, g1b, g1c, g1d, g2a, g2b, g2c, g2d, g3a, g3b, g3c, g3d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrIV', g1a, g1b, g1c, g1d, g2a, g2b, g2c, g2d, g3a, g3b, g3c, g3d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrIX', g1a, g1b, g1c, g1d, g2a, g2b, g2c, g2d, g3a, g3b, g3c, g3d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrV', g1a, g1b, g1c, g1d, g2a, g2b, g2c, g2d, g3a, g3b, g3c, g3d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrVI', g1a, g1b, g1c, g1d, g2a, g2b, g2c, g2d, g3a, g3b, g3c, g3d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrVII', g1a, g1b, g1c, g1d, g2a, g2b, g2c, g2d, g3a, g3b, g3c, g3d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrVIII', g1a, g1b, g1c, g1d, g2a, g2b, g2c, g2d, g3a, g3b, g3c, g3d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrX', g1a, g1b, g1c, g1d, g2a, g2b, g2c, g2d, g3a, g3b, g3c, g3d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrXI', g1a, g1b, g1c, g1d, g2a, g2b, g2c, g2d, g3a, g3b, g3c, g3d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrXII', g1a, g1b, g1c, g1d, g2a, g2b, g2c, g2d, g3a, g3b, g3c, g3d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrXIII', g1a, g1b, g1c, g1d, g2a, g2b, g2c, g2d, g3a, g3b, g3c, g3d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrXIV', g1a, g1b, g1c, g1d, g2a, g2b, g2c, g2d, g3a, g3b, g3c, g3d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2) 
WriteChr('chrXV', g1a, g1b, g1c, g1d, g2a, g2b, g2c, g2d, g3a, g3b, g3c, g3d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteChr('chrXVI', g1a, g1b, g1c, g1d, g2a, g2b, g2c, g2d, g3a, g3b, g3c, g3d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)

WriteAllChr('chrI',g4a, g4b, g4c, g4d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteAllChr('chrII', g4a, g4b, g4c, g4d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteAllChr('chrIII',  g4a, g4b, g4c, g4d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteAllChr('chrIV', g4a, g4b, g4c, g4d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteAllChr('chrIX', g4a, g4b, g4c, g4d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteAllChr('chrV', g4a, g4b, g4c, g4d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteAllChr('chrVI', g4a, g4b, g4c, g4d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteAllChr('chrVII', g4a, g4b, g4c, g4d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteAllChr('chrVIII', g4a, g4b, g4c, g4d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteAllChr('chrX', g4a, g4b, g4c, g4d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteAllChr('chrXI', g4a, g4b, g4c, g4d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteAllChr('chrXII', g4a, g4b, g4c, g4d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteAllChr('chrXIII', g4a, g4b, g4c, g4d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteAllChr('chrXIV', g4a, g4b, g4c, g4d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2) 
WriteAllChr('chrXV', g4a, g4b, g4c, g4d,chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)
WriteAllChr('chrXVI', g4a, g4b, g4c, g4d, chromosome2, position1_2, position2_2,
trinucleotide2, allele2, mutation2, genotype2, mutation_type_2, isolates2, strand2)



g1a.close()
g1b.close()
g1c.close()
g1d.close()
g2a.close()
g2b.close()
g2c.close()
g2d.close()
g3a.close()
g3b.close()
g3c.close()
g3d.close()
g4a.close()
g4b.close()
g4c.close()
g4d.close()




a1b = open("UVB_WT_deletions_sorted_tandem.bed")
a2b = open('UVB_WT_singlebase_deletions_Trinuc.bed', 'w+')
#a3b = open('UVB_WT_singlebase_deletions_Trinuc_dipyrimidine.bed', 'w+')
#a4b = open('UVB_WT_singlebase_deletions_Trinuc_nondipyrimidine.bed', 'w+')
a3b = open('UVB_WT_singlebase_deletions_Trinuc_homopolymer.bed', 'w+')
a4b = open('UVB_WT_singlebase_deletions_Trinuc_nonhomopolymer.bed', 'w+')
a5b = open('UVB_WT_sorted_deletions_Pentanuc.bed', 'w+')
a6b = open('UVB_WT_sorted_singledeletions_25nuc.bed', 'w+')
a7b = open("UVB_WT_allcombined_sorted_tandem.bed")
a8b = open('UVB_WT_sorted_all_25nuc.bed', 'w+')

b1b = open("UVB_rad16_deletions_sorted_tandem.bed")
b2b = open('UVB_rad16_singlebase_deletions_Trinuc.bed', 'w+')
#b3b = open('UVB_rad16_singlebase_deletions_Trinuc_dipyrimidine.bed', 'w+')
#b4b = open('UVB_rad16_singlebase_deletions_Trinuc_nondipyrimidine.bed', 'w+')
b3b = open('UVB_rad16_singlebase_deletions_Trinuc_homopolymer.bed', 'w+')
b4b = open('UVB_rad16_singlebase_deletions_Trinuc_nonhomopolymer.bed', 'w+')
b5b = open('UVB_rad16_sorted_deletions_Pentanuc.bed', 'w+')
b6b = open('UVB_rad16_sorted_singledeletions_25nuc.bed', 'w+')
b7b = open("UVB_rad16_allcombined_sorted_tandem.bed")
b8b = open('UVB_rad16_sorted_all_25nuc.bed', 'w+')

d1b = open("UVB_rad26_deletions_sorted_tandem.bed")
d2b = open('UVB_rad26_singlebase_deletions_Trinuc.bed', 'w+')
#d3b = open('UVB_rad26_singlebase_deletions_Trinuc_dipyrimidine.bed', 'w+')
#d4b = open('UVB_rad26_singlebase_deletions_Trinuc_nondipyrimidine.bed', 'w+')
d3b = open('UVB_rad26_singlebase_deletions_Trinuc_homopolymer.bed', 'w+')
d4b = open('UVB_rad26_singlebase_deletions_Trinuc_nonhomopolymer.bed', 'w+')
d5b = open('UVB_rad26_sorted_deletions_Pentanuc.bed', 'w+')
d6b = open('UVB_rad26_sorted_singledeletions_25nuc.bed', 'w+')
d7b = open("UVB_rad26_allcombined_sorted_tandem.bed")
d8b = open('UVB_rad26_sorted_all_25nuc.bed', 'w+')

e1b = open("UVB_rad30_deletions_sorted_tandem.bed")
e2b = open('UVB_rad30_singlebase_deletions_Trinuc.bed', 'w+')
#e3b = open('UVB_rad30_singlebase_deletions_Trinuc_dipyrimidine.bed', 'w+')
#e4b = open('UVB_rad30_singlebase_deletions_Trinuc_nondipyrimidine.bed', 'w+')
e3b = open('UVB_rad30_singlebase_deletions_Trinuc_homopolymer.bed', 'w+')
e4b = open('UVB_rad30_singlebase_deletions_Trinuc_nonhomopolymer.bed', 'w+')
e5b = open('UVB_rad30_sorted_deletions_Pentanuc.bed', 'w+')
e6b = open('UVB_rad30_sorted_singledeletions_25nuc.bed', 'w+')

m1 = open("UVB_WT_insertions_sorted_tandem.bed")
m2 = open('UVB_WT_singlebase_insertions_Trinuc.bed', 'w+')
#m3 = open('UVB_WT_singlebase_insertions_Trinuc_dipyrimidine.bed', 'w+')
#m4 = open('UVB_WT_singlebase_insertions_Trinuc_nondipyrimidine.bed', 'w+')
m3 = open('UVB_WT_singlebase_insertions_Trinuc_homopolymer.bed', 'w+')
m4 = open('UVB_WT_singlebase_insertions_Trinuc_nonhomopolymer.bed', 'w+')
#g5 = open('WT_sorted_insertions_Pentanuc.bed', 'w+')
m6 = open('UVB_WT_sorted_nonhomopolymer_singleinsertions_25nuc.bed', 'w+')

h1b = open("UVB_rad16_insertions_sorted_tandem.bed")
h2b = open('UVB_rad16_singlebase_insertions_Trinuc.bed', 'w+')
#h3b = open('UVB_rad16_singlebase_combined_Trinuc_dipyrimidine.bed', 'w+')
#h4b = open('UVB_rad16_singlebase_combined_Trinuc_nondipyrimidine.bed', 'w+')
h3b = open('UVB_rad16_singlebase_insertions_Trinuc_homopolymer.bed', 'w+')
h4b = open('UVB_rad16_singlebase_insertions_Trinuc_nonhomopolymer.bed', 'w+')
#h5 = open('rad16_sorted_insertions_Pentanuc.bed', 'w+')
h6b = open('UVB_rad16_sorted_nonhomopolymer_singleinsertions_25nuc.bed', 'w+')

j1b = open("UVB_rad26_insertions_sorted_tandem.bed")
j2b = open('UVB_rad26_singlebase_insertions_Trinuc.bed', 'w+')
#j3b = open('UVB_rad26_singlebase_combined_Trinuc_dipyrimidine.bed', 'w+')
#j4b = open('UVB_rad26_singlebase_combined_Trinuc_nondipyrimidine.bed', 'w+')
j3b = open('UVB_rad26_singlebase_insertions_Trinuc_homopolymer.bed', 'w+')
j4b = open('UVB_rad26_singlebase_insertions_Trinuc_nonhomopolymer.bed', 'w+')
#j5 = open('rad26_sorted_insertions_Pentanuc.bed', 'w+')
j6b = open('UVB_rad26_sorted_nonhomopolymer_singleinsertions_25nuc.bed', 'w+')
k1b = open("UVB_rad30_insertions_sorted_tandem.bed")
k2b = open('UVB_rad30_singlebase_insertions_Trinuc.bed', 'w+')
#k3b = open('UVB_rad30_singlebase_combined_Trinuc_dipyrimidine.bed', 'w+')
#k4b = open('UVB_rad30_singlebase_combined_Trinuc_nondipyrimidine.bed', 'w+')
k3b = open('UVB_rad30_singlebase_insertions_Trinuc_homopolymer.bed', 'w+')
k4b = open('UVB_rad30_singlebase_insertions_Trinuc_nonhomopolymer.bed', 'w+')
#k5 = open('rad30_sorted_insertions_Pentanuc.bed', 'w+')
k6b = open('UVB_rad30_sorted_singleinsertions_25nuc.bed', 'w+')


sequencefile7 = open('UVB_single_WT_Tdeletion_sequences', 'w+')
sequencefile8 = open('UVB_single_WT_Tinsertion_sequences', 'w+')

sequencefile9 = open('UVB_single_rad16_Tdeletion_sequences', 'w+')
sequencefile10 = open('UVB_single_rad16_Tinsertion_sequences', 'w+')

sequencefile11 = open('UVB_single_rad26_Tdeletion_sequences', 'w+')
sequencefile12 = open('UVB_single_rad26_Tinsertion_sequences', 'w+')

homopolymer_UVB_data2 = MakeMultiNucFile(a1b, a6b)
WT_del_lengths = homopolymer_UVB_data2[8]
homopolymer_UVB_data16_2 = MakeMultiNucFile(b1b, b6b)
rad16_del_lengths = homopolymer_UVB_data16_2[8]
homopolymer_UVB_data26_2 = MakeMultiNucFile(d1b, d6b)
rad26_del_lengths = homopolymer_UVB_data26_2[8]
rad30_del_lengths = MakeMultiNucFile(e1b, e6b)[8]
homopolymer_UVB_data = MakeMultiNucFile(m1, m6)
WT_ins_lengths = homopolymer_UVB_data[8]
polymer_counts_UVB = GetPolymerCounts(homopolymer_UVB_data[0],  homopolymer_UVB_data[2], homopolymer_UVB_data[3], homopolymer_UVB_data[4], homopolymer_UVB_data[5], homopolymer_UVB_data[6], homopolymer_UVB_data[7], homopolymer_UVB_data[8])[0]
homopolymerB = GetPolymerCounts(homopolymer_UVB_data[0],  homopolymer_UVB_data[2], homopolymer_UVB_data[3], homopolymer_UVB_data[4], homopolymer_UVB_data[5], homopolymer_UVB_data[6], homopolymer_UVB_data[7], homopolymer_UVB_data[8])[1]
homopolymer_UVB_data16 = MakeMultiNucFile(h1b, h6b)
rad16_ins_lengths = homopolymer_UVB_data16[8]
polymer_counts_UVB16 = GetPolymerCounts(homopolymer_UVB_data16[0],  homopolymer_UVB_data16[2], homopolymer_UVB_data16[3], homopolymer_UVB_data16[4], homopolymer_UVB_data16[5], homopolymer_UVB_data16[6], homopolymer_UVB_data16[7], homopolymer_UVB_data16[8])[0]
homopolymer16B = GetPolymerCounts(homopolymer_UVB_data16[0],  homopolymer_UVB_data16[2], homopolymer_UVB_data16[3], homopolymer_UVB_data16[4], homopolymer_UVB_data16[5], homopolymer_UVB_data16[6], homopolymer_UVB_data16[7], homopolymer_UVB_data16[8])[1]
homopolymer_UVB_data30 = MakeMultiNucFile(k1b, k6b)
polymer_counts_UVB30 = GetPolymerCounts(homopolymer_UVB_data30[0],  homopolymer_UVB_data30[2], homopolymer_UVB_data30[3], homopolymer_UVB_data30[4], homopolymer_UVB_data30[5], homopolymer_UVB_data30[6], homopolymer_UVB_data30[7], homopolymer_UVB_data30[8])[0]
homopolymer30B = GetPolymerCounts(homopolymer_UVB_data30[0],  homopolymer_UVB_data30[2], homopolymer_UVB_data30[3], homopolymer_UVB_data30[4], homopolymer_UVB_data30[5], homopolymer_UVB_data30[6], homopolymer_UVB_data30[7], homopolymer_UVB_data30[8])[1]
homopolymer_UVB_data26 = MakeMultiNucFile(j1b, j6b)
rad26_ins_lengths = homopolymer_UVB_data26[8]
polymer_counts_UVB26 = GetPolymerCounts(homopolymer_UVB_data26[0],  homopolymer_UVB_data26[2], homopolymer_UVB_data26[3], homopolymer_UVB_data26[4], homopolymer_UVB_data26[5], homopolymer_UVB_data26[6], homopolymer_UVB_data26[7], homopolymer_UVB_data26[8])[0]
homopolymer26B = GetPolymerCounts(homopolymer_UVB_data26[0],  homopolymer_UVB_data26[2], homopolymer_UVB_data26[3], homopolymer_UVB_data26[4], homopolymer_UVB_data26[5], homopolymer_UVB_data26[6], homopolymer_UVB_data26[7], homopolymer_UVB_data26[8])[1]


polymer_counts_UVB2 = GetPolymerCounts(homopolymer_UVB_data2[0],  homopolymer_UVB_data2[2], homopolymer_UVB_data2[3], homopolymer_UVB_data2[4], homopolymer_UVB_data2[5], homopolymer_UVB_data2[6], homopolymer_UVB_data2[7], homopolymer_UVB_data2[8])[0]
homopolymerB2 = GetPolymerCounts(homopolymer_UVB_data2[0],  homopolymer_UVB_data2[2], homopolymer_UVB_data2[3], homopolymer_UVB_data2[4], homopolymer_UVB_data2[5], homopolymer_UVB_data2[6], homopolymer_UVB_data2[7], homopolymer_UVB_data2[8])[1]
polymer_counts_UVB16_2 = GetPolymerCounts(homopolymer_UVB_data16_2[0],  homopolymer_UVB_data16_2[2], homopolymer_UVB_data16_2[3], homopolymer_UVB_data16_2[4], homopolymer_UVB_data16_2[5], homopolymer_UVB_data16_2[6], homopolymer_UVB_data16_2[7], homopolymer_UVB_data16_2[8])[0]
homopolymer16B2 = GetPolymerCounts(homopolymer_UVB_data16_2[0],  homopolymer_UVB_data16_2[2], homopolymer_UVB_data16_2[3], homopolymer_UVB_data16_2[4], homopolymer_UVB_data16_2[5], homopolymer_UVB_data16_2[6], homopolymer_UVB_data16_2[7], homopolymer_UVB_data16_2[8])[1]
polymer_counts_UVB26_2 = GetPolymerCounts(homopolymer_UVB_data26_2[0],  homopolymer_UVB_data26_2[2], homopolymer_UVB_data26_2[3], homopolymer_UVB_data26_2[4], homopolymer_UVB_data26_2[5], homopolymer_UVB_data26_2[6], homopolymer_UVB_data26_2[7], homopolymer_UVB_data26_2[8])[0]
homopolymer26B2 = GetPolymerCounts(homopolymer_UVB_data26_2[0],  homopolymer_UVB_data26_2[2], homopolymer_UVB_data26_2[3], homopolymer_UVB_data26_2[4], homopolymer_UVB_data26_2[5], homopolymer_UVB_data26_2[6], homopolymer_UVB_data26_2[7], homopolymer_UVB_data26_2[8])[1]

homopolymer_combined_dataB = MakeMultiNucFile(a7b, a8b)
homopolymer_combinedB = GetPolymerCounts(homopolymer_combined_dataB[0],  homopolymer_combined_dataB[2], homopolymer_combined_dataB[3], homopolymer_combined_dataB[4], homopolymer_combined_dataB[5], homopolymer_combined_dataB[6], homopolymer_combined_dataB[7], homopolymer_combined_dataB[8])[1]
homopolymer_combined_data16B = MakeMultiNucFile(b7b, b8b)
homopolymer16_combinedB = GetPolymerCounts(homopolymer_combined_data16B[0],  homopolymer_combined_data16B[2], homopolymer_combined_data16B[3], homopolymer_combined_data16B[4], homopolymer_combined_data16B[5], homopolymer_combined_data16B[6], homopolymer_combined_data16B[7], homopolymer_combined_data16B[8])[1]
homopolymer_combined_data26B = MakeMultiNucFile(d7b, d8b)
homopolymer26_combinedB = GetPolymerCounts(homopolymer_combined_data26B[0],  homopolymer_combined_data26B[2], homopolymer_combined_data26B[3], homopolymer_combined_data26B[4], homopolymer_combined_data26B[5], homopolymer_combined_data26B[6], homopolymer_combined_data26B[7], homopolymer_combined_data26B[8])[1]

WT_dataB = homopolymer_combined_dataB
WT_isolateB = WT_dataB[4]
WT_chromosomeB = WT_dataB[2]
WT_position1B = WT_dataB[9]
WT_positionB = WT_dataB[3]
WT_trinucleotideB = WT_dataB[5]
WT_mutationB = WT_dataB[6]
WT_strandB = WT_dataB[1]
WT_multinucleotideB = WT_dataB[0] 
WT_lengthsB = WT_dataB[8]
WT_alleleB = WT_dataB[7]


rad16_dataB = homopolymer_combined_data16B
rad16_isolateB = rad16_dataB[4]
rad16_chromosomeB = rad16_dataB[2]
rad16_position1B = rad16_dataB[9]
rad16_positionB = rad16_dataB[3]
rad16_trinucleotideB = rad16_dataB[5]
rad16_mutationB = rad16_dataB[6]
rad16_strandB = rad16_dataB[1]
rad16_multinucleotideB = rad16_dataB[0] 
rad16_lengthsB = rad16_dataB[8]
rad16_alleleB = rad16_dataB[7]

rad26_dataB = homopolymer_combined_data26B

rad26_isolateB = rad26_dataB[4]
rad26_chromosomeB = rad26_dataB[2]
rad26_position1B = rad26_dataB[9]
rad26_positionB = rad26_dataB[3]
rad26_trinucleotideB = rad26_dataB[5]
rad26_mutationB = rad26_dataB[6]
rad26_strandB = rad26_dataB[1]
rad26_multinucleotideB = rad26_dataB[0] 
rad26_lengthsB = rad26_dataB[8]
rad26_alleleB = rad26_dataB[7]

rad30_dataB = homopolymer_UVB_data30
rad30_isolateB = rad30_dataB[4]
rad30_chromosomeB = rad30_dataB[2]
rad30_position1B = rad30_dataB[9]
rad30_positionB = rad30_dataB[3]
rad30_trinucleotideB = rad30_dataB[5]
rad30_mutationB = rad30_dataB[6]
rad30_strandB = rad30_dataB[1]
rad30_multinucleotideB = rad30_dataB[0] 
rad30_lengthsB = rad30_dataB[8]




MakeTriNucFile(a1b, a2b, a3b, a4b, WT_del_lengths)
MakeTriNucFile(b1b, b2b, b3b, b4b, rad16_del_lengths)   
MakeTriNucFile(d1b, d2b, d3b, d4b, rad26_del_lengths)
MakeInsertionFile(m1, m2, m3, m4, WT_ins_lengths)
MakeInsertionFile(h1b, h2b, h3b, h4b, rad16_ins_lengths)
MakeInsertionFile(j1b, j2b, j3b, j4b, rad26_ins_lengths)

m2.close()
m3.close()
m4.close()
h2b.close()
h3b.close()
h4b.close()
j2b.close()
j3b.close()
j4b.close()
k2b.close()
k3b.close()
k4b.close()
a1b.close()
b1b.close()
b2b.close()
b3b.close()
b4b.close()
e1b.close()
h1b.close()
j1b.close()
k1b.close()
a1b.close()
a2b.close()
a3b.close()
a4b.close()
m1.close()
a7.close()
a8.close()
b7.close()
b8.close()
d7.close()
d8.close()
a7b.close()
a8b.close()
b7b.close()
b8b.close()
d7b.close()
d8b.close()



MakeSequenceFile(sequencefile7, sequencefile8, WT_multinucleotideB, WT_mutationB, WT_alleleB, WT_strandB)

MakeSequenceFile(sequencefile9, sequencefile10, rad16_multinucleotideB, rad16_mutationB, rad16_alleleB, rad16_strandB)

MakeSequenceFile(sequencefile11, sequencefile12, rad26_multinucleotideB, rad26_mutationB, rad26_alleleB, rad26_strandB)


WT_isolate_listB = {}
for i1 in range(0, len(WT_isolateB)):

    if WT_isolateB[i1] not in WT_isolate_listB.keys():
            WT_isolate_listB[WT_isolateB[i1]] = [(WT_positionB[i1],  str(WT_chromosomeB[i1]), WT_trinucleotideB[i1], WT_mutationB[i1], WT_strandB[i1], WT_multinucleotideB[i1], WT_isolateB[i1], WT_position1B[i1])]
    else:
        WT_isolate_listB[WT_isolateB[i1]] .append((WT_positionB[i1],  str(WT_chromosomeB[i1]), WT_trinucleotideB[i1], WT_mutationB[i1], WT_strandB[i1], WT_multinucleotideB[i1], WT_isolateB[i1], WT_position1B[i1]))
rad16_isolate_listB = {}
for i1 in range(0, len(rad16_isolateB)):

    if rad16_isolateB[i1] not in rad16_isolate_listB.keys():
            rad16_isolate_listB[rad16_isolateB[i1]] = [(rad16_positionB[i1],  str(rad16_chromosomeB[i1]), rad16_trinucleotideB[i1], rad16_mutationB[i1], rad16_strandB[i1], rad16_multinucleotideB[i1], rad16_isolateB[i1], rad16_position1B[i1])]
    else:
        rad16_isolate_listB[rad16_isolateB[i1]] .append((rad16_positionB[i1],  str(rad16_chromosomeB[i1]), rad16_trinucleotideB[i1], rad16_mutationB[i1], rad16_strandB[i1], rad16_multinucleotideB[i1], rad16_isolateB[i1], rad16_position1B[i1]))

rad26_isolate_listB = {}
for i1 in range(0, len(rad26_isolateB)):

    if rad26_isolateB[i1] not in rad26_isolate_listB.keys():
            rad26_isolate_listB[rad26_isolateB[i1]] = [(rad26_positionB[i1],  str(rad26_chromosomeB[i1]), rad26_trinucleotideB[i1], rad26_mutationB[i1], rad26_strandB[i1], rad26_multinucleotideB[i1], rad26_isolateB[i1], rad26_position1B[i1])]
    else:
        rad26_isolate_listB[rad26_isolateB[i1]] .append((rad26_positionB[i1],  str(rad26_chromosomeB[i1]), rad26_trinucleotideB[i1], rad26_mutationB[i1], rad26_strandB[i1], rad26_multinucleotideB[i1], rad26_isolateB[i1], rad26_position1B[i1]))

rad30_isolate_listB = {}
for i1 in range(0, len(rad30_isolateB)):

    if rad30_isolateB[i1] not in rad30_isolate_listB.keys():
            rad30_isolate_listB[rad30_isolateB[i1]] = [(rad30_positionB[i1],  str(rad30_chromosomeB[i1]), rad30_trinucleotideB[i1], rad30_mutationB[i1], rad30_strandB[i1], rad30_multinucleotideB[i1], rad30_isolateB[i1], rad30_position1B[i1])]
    else:
        rad30_isolate_listB[rad30_isolateB[i1]] .append((rad30_positionB[i1],  str(rad30_chromosomeB[i1]), rad30_trinucleotideB[i1], rad30_mutationB[i1], rad30_strandB[i1], rad30_multinucleotideB[i1], rad30_isolateB[i1], rad30_position1B[i1]))

sequencefile7.close()
sequencefile8.close()
sequencefile9.close()
sequencefile10.close()
sequencefile11.close()
sequencefile12.close()


a1b.close()
a2b.close()
a3b.close()
a4b.close()
a5b.close()
a6b.close()
b1b.close()
b2b.close()
b3b.close()
b4b.close()
b5b.close()
b6b.close()
d1b.close()
d2b.close()
d3b.close()
d4b.close()
d5b.close()
d6b.close()
e1b.close()
e2b.close()
e3b.close()
e4b.close()
e5b.close()
e6b.close()
m1.close()
m6.close()
h1b.close()
h6b.close()
j1b.close()
j6b.close()
k1b.close()
k6b.close()




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

    # for l in range(0, len(lengths)):
    #     if lengths[l] not in homopolymer_counts.keys():
    #         homopolymer_counts[lengths[l]] = 1
    #         homopolymer_dict[lengths[l]] = 0
    #     else:
    #         homopolymer_counts[lengths[l]] = homopolymer_counts[lengths[l]] + 1

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
        # for  c in range(0, len(coordinates)):
        #     for p in range(0, len(pos1)):
        #         if chr[p] == key:
        #             if int(coordinates[c][0]) >= int(pos1[p]) or int(coordinates[c][1]) <= int(pos2[p]):
        #                 if strand_data[p] == '-':
        #                     for i in range(max(coordinates[c][0], int(pos1[p])), min(coordinates[c][1] + 1, int(pos2[p]))):
        #                         if sequence[i- 1] == 'C' or sequence[i- 1] == 'T':
        #                             homopolymer_TS = homopolymer_TS + 1
        #                 else:
        #                      for j in range(max(coordinates[c][0], int(pos1[p])), min(coordinates[c][1] + 1, int(pos2[p]))):
        #                         if sequence[j- 1] == 'C' or sequence[j- 1] == 'T':
        #                             homopolymer_NTS = homopolymer_NTS + 1

                       
                         
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




dicts = NormalizeHomopolymer(sequences, polymer_counts, homopolymer_file)
dicts2 = NormalizeHomopolymer(sequences, polymer_counts16, homopolymer_file)
dicts3 = NormalizeHomopolymer(sequences, polymer_counts26, homopolymer_file)
dicts4 = NormalizeHomopolymer(sequences, polymer_counts_UVB, homopolymer_file)
dicts5 = NormalizeHomopolymer(sequences, polymer_counts_UVB16, homopolymer_file)
dicts6 = NormalizeHomopolymer(sequences, polymer_counts_UVB26, homopolymer_file)




def index(tuple):
    return tuple[1]

def FindClusters( isolate_list, list2, strain, homopolymers): #Find indels that are <= 10 base pairs from another indel in the same isolate
    lengths = {}
    clusters = {}
    nonclusters = {}
    total_clusters = 0
    indel_clusters = 0
    list2_key = {}
    for key in isolate_list.keys():
        if strain in key[1]:
            sorted_isolate_list = sorted(isolate_list[key])
            sorted_list = sorted(sorted_isolate_list, key = index)
            if key[0] in list2.keys():
                list2_key = list2[key[0]]
            new_list2 = {}
            for i in range(0, len(list2_key)):
                for j in range(0, len(sorted_list)):
                    if sorted_list[j][1] == list2_key[i][1] and sorted_list[j][0] == list2_key[i][0]:
                        new_list2[j] = list2_key[i]

            if sorted_list[0][1] == sorted_list[1][1] and (int(sorted_list[1][0]) <= int((sorted_list[ 0][0])) + 10):# and (int(sorted_list[1][0]) != int((sorted_list[ 0][0]))):
                if (sorted_list[0][6], sorted_list[0][1], sorted_list[0][0]) not in clusters.keys() and (sorted_list[0][5] == "DIP" or sorted_list[0][5] == "Deletion" or sorted_list[0][5] == "Insertion"):
                #if (sorted_list[0][6], sorted_list[0][1], sorted_list[0][0]) not in clusters.keys() and (sorted_list[0][5] == "SNV" or sorted_list[0][5] == "SNP"):    
                    clusters[(sorted_list[0][6], sorted_list[0][1], sorted_list[0][0])] = [sorted_list[0]]
                 

                if (sorted_list[1][6], sorted_list[1][1], sorted_list[1][0]) not in clusters.keys() and (sorted_list[1][5] == "DIP" or sorted_list[1][5] == "Deletion" or sorted_list[1][5] == "Insertion"):
                #if (sorted_list[1][6], sorted_list[1][1], sorted_list[1][0]) not in clusters.keys() and (sorted_list[1][5] == "SNV" or sorted_list[1][5] == "SNP"):
                    clusters[(sorted_list[1][6], sorted_list[1][1], sorted_list[1][0])] = [sorted_list[1]]
                    
            else:#elif (int(sorted_list[1][0]) != int((sorted_list[ 0][0]))):
                if ('-' not in sorted_list[0][3] and sorted_list[0][5] == "DIP") or sorted_list[0][5] == "Insertion":
                    if sorted_list[0][4] == '-':
                        DNA_strand = '-'
                        end_position = GetHomopolymerLengthIns(GetReverseComplement(new_list2[0][5]), int(new_list2[0][7]), int(new_list2[0][0]), int(len(new_list2[0][5])/2) + 1, DNA_strand)[1][1] 
                    else:
                        DNA_strand = '+'
                        end_position = GetHomopolymerLengthIns(new_list2[0][5], int(new_list2[0][7]), int(new_list2[0][0]), int(len(new_list2[0][5])/2) + 1, DNA_strand)[1][1] 
                else:#if (sorted_list[0][5] == "DIP" or sorted_list[0][5] == "Complex DIP" or sorted_list[0][5] == "Deletion" or sorted_list[0][5] == "Insertion") and 0 in new_list2.keys():
                    if sorted_list[0][4] == '-':#Move position to end of homopolymer since it could be anywhere in the homopolymer
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
                    #if (sorted_list[1][5] == "SNV" or sorted_list[1][5] == "SNP"):
                        if (sorted_list[1][6], sorted_list[1][1], sorted_list[1][0]) not in clusters.keys():
                            clusters[(sorted_list[1][6], sorted_list[1][1], sorted_list[1][0])] = [sorted_list[1]]
                    if (sorted_list[0][5] == "DIP" or sorted_list[0][5] == "Deletion" or sorted_list[0][5] == "Insertion"):
                    #if (sorted_list[0][5] == "SNV" or sorted_list[0][5] == "SNP"):
                        if (sorted_list[0][6], sorted_list[0][1], sorted_list[0][0]) not in clusters.keys():
                            clusters[(sorted_list[0][6], sorted_list[0][1], sorted_list[0][0])] = [sorted_list[0]]
            
            if sorted_list[len (sorted_list) - 2][1] == sorted_list[len (sorted_list) - 1][1] and (int(sorted_list[len (sorted_list) - 1][0]) <= int((sorted_list[ len (sorted_list) - 2][0])) + 10):#and (int(sorted_list[len (sorted_list) - 1][0]) != int((sorted_list[ len (sorted_list) - 2][0]))):
                if (sorted_list[len(sorted_list) - 1][6], sorted_list[len(sorted_list) - 1][1], sorted_list[len(sorted_list) - 1][0]) not in clusters.keys() and (sorted_list[len(sorted_list) - 1][5] == "DIP" or sorted_list[len(sorted_list) - 1][5] == "Deletion" or sorted_list[len(sorted_list) - 1][5] == "Insertion"):
                #if (sorted_list[len(sorted_list) - 1][6], sorted_list[len(sorted_list) - 1][1], sorted_list[len(sorted_list) - 1][0]) not in clusters.keys() and (sorted_list[len(sorted_list) - 1][5] == "SNV" or sorted_list[len(sorted_list) - 1][5] == "SNP"):    
                    clusters[(sorted_list[len(sorted_list) - 1][6], sorted_list[len(sorted_list) - 1][1], sorted_list[len(sorted_list) - 1][0])] = [sorted_list[len(sorted_list) - 1]]
                   
                if (sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][0]) not in clusters.keys() and (sorted_list[len(sorted_list) - 2][5] == "DIP" or sorted_list[len(sorted_list) - 2][5] == "Deletion" or sorted_list[len(sorted_list) - 2][5] == "Insertion"):
                #if (sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][0]) not in clusters.keys() and (sorted_list[len(sorted_list) - 2][5] == "SNV" or sorted_list[len(sorted_list) - 2][5] == "SNP"): 
                    clusters[(sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][0])] = [sorted_list[len(sorted_list) - 2]]
                   
            elif (int(sorted_list[len (sorted_list) - 1][0]) != int((sorted_list[ len (sorted_list) - 2][0]))):#if strand[0] == '+':
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
                    #if (sorted_list[len(sorted_list) - 2][5] == "SNV" or sorted_list[len(sorted_list) - 2][5] == "SNP"):
                        if (sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][6]) not in clusters.keys():
                            clusters[(sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][6])] = [sorted_list[len(sorted_list) - 2]]
                    if (sorted_list[len(sorted_list) - 1][5] == "DIP" or sorted_list[len(sorted_list) - 1][5] == "Deletion" or sorted_list[len(sorted_list) - 1][5] == "Insertion"):  
                    #if (sorted_list[len(sorted_list) - 1][5] == "SNV" or sorted_list[len(sorted_list) - 1][5] == "SNP"): 
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
                    #if (sorted_list[len(sorted_list) - 1][5] == "SNV" or sorted_list[len(sorted_list) - 1][5] == "SNP"):
                        if (sorted_list[len(sorted_list) - 1][6], sorted_list[len(sorted_list) - 1][1], sorted_list[len(sorted_list) - 1][0]) not in clusters.keys():
                            clusters[(sorted_list[len(sorted_list) - 1][6], sorted_list[len(sorted_list) - 1][1], sorted_list[len(sorted_list) - 1][0])] = [sorted_list[len(sorted_list) - 1]]
                    if (sorted_list[len(sorted_list) - 2][5] == "DIP" or sorted_list[len(sorted_list) - 2][5] == "Deletion" or sorted_list[len(sorted_list) - 2][5] == "Insertion"):
                    #if (sorted_list[len(sorted_list) - 2][5] == "SNV" or sorted_list[len(sorted_list) - 2][5] == "SNP"): 
                        if (sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][0]) not in clusters.keys():
                            clusters[(sorted_list[len(sorted_list) - 2][6], sorted_list[len(sorted_list) - 2][1], sorted_list[len(sorted_list) - 2][0])] = [sorted_list[len(sorted_list) - 2]]
                    
            for j in range (1,len (sorted_list) - 1):
                #if int(sorted_list[j][0]) != int((sorted_list[j + 1][0])):
                    if (sorted_list[j][1] == sorted_list[j - 1][1] and (int(sorted_list[j][0]) <= int((sorted_list[j - 1][0])) + 10)) or (sorted_list[j + 1][1] == sorted_list[j][1] and (int(sorted_list[j + 1][0]) <= int((sorted_list[j][0])) + 10)): #or ((isolate_list[key][j + n][6], int(isolate_list[key][j + n][0]), isolate_list[key][j + n][1]) in homopolymers.keys() and (int(isolate_list[key][j + n][0]) <= int((isolate_list[key][j + n - 1][0])) + 10 + homopolymers[key])) :  
                            if (sorted_list[j][6], sorted_list[j][1], sorted_list[j][0]) not in clusters.keys() and (sorted_list[j][5] == "DIP" or sorted_list[j][5] == "Deletion" or sorted_list[j][5] == "Insertion"):
                            #if (sorted_list[j][6], sorted_list[j][1], sorted_list[j][0]) not in clusters.keys() and (sorted_list[j][5] == "SNV" or sorted_list[j][5] == "SNP"):    
                                #if j == 0 or isolate_list[key][j - n] not in clusters.keys():
                                clusters[(sorted_list[j][6], sorted_list[j][1], sorted_list[j][0])] = [sorted_list[j]]
                                #total_clusters = total_clusters + 1
                    else:#elif int(sorted_list[j][0]) != int((sorted_list[j + 1][0])) and int(sorted_list[j][0]) != int((sorted_list[j - 1][0])):
                        if (sorted_list[j][5] == "DIP"  and '-' not in sorted_list[j][3])or sorted_list[j][5] == "Insertion": 
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
                            #if (sorted_list[j][5] == "SNV"or sorted_list[j][5] == "SNP"):
                                    #if sorted_list[j][4] == '-':
                                if (sorted_list[j][6], sorted_list[j][1], sorted_list[j][0]) not in clusters.keys():
                                    clusters[(sorted_list[j][6], sorted_list[j][1], sorted_list[j][0])] = [sorted_list[j]]
                            if (sorted_list[j + 1][5] == "DIP"  or sorted_list[j + 1][5] == "Deletion" or sorted_list[j + 1][5] == "Insertion"): 
                            #if (sorted_list[j + 1][5] == "SNV"  or sorted_list[j + 1][5] == "SNP"):
                                if (sorted_list[j + 1][6], sorted_list[j + 1][1], sorted_list[j + 1][0]) not in clusters.keys():
                                    clusters[(sorted_list[j + 1][6], sorted_list[j + 1][1], sorted_list[j + 1][0])] = [sorted_list[j + 1]]
                        if (sorted_list[j][5] == "DIP"  and '-' not in sorted_list[j][3])or sorted_list[j][5] == "Insertion":
                            if sorted_list[j][4] == '-':
                                DNA_strand = "-"
                                start_position = GetHomopolymerLengthIns(GetReverseComplement(new_list2[j][5]), int(new_list2[j][7]), int(new_list2[j][0]), int(len(new_list2[j][5])/2) + 1, DNA_strand)[1][0] 
                            else:
                                DNA_strand = "+"
                                start_position = GetHomopolymerLengthIns(new_list2[j][5], int(new_list2[j][7]), int(new_list2[j][0]), int(len(new_list2[j][5])/2) + 1, DNA_strand)[1][0] 
                        elif int(sorted_list[j][0]) != int((sorted_list[j + 1][0])) and int(sorted_list[j][0]) != int((sorted_list[j - 1][0])):    
                            if sorted_list[j][4] == '-':
                                start_position = GetHomopolymerLength(GetReverseComplement(new_list2[j][5]), int(new_list2[j][7]), int(new_list2[j][0]), int(len(new_list2[j][5])/2) + 1)[1][0] 
                            else:
                                start_position = GetHomopolymerLength(new_list2[j][5], int(new_list2[j][7]), int(new_list2[j][0]), int(len(new_list2[j][5])/2) + 1)[1][0] 
                        if start_position <= int((sorted_list[j - 1][0] + 10)) and sorted_list[j - 1][1] == sorted_list[j][1]:
                            if (sorted_list[j][5] == "DIP" or sorted_list[j][5] == "Deletion" or sorted_list[j][5] == "Insertion"): 
                            #if (sorted_list[j][5] == "SNV" or sorted_list[j][5] == "SNP"): 
                                if (sorted_list[j][6], sorted_list[j][1], sorted_list[j][0]) not in clusters.keys():
                                    clusters[(sorted_list[j][6], sorted_list[j][1], sorted_list[j][0])] = [sorted_list[j]]  #indel_clusters = indel_clusters + 1
                            if (sorted_list[j - 1][5] == "DIP" or sorted_list[j - 1][5] == "Deletion" or sorted_list[j - 1][5] == "Insertion"): 
                            #if (sorted_list[j - 1][5] == "SNV" or sorted_list[j - 1][5] == "SNP"):    
                                if (sorted_list[j - 1][6], sorted_list[j - 1][1], sorted_list[j - 1][0]) not in clusters.keys():
                                    clusters[(sorted_list[j - 1][6], sorted_list[j - 1][1], sorted_list[j - 1][0])] = [sorted_list[j - 1]]
                                #elif isolate_list[key][j] not in clusters.keys():
                                #total_clusters = total_clusters + 1
                    
                        
            for k in range(0, len(sorted_list)):
                if (sorted_list[k][5] == "DIP" or sorted_list[k][5] == "Deletion" or sorted_list[k][5] == "Insertion"):
                #if (sorted_list[k][5] == "SNV" or sorted_list[k][5] == "SNP"):
                    if (sorted_list[k][6], sorted_list[k][1], sorted_list[k][0]) not in clusters.keys() and (sorted_list[k][6], sorted_list[k][1], sorted_list[k][0]) not in nonclusters.keys():
                        nonclusters[(sorted_list[k][6], sorted_list[k][1], sorted_list[k][0])] = [sorted_list[k]]

#Count the number of complex homopolymer and non-homopolymer deletions and insertions
    homopolymer_deletions = {}
    homopolymer_insertions = {}
    homopolymer_subs = {}
    complex_deletions = []
    complex_insertions = []
    complex_subs = []
    for key in clusters.keys():
        if len(clusters[key][0][7]) == 1:
            if homopolymers[key] not in homopolymer_subs.keys():
                homopolymer_subs[homopolymers[key]] = 1
            else:
                homopolymer_subs[homopolymers[key]] = homopolymer_subs[homopolymers[key]] + 1
        if '-' in clusters[key][0][3] and len(clusters[key][0][7]) == 1:
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
        if clusters[key][0][3] == '-' and len(clusters[key][0][7]) == 1:
            single_deletion.append(clusters[key])
        elif '-' in clusters[key][0][3]:
            multi_deletion.append(clusters[key])

        elif len(clusters[key][0][3]) == 1:
            single_insertion.append(clusters[key])
        else:
            multi_insertion.append(clusters[key])
    return clusters, nonclusters

clusters = FindClusters(isolate_list, WT_isolate_list, 'WT', homopolymer_combined)
clusters16 = FindClusters(isolate_list, rad16_isolate_list, 'rad16', homopolymer16_combined)
clusters26 = FindClusters(isolate_list, rad26_isolate_list, 'rad26', homopolymer26_combined)
#clusters30 = FindClusters(isolate_list, rad30_isolate_list, 'rad30', homopolymer30)
clustersB = FindClusters(isolate_list2, WT_isolate_listB, 'WT', homopolymer_combinedB)
clusters16B = FindClusters(isolate_list2, rad16_isolate_listB, 'Rad16', homopolymer16_combinedB)
clusters26B = FindClusters(isolate_list2, rad26_isolate_listB, 'Rad26', homopolymer26_combinedB)
#clusters30B = FindClusters(isolate_list2, rad30_isolate_listB, 'Rad30', homopolymer30B)

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
    
a1 = open('WT_indel_clustersv2.bed', 'w+')
a2 = open('WT_indel_sequences', 'w+')
a3 = open('WT_indel_nonclusters.bed', 'w+')
a4 = open('WT_indel_sequences2', 'w+')
b1 = open('rad16_indel_clustersv2.bed', 'w+')
b2 = open('rad16_indel_sequences', 'w+')
b3 = open('rad16_indel_nonclusters.bed', 'w+')
b4 = open('rad16_indel_sequences2', 'w+')
c1 = open('rad26_indel_clustersv2.bed', 'w+')
c2 = open('rad26_indel_sequences', 'w+')
c3 = open('rad26_indel_nonclusters.bed', 'w+')
c4 = open('rad26_indel_sequences2', 'w+')
d1 = open('rad30_indel_clusters.bed', 'w+')
d2 = open('rad30_indel_sequences', 'w+')
d3 = open('rad30_indel_nonclusters.bed', 'w+')
d4 = open('rad30_indel_sequences2', 'w+')

e1 = open('WT_indel_clustersUVB.bed', 'w+')
e2 = open('WT_indel_sequencesUVB', 'w+')
e3 = open('WT_indel_nonclustersUVB.bed', 'w+')
e4 = open('WT_indel_sequencesUVB2', 'w+')
f1 = open('rad16_indel_clustersUVB.bed', 'w+')
f2 = open('rad16_indel_sequencesUVB', 'w+')
f3 = open('rad16_indel_nonclustersUVB.bed', 'w+')
f4 = open('rad16_indel_sequencesUVB2', 'w+')
g1 = open('rad26_indel_clustersUVB.bed', 'w+')
g2 = open('rad26_indel_sequencesUVB', 'w+')
g3 = open('rad26_indel_nonclustersUVB.bed', 'w+')
g4 = open('rad26_indel_sequencesUVB2', 'w+')
h1 = open('rad30_indel_clustersUVB.bed', 'w+')
h2 = open('rad30_indel_sequencesUVB', 'w+')
h3 = open('rad30_indel_nonclustersUVB.bed', 'w+')
h4 = open('rad30_indel_sequencesUVB2', 'w+')



def MakeClusterFiles(clusters, nonclusters, file1, file2, file3, file4):
    cluster_keys = sorted(clusters.keys(), key = index2)
    sorted_list = sorted(cluster_keys, key = index)
    file2.write("Deletions: ")
    file2.write('\n')
    for c in sorted_list:
        if '-' in clusters[c][0][3] and len(clusters[c][0][7]) == 1:
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
                    if clusters[c][0][4] == '-' and  len(clusters[c][0][7]) == 1:
                        file2.write(GetReverseComplement(seq[::-1]))
                    elif len(clusters[c][0][7]) == 1:
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
        if '-' in nonclusters[d][0][3] and len(nonclusters[d][0][7]) == 1:
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

                    if nonclusters[d][0][4] == '-' and len(nonclusters[d][0][7]) == 1:
                        file4.write(GetReverseComplement(seq[::-1]))
                    elif len(nonclusters[d][0][7]) == 1:
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

                
MakeClusterFiles(clusters[0], clusters[1], a1, a2, a3, a4)
MakeClusterFiles(clusters16[0], clusters16[1],b1, b2, b3, b4)
MakeClusterFiles(clusters26[0], clusters26[1], c1, c2, c3, c4)
#MakeClusterFiles(clusters30[0], clusters30[1], d1, d2, d3, d4)
MakeClusterFiles(clustersB[0], clustersB[1], e1, e2, e3, e4)
MakeClusterFiles(clusters16B[0], clusters16B[1], f1, f2, f3, f4)
MakeClusterFiles(clusters26B[0], clusters26B[1],g1, g2, g3, g4)
#MakeClusterFiles(clusters30B[0], clusters30B[1], h1, h2, h3, h4)

a1.close()
a2.close()
a3.close()
a4.close()
b1.close()
b2.close()
b3.close()
b4.close()
c1.close()
c2.close()
c3.close()
c4.close()
d1.close()
d2.close()
d3.close()
e1.close()
e2.close()
e3.close()
e4.close()
f1.close()
f2.close()
f3.close()
f4.close()
g1.close()
g2.close()
g3.close()
g4.close()
h1.close()
h2.close()
h3.close()
h4.close()








