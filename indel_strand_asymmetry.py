


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

def GetHomopolymerLength(sequence, mut_pos1, mut_pos, mut_position):
    #mut_position = int(len(sequence)/2) + 1
    mut_pos1 = int(mut_pos1)
    mut_pos = int(mut_pos) 
    difference = mut_pos - mut_pos1 - 1
    mut_position2 = mut_position - difference
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


def GetSequence(file, chr_name):
    f1 = open(file)
    
    sequence_lines = []
    empty = False
    start_line = -1
    line_number = 0
    while not empty:
        line = f1.readline()

        if chr_name in line:
            start_line = line_number
        elif 'chr' in line and start_line != -1:
            end_line = line_number
            break
        if start_line != -1: 
            sequence_lines.append(str(line))
        if line == "":
            empty = True
        line_number = line_number + 1
    sequence = ''.join(sequence_lines).upper().strip()
    for n in range(6, len(sequence)):
        if sequence[n] == 'N' or sequence[n] == 'A' or sequence[n] == 'C' or sequence[n] == 'G' or sequence[n] == 'T':
            sequence = sequence[n:len(sequence)]
            break
    f1.close()
    print("sequence finished")
    return sequence.replace('\n', '')

sequences = []
sequenceX = GetSequence('hg19.fa', 'chrX')
sequences.append(sequenceX)
sequence1 = GetSequence('hg19.fa', 'chr1')
sequences.append(sequence1)
sequence2 = GetSequence('hg19.fa', 'chr2')
sequences.append(sequence2)
sequence3 = GetSequence('hg19.fa', 'chr3')
sequences.append(sequence3)
sequence4 = GetSequence('hg19.fa', 'chr4')
sequences.append(sequence4)
sequence5 = GetSequence('hg19.fa', 'chr5')
sequences.append(sequence5)
sequence6 = GetSequence('hg19.fa', 'chr6')
sequences.append(sequence6)
sequence7 = GetSequence('hg19.fa', 'chr7')
sequences.append(sequence7)
sequence8 = GetSequence('hg19.fa', 'chr8')
sequences.append(sequence8)
sequence9 = GetSequence('hg19.fa', 'chr9')
sequences.append(sequence9)
sequence10 = GetSequence('hg19.fa', 'chr10')
sequences.append(sequence10)
sequence11 = GetSequence('hg19.fa', 'chr11')
sequences.append(sequence11)
sequence12 = GetSequence('hg19.fa', 'chr12')
sequences.append(sequence12)
sequence13 = GetSequence('hg19.fa', 'chr13')
sequences.append(sequence13)
sequence14 = GetSequence('hg19.fa', 'chr14')
sequences.append(sequence14)
sequence15 = GetSequence('hg19.fa', 'chr15')
sequences.append(sequence15)
sequence16 = GetSequence('hg19.fa', 'chr16')
sequences.append(sequence16)
sequence17 = GetSequence('hg19.fa', 'chr17')
sequences.append(sequence17)
sequence18 = GetSequence('hg19.fa', 'chr18')
sequences.append(sequence18)
sequence19 = GetSequence('hg19.fa', 'chr19')
sequences.append(sequence19)
sequence20 = GetSequence('hg19.fa', 'chr20')
sequences.append(sequence20)
sequence21 = GetSequence('hg19.fa', 'chr21')
sequences.append(sequence21)
sequence22 = GetSequence('hg19.fa', 'chr22')
sequences.append(sequence22)
sequence23 = GetSequence('hg19.fa', 'chrY')
sequences.append(sequence23)

with open('lowrpkm_melanocyte01_E059_gencodev10_hg19_sorted.bed') as f:
    data = f.read()
lines = [s.strip().split() for s in data.splitlines()]
chr = []
pos1 = []
pos2 = []
strand_data = []
for line in lines:
    chr.append(line[0])
    pos1.append(int(line[1]))
    pos2.append(int(line[2]))
    strand_data.append(line[5])

with open('highrpkm_melanocyte01_E059_gencodev10_hg19_sorted.bed') as f1:
    data2 = f1.read()
lines2 = [s.strip().split() for s in data2.splitlines()]
chr2 = []
pos1_2 = []
pos2_2 = []
strand_data2 = []
for line in lines2:
    chr2.append(line[0])
    pos1_2.append(int(line[1]))
    pos2_2.append(int(line[2]))
    strand_data2.append(line[5])

with open("Acral_singlebase_deletions_all.bed") as a:  
    acral_data = a.read()
with open("Cutaneous_singlebase_deletions_all.bed") as c:
    cutaneous_data = c.read()

with open("Acral_singlebase_insertions.bed") as ai:  
    acral_data_ins = ai.read()
with open("Cutaneous_singlebase_insertions.bed") as i:
    cutaneous_data_ins = i.read()
i
with open("Acral_singlebase_deletions_homopolymer.bed") as ah:
     a_homopolymer_data = ah.read()
with open("Acral_singlebase_insertions_homopolymer.bed") as an:
     a_homopolymer_data2 = an.read()

with open("Cutaneous_singlebase_deletions_homopolymer.bed") as ch:
     c_homopolymer_data = ch.read()
with open("Cutaneous_singlebase_insertions_homopolymer.bed") as cn:
     c_homopolymer_data2 = cn.read()

def getChromosome(data):       
    lines = [s.strip().split() for s in data.splitlines()]
    chr_counts = {}
    chromosome = []
    for i in range (0, len(lines)):
        #if lines[i][0] != "chrXII" or (int(lines[i][2]) < 451430) or (int(lines[i][2]) > 468920 and int(lines[i][2]) < 489928) or (int(lines[i][2]) > 490455):
            chromosome.append(lines[i][0])
            if lines[i][0] not in chr_counts.keys():
                chr_counts[lines[i][0]] = 1
            else:
                chr_counts[lines[i][0]] = chr_counts[lines[i][0]] + 1
    for key in chr_counts.keys():
        print(key + ": " + str(chr_counts[key]))
    return chromosome
chromosome = getChromosome(acral_data)
chromosome_cutaneous = getChromosome(cutaneous_data)
a_homopolymer_chromosome = getChromosome(a_homopolymer_data)
a_homopolymer_chromosome2 = getChromosome(a_homopolymer_data2)
c_homopolymer_chromosome = getChromosome(c_homopolymer_data)
c_homopolymer_chromosome2 = getChromosome(c_homopolymer_data2)
chromosome_ins = getChromosome(acral_data_ins)
chromosome_cutaneous_ins = getChromosome(cutaneous_data_ins)



def getPosition1(data):
    lines = [s.strip().split() for s in data.splitlines()]   
    position1 = []
    for j in range (0, len(lines)):
        #if lines[j][0] != "chrXII" or (int(lines[j][2]) < 451430) or (int(lines[j][2]) > 468920 and int(lines[j][2]) < 489928) or (int(lines[j][2]) > 490455):
            position1.append(float(lines[j][1]))
    return position1
position1 = getPosition1(acral_data)
position1_cutaneous = getPosition1(cutaneous_data)
position1_a_homopolymer = getPosition1(a_homopolymer_data)
position1_a_homopolymer2 = getPosition1(a_homopolymer_data2)
position1_c_homopolymer = getPosition1(c_homopolymer_data)
position1_c_homopolymer2 = getPosition1(c_homopolymer_data2)
position1_ins = getPosition1(acral_data_ins)
position1_cutaneous_ins = getPosition1(cutaneous_data_ins)


def getPosition2(data):
    lines = [s.strip().split() for s in data.splitlines()]   
    position2 = []
    for k in range (0, len(lines) ):
        #if lines[k][0] != "chrXII" or (int(lines[k][2]) < 451430) or (int(lines[k][2]) > 468920 and int(lines[k][2]) < 489928) or (int(lines[k][2]) > 490455):
            position2.append(float(lines[k][2]))
    return position2
position2 = getPosition2(acral_data)
position2_cutaneous = getPosition2(cutaneous_data)
position2_a_homopolymer = getPosition2(a_homopolymer_data)
position2_a_homopolymer2 = getPosition2(a_homopolymer_data2)
position2_c_homopolymer = getPosition2(c_homopolymer_data)
position2_c_homopolymer2 = getPosition2(c_homopolymer_data2)
position2_ins = getPosition2(acral_data_ins)
position2_cutaneous_ins = getPosition2(cutaneous_data_ins)


def getTrinucleotide(data):
    lines = [s.strip().split() for s in data.splitlines()]   
    trinucleotide = []
    for m in range (0, len(lines)):
        #if lines[m][0] != "chrXII" or (int(lines[m][2]) < 451430) or (int(lines[m][2]) > 468920 and int(lines[m][2]) < 489928) or (int(lines[m][2]) > 490455):
            trinucleotide.append(lines[m][3].replace('-', ''))
    return trinucleotide
trinucleotide = getTrinucleotide(acral_data)
trinucleotide_cutaneous = getTrinucleotide(cutaneous_data)
trinucleotide_a_homopolymer = getTrinucleotide(a_homopolymer_data)
trinucleotide_a_homopolymer2 = getTrinucleotide(a_homopolymer_data2)
trinucleotide_c_homopolymer = getTrinucleotide(c_homopolymer_data)
trinucleotide_c_homopolymer2 = getTrinucleotide(c_homopolymer_data2)
trinucleotide_ins = getTrinucleotide(acral_data_ins)
trinucleotide_cutaneous_ins = getTrinucleotide(cutaneous_data_ins)

def getMutation(data):
    lines = [s.strip().split() for s in data.splitlines()]   
    mutation = []
    for n in range(0, len(lines)):
            mutation.append(lines[n][4])
    return mutation

mutation = getMutation(acral_data)
mutation_cutaneous = getMutation(cutaneous_data)
mutation_a_homopolymer = getMutation(a_homopolymer_data)
mutation_a_homopolymer2 = getMutation(a_homopolymer_data2)
mutation_c_homopolymer = getMutation(c_homopolymer_data)
mutation_c_homopolymer2 = getMutation(c_homopolymer_data2)
mutation_ins = getMutation(acral_data_ins)
mutation_cutaneous_ins = getMutation(cutaneous_data_ins)

def getStrand(data):
    lines = [s.strip().split() for s in data.splitlines()]   
    strand = []
    for n in range (0, len(lines)):
            strand.append(lines[n][5])
    return strand

strand = getStrand(acral_data)
strand_cutaneous = getStrand(cutaneous_data)
strand_a_homopolymer = getStrand(a_homopolymer_data)
strand_a_homopolymer2 = getStrand(a_homopolymer_data2)
strand_c_homopolymer = getStrand(c_homopolymer_data)
strand_c_homopolymer2 = getStrand(c_homopolymer_data2)
strand_ins = getStrand(acral_data_ins)
strand_cutaneous_ins = getStrand(cutaneous_data_ins)

high_file = open('high_expression.bed', 'w+')
low_file = open('low_expression.bed', 'w+')
acral_high_file = open('acral_highexpression.bed', 'w+')
acral_low_file = open('acral_lowexpression.bed', 'w+')

high_homopolymer_file = open('high_expression_homopolymer.bed', 'w+')
low_homopolymer_file = open('low_expression_homopolymer.bed', 'w+')
acral_high_homopolymer_file = open('acral_highexpression_homopolymer.bed', 'w+')
acral_low_homopolymer_file = open('acral_lowexpression_homopolymer.bed', 'w+')

high_homopolymer_file2 = open('high_expression_homopolymer2.bed', 'w+')
low_homopolymer_file2 = open('low_expression_homopolymer2.bed', 'w+')
acral_high_homopolymer_file2 = open('acral_highexpression_homopolymer2.bed', 'w+')
acral_low_homopolymer_file2 = open('acral_lowexpression_homopolymer2.bed', 'w+')

high_file_ins = open('high_expression_ins.bed', 'w+')
low_file_ins = open('low_expression_ins.bed', 'w+')
acral_high_file_ins = open('acral_highexpression_ins.bed', 'w+')
acral_low_file_ins = open('acral_lowexpression_ins.bed', 'w+')


def FindAsymmetry(file, chromosome, position1, position2, trinucleotide, strand, chr, pos1, pos2, strand_data, mutation):

    intergenic = []
    intergenic_count = 0
    NTS = []
    NTS_count = 0
    C_NTS = 0
    T_NTS = 0
    TS = []
    C_TS = 0
    T_TS = 0
    TS_count = 0
    NTS_sequence_contexts = {}
    TS_sequence_contexts = {}


    po = 0
    p2 = 1
    #Goes through mutation positions and gene positions in parallel to check whether a mutation is in a gene
    while po <= len(position2) - 1  and p2 <= len(pos2): 
            if position2[po] > pos1[p2 - 1] and position2[po] <= pos2[p2 -1] and chromosome[po] == chr[p2 - 1]:
                if chromosome[po] == 'chrY':
                     print('Y')
            #falls within gene 
            #The template will be in the opposite orientation as the transcript.
                if strand[po] == strand_data[p2 - 1]:
                    NTS.append([chromosome[po], position1[po], position2[po], trinucleotide[po], strand[po]])
                    NTS_count = NTS_count + 1
                    if trinucleotide[po][1] == 'C' or trinucleotide[po][1] == 'G':
                         C_NTS = C_NTS + 1
                    if trinucleotide[po][1] == 'T' or trinucleotide[po][1]== 'A':
                         T_NTS = T_NTS + 1
                    if trinucleotide[po] in NTS_sequence_contexts.keys():
                         NTS_sequence_contexts[trinucleotide[po]] = NTS_sequence_contexts[trinucleotide[po]] + 1
                    else:
                         NTS_sequence_contexts[trinucleotide[po]] = 1
        
                else:
                    TS.append([chromosome[po], position1[po], position2[po], trinucleotide[po],  strand[po]])
                    TS_count = TS_count + 1
                    if trinucleotide[po][1] == 'C' or trinucleotide[po][1] == 'G':
                         C_TS = C_TS + 1
                    if trinucleotide[po][1] == 'T' or trinucleotide[po][1] == 'A':
                         T_TS = T_TS + 1
                    if trinucleotide[po] in TS_sequence_contexts.keys():
                         TS_sequence_contexts[trinucleotide[po]] = TS_sequence_contexts[trinucleotide[po]] + 1
                    else:
                         TS_sequence_contexts[trinucleotide[po]] = 1

                po = po + 1
                file.write(chromosome[po] + '\t' + str(position1[po]) + '\t' + str(position2[po]) + '\t' + trinucleotide[po] + '\t' + mutation[po] + '\t' + strand[po])
                file.write('\n')
        #Falls between genes
            elif position2[po] > pos2[p2 - 1] and position2[po] <= pos1[p2 - 1] and chromosome[po] == chr[p2 - 1]:
                intergenic.append([chromosome[po], position1[po], position2[po], trinucleotide[po], strand[po]])
                intergenic_count = intergenic_count + 1
                po = po + 1

            elif chromosome[po] < chr[p2 - 1]:
                intergenic.append([chromosome[po], position1[po], position2[po], trinucleotide[po], strand[po]])
                intergenic_count = intergenic_count + 1
                po = po + 1

            #Move on to the next chromosome    
            elif chromosome[po] > chr[p2 - 1]:
                p2 = p2 + 1
        #Mutation position is greater than the end position of the current gene, so we move on to the next one.
            elif position2[po] > pos2[p2 - 1] :
                p2 = p2 + 1

            elif position2[po] < pos2[p2-1]:
                intergenic.append([chromosome[po], position1[po], position2[po], trinucleotide[po],  strand[po]])
                intergenic_count = intergenic_count + 1
                po = po + 1
    return T_NTS, C_NTS, T_TS, C_TS, NTS_sequence_contexts, TS_sequence_contexts



low_acral = FindAsymmetry(acral_low_file, chromosome, position1, position2, trinucleotide, strand, chr, pos1, pos2, strand_data, mutation)
high_acral = FindAsymmetry(acral_high_file,chromosome, position1, position2, trinucleotide, strand, chr2, pos1_2, pos2_2, strand_data2, mutation)
low_cutaneous = FindAsymmetry(low_file, chromosome_cutaneous, position1_cutaneous, position2_cutaneous, trinucleotide_cutaneous, strand_cutaneous, chr, pos1, pos2, strand_data, mutation_cutaneous)
high_cutaneous = FindAsymmetry(high_file, chromosome_cutaneous, position1_cutaneous, position2_cutaneous, trinucleotide_cutaneous, strand_cutaneous, chr2, pos1_2, pos2_2, strand_data2, mutation_cutaneous)
         
low_acral_homopolymer = FindAsymmetry(acral_low_homopolymer_file, a_homopolymer_chromosome, position1_a_homopolymer, position2_a_homopolymer, trinucleotide_a_homopolymer, strand_a_homopolymer, chr, pos1, pos2, strand_data, mutation_a_homopolymer)
high_acral_homopolymer = FindAsymmetry(acral_high_homopolymer_file,a_homopolymer_chromosome, position1_a_homopolymer, position2_a_homopolymer, trinucleotide_a_homopolymer, strand_a_homopolymer, chr2, pos1_2, pos2_2, strand_data2, mutation_a_homopolymer)
low_cutaneous_homopolymer = FindAsymmetry(low_homopolymer_file, c_homopolymer_chromosome, position1_c_homopolymer, position2_c_homopolymer, trinucleotide_c_homopolymer, strand_c_homopolymer, chr, pos1, pos2, strand_data, mutation_c_homopolymer)
high_cutaneous_homopolymer = FindAsymmetry(high_homopolymer_file, c_homopolymer_chromosome, position1_c_homopolymer, position2_c_homopolymer, trinucleotide_c_homopolymer, strand_c_homopolymer, chr2, pos1_2, pos2_2, strand_data2, mutation_c_homopolymer)

low_acral_homopolymer2 = FindAsymmetry(acral_low_homopolymer_file2, a_homopolymer_chromosome2, position1_a_homopolymer2, position2_a_homopolymer2, trinucleotide_a_homopolymer2, strand_a_homopolymer2, chr, pos1, pos2, strand_data, mutation_a_homopolymer2)
high_acral_homopolymer2 = FindAsymmetry(acral_high_homopolymer_file2,a_homopolymer_chromosome2, position1_a_homopolymer2, position2_a_homopolymer2, trinucleotide_a_homopolymer2, strand_a_homopolymer2, chr2, pos1_2, pos2_2, strand_data2, mutation_a_homopolymer2)
low_cutaneous_homopolymer2 = FindAsymmetry(low_homopolymer_file2, c_homopolymer_chromosome2, position1_c_homopolymer2, position2_c_homopolymer2, trinucleotide_c_homopolymer2, strand_c_homopolymer2, chr, pos1, pos2, strand_data, mutation_c_homopolymer2)
high_cutaneous_homopolymer2 = FindAsymmetry(high_homopolymer_file2, c_homopolymer_chromosome2, position1_c_homopolymer2, position2_c_homopolymer2, trinucleotide_c_homopolymer2, strand_c_homopolymer2, chr2, pos1_2, pos2_2, strand_data2, mutation_c_homopolymer2)



#low_acral_ins = FindAsymmetry(acral_low_file_ins, chromosome_ins, position1_ins, position2_ins, trinucleotide_ins, strand_ins, chr, pos1, pos2, strand_data, mutation)
#high_acral_ins = FindAsymmetry(acral_high_file_ins,chromosome_ins, position1_ins, position2_ins, trinucleotide_ins, strand_ins, chr2, pos1_2, pos2_2, strand_data2, mutation)
#low_cutaneous_ins = FindAsymmetry(low_file_ins, chromosome_cutaneous_ins, position1_cutaneous_ins, position2_cutaneous_ins, trinucleotide_cutaneous_ins, strand_cutaneous_ins, chr, pos1, pos2, strand_data, mutation_cutaneous)
#high_cutaneous_ins = FindAsymmetry(high_file_ins, chromosome_cutaneous_ins, position1_cutaneous_ins, position2_cutaneous_ins, trinucleotide_cutaneous_ins, strand_cutaneous_ins, chr2, pos1_2, pos2_2, strand_data2, mutation_cutaneous)




file1 = open("lowrpkm_melanocyte01_E059_gencodev10_hg19_sorted_name_strand_trinucfreq.txt")
file2 = open("lowrpkm_acral_strandasymmetry_deletions.txt", 'w+')
file1a = open("highrpkm_melanocyte01_E059_gencodev10_hg19_sorted_name_strand_trinucfreq.txt")
file2a = open("highrpkm_acral_strandasymmetry_deletions.txt", 'w+')
file3 = open("lowrpkm_cutaneous_strandasymmetry_deletions.txt", 'w+')
file3a = open("highrpkm_cutaneous_strandasymmetry_deletions.txt", 'w+')

def FindNormalizedRatios(file1, file2, NTS_count, TS_count, NTS_sequence_contexts, TS_sequence_contexts):
    Percentages = {}
    NTS_numbers = {}
    TS_numbers = {}
    Ratios = {}
    counts = file1.read()
    context_lines = [s.strip().split() for s in counts.splitlines()] 
    context_lines.pop(0)
    file2.write('NTS Count: ' + str(NTS_count))
    file2.write('\n')
    C_NTS = 0
    T_NTS = 0
    for key in NTS_sequence_contexts.keys():
        file2.write(key + ": " + str(NTS_sequence_contexts[key]))
        file2.write('\n')
        if key[1] == 'C':
              C_NTS = C_NTS + int(NTS_sequence_contexts[key])
        elif key[1] == 'T':
              T_NTS = T_NTS + int(NTS_sequence_contexts[key])
    file2.write('NTS C: ' + str(C_NTS))
    file2.write('\n')
    file2.write('NTS T: ' + str(T_NTS))
    file2.write('\n')

    file2.write('TS Count: ' + str(TS_count))
    file2.write("\n")
    C_TS = 0
    T_TS = 0
    for key in TS_sequence_contexts.keys():
        file2.write(key + ": " + str(TS_sequence_contexts[key]))
        file2.write('\n')
        if key[1] == 'C':
              C_TS = C_TS + int(TS_sequence_contexts[key])
        elif key[1] == 'T':
              T_TS = T_TS + int(TS_sequence_contexts[key])
    file2.write('TS C: ' + str(C_TS))
    file2.write('\n')
    file2.write('TS T: ' + str(T_TS))
    file2.write('\n')
    for line in context_lines:
         NTS_numbers[line[0]] = float(line[2])
         TS_numbers[line[0]] = float(line[1])
    NTS_total = 0
    TS_total = 0
    for key in NTS_numbers.keys():
         NTS_total = NTS_total + NTS_numbers[key]
    for key in TS_numbers.keys():
         TS_total = TS_total + TS_numbers[key]
    for key in NTS_numbers.keys():
        print(key)
        print(NTS_numbers[key]/NTS_total)
    for key in TS_numbers.keys():
        print(key)
        print(TS_numbers[key]/TS_total)

    for key in NTS_sequence_contexts.keys():
         if key in NTS_numbers.keys():
            NTS_sequence_contexts[key] = NTS_sequence_contexts[key]/NTS_numbers[key]
    for key in TS_sequence_contexts.keys():
         if key in TS_numbers.keys():
            TS_sequence_contexts[key] = TS_sequence_contexts[key]/TS_numbers[key]
    file2.write("Ratios:")
    file2.write('\n')
    for key in NTS_numbers.keys():
         if key in NTS_sequence_contexts.keys() and key in TS_sequence_contexts.keys():
            Ratios[key] = NTS_sequence_contexts[key]/TS_sequence_contexts[key]
            file2.write(key + ': ' + str(Ratios[key]))
            file2.write('\t')




FindNormalizedRatios(file1, file2, low_acral[0] + low_acral[2], low_acral[1] + low_acral[3], low_acral[4], low_acral[5])

FindNormalizedRatios(file1a, file2a, high_acral[0] + high_acral[2], high_acral[1] + high_acral[3], high_acral[4], high_acral[5])
file1.close()
file1 = open("lowrpkm_melanocyte01_E059_gencodev10_hg19_sorted_name_strand_trinucfreq.txt")
file1a.close()
file1a = open("highrpkm_melanocyte01_E059_gencodev10_hg19_sorted_name_strand_trinucfreq.txt")
FindNormalizedRatios(file1, file3, low_cutaneous[0] + low_cutaneous[2], low_cutaneous[1] + low_cutaneous[3], low_cutaneous[4], low_cutaneous[5])
FindNormalizedRatios(file1a, file3a, high_cutaneous[0] + high_cutaneous[2], high_cutaneous[1] + high_cutaneous[3], high_cutaneous[4], high_cutaneous[5])


file1.close()
file1a.close()
file2.close()
file2a.close()
file3.close()
file3a.close()
acral_low_file.close()
acral_high_file.close()
low_file.close()
high_file.close()
acral_low_homopolymer_file.close()
acral_high_homopolymer_file.close()
low_homopolymer_file.close()
high_homopolymer_file.close()
acral_low_homopolymer_file2.close()
acral_high_homopolymer_file2.close()
low_homopolymer_file2.close()
high_homopolymer_file2.close()
acral_low_file_ins.close()
acral_high_file_ins.close()
low_file_ins.close()
high_file_ins.close()

with open('high_expression.bed') as h:
     high_data = h.read()
with open('low_expression.bed') as l:
     low_data = l.read()
with open('acral_highexpression.bed') as ah:
    acral_high_data = ah.read()
with open('acral_lowexpression.bed') as al:
     acral_low_data = al.read()


with open('high_expression_ins.bed') as hi:
     high_data_ins = hi.read()
with open('low_expression_ins.bed') as li:
     low_data_ins = li.read()
with open('acral_highexpression_ins.bed') as ahi:
    acral_high_data_ins = ahi.read()
with open('acral_lowexpression_ins.bed') as ali:
     acral_low_data_ins = ali.read()

high_trinucleotide = getTrinucleotide(high_data)
low_trinucleotide = getTrinucleotide(low_data)
acral_high_trinucleotide = getTrinucleotide(acral_high_data)
acral_low_trinucleotide = getTrinucleotide(acral_low_data)



high_chromosome = getChromosome(high_data)
low_chromosome = getChromosome(low_data)
acral_high_chromosome = getChromosome(acral_high_data)
acral_low_chromosome = getChromosome(acral_low_data)

high_position = getPosition2(high_data)
low_position = getPosition2(low_data)
acral_high_position = getPosition2(acral_high_data)
acral_low_position = getPosition2(acral_low_data)

high_strand = getStrand(high_data)
low_strand = getStrand(low_data)
acral_high_strand = getStrand(acral_high_data)
acral_low_strand = getStrand(acral_low_data)


low_trinucleotide_ins = getTrinucleotide(low_data_ins)
acral_high_trinucleotide_ins = getTrinucleotide(acral_high_data_ins)
acral_low_trinucleotide_ins = getTrinucleotide(acral_low_data_ins)



high_chromosome_ins = getChromosome(high_data_ins)
low_chromosome_ins = getChromosome(low_data_ins)
acral_high_chromosome_ins = getChromosome(acral_high_data_ins)
acral_low_chromosome_ins = getChromosome(acral_low_data_ins)

high_position_ins = getPosition2(high_data_ins)
low_position_ins = getPosition2(low_data_ins)
acral_high_position_ins = getPosition2(acral_high_data_ins)
acral_low_position_ins = getPosition2(acral_low_data_ins)

high_strand_ins = getStrand(high_data_ins)
low_strand_ins = getStrand(low_data_ins)
acral_high_strand_ins = getStrand(acral_high_data_ins)
acral_low_strand_ins = getStrand(acral_low_data_ins)


def GetSignature(trinucleotide):
    signatures = {}
    for t in trinucleotide:
        if t not in signatures.keys():
            signatures[t] = 1
        else:
            signatures[t] = signatures[t] + 1
    return signatures

GetSignature(high_trinucleotide)
GetSignature(low_trinucleotide)
GetSignature(acral_high_trinucleotide)
GetSignature(acral_low_trinucleotide)


he = open('high_expression_deletion_sequences', 'w+')
le = open('low_exp_deletion_sequences', 'w+')
ahe = open('acral_high_expression_deletion_sequences', 'w+')
lea = open('acral_low_expression_deletion_sequences', 'w+')
hei = open('high_expression_insertion_sequences', 'w+')
lei = open('low_exp_insertion_sequences', 'w+')
ahei = open('acral_high_expression_insertion_sequences', 'w+')
leai = open('acral_low_expression_insertion_sequences', 'w+')

def MakeSequenceFile(file, chromosome, position, strand):
     for p in range(0, len(position)):
        if chromosome[p] == 'chrX':
            sequence = sequences[0]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr1':
            sequence = sequences[1]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr2':
            sequence = sequences[2]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr3':
            sequence = sequences[3]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr4':
            sequence = sequences[4]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr5':
            sequence = sequences[5]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr6':
            sequence = sequences[6]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr7':
            sequence = sequences[7]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr8':
            sequence = sequences[8]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr9':
            sequence = sequences[9]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr10':
            sequence = sequences[10]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr11':
            sequence = sequences[11]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr12':
            sequence = sequences[12]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr13':
            sequence = sequences[13]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr14':
            sequence = sequences[14]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr15':
            sequence = sequences[15]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr16':
            sequence = sequences[16]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr17':
            sequence = sequences[17]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr18':
            sequence = sequences[18]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr19':
            sequence = sequences[19]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr20':
            sequence = sequences[20]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr21':
            sequence = sequences[21]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
        elif chromosome[p] == 'chr22':
            sequence = sequences[22]
            seq = sequence[int(position[p]) - 11: int(position[p]) + 10]
            if strand[p] == '+':
                 file.write(seq)
                 file.write('\n')
            else:
                 file.write(GetComplement(seq))
                 file.write('\n')
          
     
MakeSequenceFile(he, high_chromosome, high_position, high_strand)
MakeSequenceFile(le, low_chromosome, low_position, low_strand)
MakeSequenceFile(ahe, acral_high_chromosome, acral_high_position, acral_high_strand)
MakeSequenceFile(lea, acral_low_chromosome, acral_low_position, acral_low_strand)

MakeSequenceFile(hei, high_chromosome_ins, high_position_ins, high_strand_ins)
MakeSequenceFile(lei, low_chromosome_ins, low_position_ins, low_strand_ins)
MakeSequenceFile(ahei, acral_high_chromosome_ins, acral_high_position_ins, acral_high_strand_ins)
MakeSequenceFile(leai, acral_low_chromosome_ins, acral_low_position_ins, acral_low_strand_ins)

he.close()
le.close()
ahe.close()
lea.close()

hei.close()
lei.close()
ahei.close()
leai.close()


file1 = open("lowrpkm_melanocyte01_E059_gencodev10_hg19_sorted_name_strand_trinucfreq.txt")
file1a = open("highrpkm_melanocyte01_E059_gencodev10_hg19_sorted_name_strand_trinucfreq.txt")
file = open("Iyer_saccer3_genecoord_noverlap_expand_strand_trinucfreq.txt")

def FindTC(file):
    AT_NTS = 0
    CG_NTS = 0
    AT_TS = 0
    CG_TS = 0
    counts = file.read()
    context_lines = [s.strip().split() for s in counts.splitlines()] 
    context_lines.pop(0)
    for line in context_lines:
        #for i in range(0, len(line[0])):
        if line[0][1] == 'T':
                 AT_TS = AT_TS + int(line[1])
                 AT_NTS = AT_NTS + int(line[2])
        elif  line[0][1] == 'C':
                 CG_TS = CG_TS + int(line[1])
                 CG_NTS = CG_NTS + int(line[2])
    return AT_NTS,CG_NTS, AT_TS,CG_TS

result = FindTC(file1)
result2 = FindTC(file1a)
FindTC(file)

T_NTS_base2 = result[0]
C_NTS_base2 = result[1]
T_TS_base2 = result[2]
C_TS_base2 = result[3]

T_NTS_base3 = result2[0]
C_NTS_base3 = result2[1]
T_TS_base3 = result2[2]
C_TS_base3 = result2[3]

file1.close()
file1a.close()
file.close()

file1 = open("lowrpkm_melanocyte01_E059_gencodev10_hg19_sorted_name_strand_trinucfreq.txt")
file1a = open("highrpkm_melanocyte01_E059_gencodev10_hg19_sorted_name_strand_trinucfreq.txt")
file = open("Iyer_saccer3_genecoord_noverlap_expand_strand_trinucfreq.txt")



def GetSequence(file):
    f1 = open(file)
    sequence1 = f1.read()
    sequence1 = str(sequence1.strip().upper())
    sequence1 = ''.join(sequence1.split('\n'))
    return sequence1


sequences2 = []
sequence1 = GetSequence('chr1.txt')
sequences2.append(sequence1)
print(len(sequence1))
sequence2 = GetSequence('chr2.txt')
sequences2.append(sequence2)
print(len(sequence2))
sequence3 = GetSequence('chr3.txt')
sequences2.append(sequence3)
print(len(sequence3))
sequence4 = GetSequence('chr4.txt')
sequences2.append(sequence4)
print(len(sequence4))
sequence5 = GetSequence('chr5.txt')
sequences2.append(sequence5)
print(len(sequence5))
sequence6 = GetSequence('chr6.txt')
print(len(sequence6))
sequences2.append(sequence6)
sequence7 = GetSequence('chr7.txt')
sequences2.append(sequence7)
print(len(sequence7))
sequence8 = GetSequence('chr8.txt')
sequences2.append(sequence8)
print(len(sequence8))
sequence9 = GetSequence('chr9.txt')
sequences2.append(sequence9)
print(len(sequence9))
sequence10 = GetSequence('chr10.txt')
sequences2.append(sequence10)
print(len(sequence10))
sequence11 = GetSequence('chr11.txt')
print(len(sequence11))
sequences2.append(sequence11)
sequence12 = GetSequence('chr12.txt')
print(len(sequence12))
sequences2.append(sequence12)
sequence13 = GetSequence('chr13.txt')
sequences2.append(sequence13)
print(len(sequence13))
sequence14 = GetSequence('chr14.txt')
sequences2.append(sequence14)
print(len(sequence14))
sequence15 = GetSequence('chr15.txt')
sequences2.append(sequence15)
print(len(sequence15))
sequence16 = GetSequence('chr16.txt')
sequences2.append(sequence16)
print(len(sequence16))



# with open('Iyer_saccer3_genecoord_noverlap.bed') as f2:
#     data3 = f2.read()
# lines3 = [s.strip().split() for s in data3.splitlines()]
# chr3 = []
# pos1_3 = []
# pos2_3 = []
# strand_data3 = []
# for line in lines3:
#     chr3.append(line[0])
#     pos1_3.append(int(line[1]))
#     pos2_3.append(int(line[2]))
#     strand_data3.append(line[5])



# TS_base = 0
# NTS_base = 0
# for x in range (0, len(chr3)):
#     if chr3[x] == 'chrI':
#         sequence = sequence1
#     if chr3[x] == 'chrII':
#         sequence = sequence2
#     if chr3[x] == 'chrIII':
#         sequence = sequence3
#     if chr3[x] == 'chrIV':
#         sequence = sequence4
#     if chr3[x] == 'chrV':
#         sequence = sequence5
#     if chr3[x] == 'chrVI':
#         sequence = sequence6
#     if chr3[x] == 'chrVII':
#         sequence = sequence7
#     if chr3[x] == 'chrVIII':
#         sequence = sequence8
#     if chr3[x] == 'chrIX':
#         sequence = sequence9
#     if chr3[x] == 'chrX':
#         sequence = sequence10
#     if chr3[x] == 'chrXI':
#         sequence = sequence11
#     if chr3[x] == 'chrXII':
#         sequence = sequence12
#     if chr3[x] == 'chrXIII':
#         sequence = sequence13
#     if chr3[x] == 'chrXIV':
#         sequence = sequence14
#     if chr3[x] == 'chrXV':
#         sequence = sequence15
#     if chr3[x] == 'chrXVI':
#         sequence = sequence16
#     #for i in range(0, (pos2_3[x] - pos1_3[x]) + 1):
#     b = pos1_3[x]
#     while b <= (pos2_3[x]) - 1:
#         #for b in range(0, len(sequence) - 1):
#             polymer_range = {}
#             length = 1
#             if sequence[b] != 'N':# and b not in coordinates:
#                 length = 1
#                 i = 1
#                 polymer_range[0] = b
#                 while sequence[b] == sequence[b + i]:
#                     if (b + i) <= (pos2_3[x]) - 1:
#                         polymer_range[1] = b + i
#                         length = length + 1
#                         i = i + 1
#                     else:
#                         break
                   
              

#                 if length >= 4:
                
#                     if strand_data3[x] == '+':
#                                 if sequence[b] == 'T' or sequence[b] == 'C':
#                                         NTS_base = NTS_base + length
                                        
#                                 if sequence[b] == 'A' or sequence[b] == 'G':
#                                         TS_base = TS_base + length
#                     elif strand_data3[x] == '-':
#                                 if sequence[b] == 'A' or sequence[b] == 'G':
#                                         NTS_base = NTS_base + length
#                                 if sequence[b] == 'T' or sequence[b] == 'C':
#                                         TS_base = TS_base + length
#                 b = b + length
# print(NTS_base)
# print(TS_base)


# with open('highrpkm_melanocyte01_E059_gencodev10_hg19_sorted_name_noverlap.bed') as f3:
#     data4 = f3.read()
# lines4 = [s.strip().split() for s in data4.splitlines()]
# chr4 = []
# pos1_4 = []
# pos2_4 = []
# strand_data4 = []
# for line in lines4:
#     chr4.append(line[0])
#     pos1_4.append(int(line[1]))
#     pos2_4.append(int(line[2]))
#     strand_data4.append(line[5])



# TS_base2 = 0
# NTS_base2 = 0
# for x in range (0, len(chr4)):
#     if chr4[x] == 'chrX':
#         sequence = sequences[0]
#     if chr4[x] == 'chr1':
#         sequence = sequences[1]
#     if chr4[x] == 'chr2':
#         sequence = sequences[2]
#     if chr4[x] == 'chr3':
#         sequence = sequences[3]
#     if chr4[x] == 'chr4':
#         sequence = sequences[4]
#     if chr4[x] == 'chr5':
#         sequence = sequences[5]
#     if chr4[x] == 'chr6':
#         sequence = sequences[6]
#     if chr4[x] == 'chr7':
#         sequence = sequences[7]
#     if chr4[x] == 'chr8':
#         sequence = sequences[8]
#     if chr4[x] == 'chr9':
#         sequence = sequences[9]
#     if chr4[x] == 'chr10':
#         sequence = sequences[10]
#     if chr4[x] == 'chr11':
#         sequence = sequences[11]
#     if chr4[x] == 'chr12':
#         sequence = sequences[12]
#     if chr4[x] == 'chr13':
#         sequence = sequences[13]
#     if chr4[x] == 'chr14':
#         sequence = sequences[14]
#     if chr4[x] == 'chr15':
#         sequence = sequences[15]
#     if chr4[x] == 'chr16':
#         sequence = sequences[16]
#     if chr4[x] == 'chr17':
#         sequence = sequences[17]
#     if chr4[x] == 'chr18':
#         sequence = sequences[18]
#     if chr4[x] == 'chr19':
#         sequence = sequences[19]
#     if chr4[x] == 'chr20':
#         sequence = sequences[20]
#     if chr4[x] == 'chr21':
#         sequence = sequences[21]
#     if chr4[x] == 'chr22':
#         sequence = sequences[22]
#     if chr4[x] == 'chrY':
#         sequence = sequences[23]
#     #for i in range(0, (pos2_3[x] - pos1_3[x]) + 1):
#     b = pos1_4[x]
#     while b <= (pos2_4[x]) - 1:
#         #for b in range(0, len(sequence) - 1):
#             polymer_range = {}
#             length = 1
#             if sequence[b] != 'N':# and b not in coordinates:
#                 length = 1
#                 i = 1
#                 polymer_range[0] = b
               
#                 while sequence[b] == sequence[b + i]:
#                     if (b + i) < (pos2_4[x]) - 1:
#                         polymer_range[1] = b + i
#                         length = length + 1
#                         i = i + 1
#                     else:
#                         break
                  

#                 if length >= 4:
                
#                     if strand_data4[x] == '+':
#                             if sequence[b] == 'T' or sequence[b] == 'C':
#                                     NTS_base2 = NTS_base2 + length
                                    
#                             if sequence[b] == 'A' or sequence[b] == 'G':
#                                     TS_base2 = TS_base2 + length
#                     elif strand_data4[x] == '-':
#                             if sequence[b] == 'A' or sequence[b] == 'G':
#                                     NTS_base2 = NTS_base2 + length
#                             if sequence[b] == 'T' or sequence[b] == 'C':
#                                     TS_base2 = TS_base2 + length
            
#                 b = b + length
#             else:
#                  b = b + 1
# print(NTS_base2)
# print(TS_base2)
                    
        
# with open('lowrpkm_melanocyte01_E059_gencodev10_hg19_sorted.bed') as f4:
#     data5 = f4.read()
# lines5 = [s.strip().split() for s in data5.splitlines()]
# chr5 = []
# pos1_5 = []
# pos2_5 = []
# previous_plus2 = []
# previous_minus2 = []
# strand_data5 = []
# for line in lines5:
#     chr5.append(line[0])
#     pos1_5.append(int(line[1]))
#     pos2_5.append(int(line[2]))
#     strand_data5.append(line[5])



# TS_base3 = 0
# NTS_base3 = 0
# for x in range (0, len(chr5)):
#     if chr5[x] == 'chrX':
#         sequence = sequences[0]
#     if chr5[x] == 'chr1':
#         sequence = sequences[1]
#     if chr5[x] == 'chr2':
#         sequence = sequences[2]
#     if chr5[x] == 'chr3':
#         sequence = sequences[3]
#     if chr5[x] == 'chr4':
#         sequence = sequences[4]
#     if chr5[x] == 'chr5':
#         sequence = sequences[5]
#     if chr5[x] == 'chr6':
#         sequence = sequences[6]
#     if chr5[x] == 'chr7':
#         sequence = sequences[7]
#     if chr5[x] == 'chr8':
#         sequence = sequences[8]
#     if chr5[x] == 'chr9':
#         sequence = sequences[9]
#     if chr5[x] == 'chr10':
#         sequence = sequences[10]
#     if chr5[x] == 'chr11':
#         sequence = sequences[11]
#     if chr5[x] == 'chr12':
#         sequence = sequences[12]
#     if chr5[x] == 'chr13':
#         sequence = sequences[13]
#     if chr5[x] == 'chr14':
#         sequence = sequences[14]
#     if chr5[x] == 'chr15':
#         sequence = sequences[15]
#     if chr5[x] == 'chr16':
#         sequence = sequences[16]
#     if chr5[x] == 'chr17':
#         sequence = sequences[17]
#     if chr5[x] == 'chr18':
#         sequence = sequences[18]
#     if chr5[x] == 'chr19':
#         sequence = sequences[19]
#     if chr5[x] == 'chr20':
#         sequence = sequences[20]
#     if chr5[x] == 'chr21':
#         sequence = sequences[21]
#     if chr5[x] == 'chr22':
#         sequence = sequences[22]
#     if chr5[x] == 'chrY':
#         sequence = sequences[23]
#     #for i in range(0, (pos2_3[x] - pos1_3[x]) + 1):
#     b = pos1_5[x]
#     while b <= (pos2_5[x]) - 1:
#         #for b in range(0, len(sequence) - 1):
#             polymer_range = {}
#             length = 1
#             if sequence[b] != 'N':# and b not in coordinates:
#                 if (chr5[x], b) not in previous_plus2 and (chr5[x], b) not in previous_minus2:

#                     length = 1
#                     i = 1
#                     polymer_range[0] = b
#                     while sequence[b] == sequence[b + i]:
#                         if (b + i) <= (pos2_5[x]) - 1:
#                             polymer_range[1] = b + i
#                             length = length + 1
#                             i = i + 1
#                         else:
#                             break
              

#                 #if length >= 4:
#                 #if (chr5[x], b) not in previous_plus2 and (chr5[x], b) not in previous_minus2:

#                     if strand_data5[x] == '+':
#                             if sequence[b] == 'T' or sequence[b] == 'C':
#                                     NTS_base3 = NTS_base3 + length
                                    
#                             if sequence[b] == 'A' or sequence[b] == 'G':
#                                     TS_base3 = TS_base3 + length
#                             previous_plus2.append((chr5[x], b))
#                     elif strand_data5[x] == '-':
#                             if sequence[b] == 'A' or sequence[b] == 'G':
#                                     NTS_base3 = NTS_base3 + length
#                             if sequence[b] == 'T' or sequence[b] == 'C':
#                                     TS_base3 = TS_base3 + length
#                             previous_minus2.append((chr5[x], b))
                
#                     b = b + length

#                 else:
#                     b = b + 1
                 
            
#             else:
#                  b = b + 1
# print(NTS_base3)
# print(TS_base3)

#Caculates normalized percentages on the NTS and TS based on sequence context
def FindNormalizedRatio(NTS_count, TS_count, NTS_numbers, TS_numbers):
    base_ratio = TS_numbers/NTS_numbers
    mutation_ratio = NTS_count/TS_count
    ratio = base_ratio * mutation_ratio
    NTS_percentage = (ratio * 100)/ (ratio + 1)
    TS_percentage = 100/ (ratio + 1)
    return NTS_percentage, TS_percentage



FindNormalizedRatio(low_cutaneous[0], low_cutaneous[2], T_NTS_base2, T_TS_base2)
FindNormalizedRatio(low_cutaneous[1], low_cutaneous[3], C_NTS_base2, C_TS_base2)
FindNormalizedRatio (high_cutaneous[0], high_cutaneous[2], T_NTS_base3, T_TS_base3,)
FindNormalizedRatio( high_cutaneous[1], high_cutaneous[3], C_NTS_base3, C_TS_base3,)