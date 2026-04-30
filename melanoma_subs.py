def GetReverseComplement(str):
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


#Returns the start and end positions and the length of the homopolymer containing a base substitution
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
        
      

#Finds and extracts the sequence of a chromosome from the human genome hg19 fasta file
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
sequenceY = GetSequence('hg19.fa', 'chrY')
sequences.append(sequenceY)


with open("Cutaneous_somatic_mutation_MELA_AU_SNVs_sorted_trinuc.bed") as f:  
    data = f.read()

lines = [s.strip().split() for s in data.splitlines()]
chromosome = []
position1 = []
position2 = []
allele = []
mutation = []
length = []
isolate = []
trinucleotide = []
strand = []
for line in lines:
    chromosome.append(line[0])
    position1.append(int(line[1]) + 1)
    position2.append(int(line[2]) - 1)
    #allele.append(line[3])
    mutation.append(line[4][2])
    strand.append(line[5])


for pos in range(0, len(position2)):
        if chromosome[pos] == 'chrX':
            sequence = sequenceX
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])

        elif chromosome[pos] == 'chrY':
            sequence = sequenceY
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])
        
        
            
        elif chromosome[pos] == 'chr1':
            sequence = sequence1
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])
        
        elif chromosome[pos] == 'chr2':
            sequence = sequence2
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])
        
        
        elif chromosome[pos] == 'chr3':
            sequence = sequence3
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)

            allele.append(trinucleotide[len(trinucleotide) - 1][1])

        elif chromosome[pos] == 'chr4':
            sequence = sequence4
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])


        elif chromosome[pos] == 'chr5':
            sequence = sequence5
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])


        elif chromosome[pos] == 'chr6':
            sequence = sequence6
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])


        elif chromosome[pos] == 'chr7':
            sequence = sequence7
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])


        elif chromosome[pos] == 'chr8':
            sequence = sequence8
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])


        elif chromosome[pos] == 'chr9':
            sequence = sequence9
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])


        elif chromosome[pos] == 'chr10':
            sequence = sequence10
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])
            
        elif chromosome[pos] == 'chr11':
            sequence = sequence11
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])


        elif chromosome[pos] == 'chr12':
            sequence = sequence12
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])

        elif chromosome[pos] == 'chr13':
            sequence = sequence13
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])

        
        elif chromosome[pos] == 'chr14':
            sequence = sequence14
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])

        elif chromosome[pos] == 'chr15':
            sequence = sequence15
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])

        elif chromosome[pos] == 'chr16':
            sequence = sequence16
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])

        elif chromosome[pos] == 'chr17':
            sequence = sequence17
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])

        elif chromosome[pos] == 'chr18':
            sequence = sequence18
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])

        elif chromosome[pos] == 'chr19':
            sequence = sequence19
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])

        elif chromosome[pos] == 'chr20':
            sequence = sequence20
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])

        elif chromosome[pos] == 'chr21':
            sequence = sequence21
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])

        elif chromosome[pos] == 'chr22':
            sequence = sequence22
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]

            if strand[pos] == '-':
                if '-' not in tri and tri.isalpha() == False:
                    print(tri)
                trinucleotide.append(GetReverseComplement(tri))
                
            else:
                trinucleotide.append(tri)
            allele.append(trinucleotide[len(trinucleotide) - 1][1])
        
        else:
            trinucleotide.append('NA')
            strand.append('+')
            allele.append('NA')
    

file_1 = open('cutaneous_substitutions.bed', 'w+')
file_2 = open('cutaneous_substitution_sequences', 'w+')


#Makes bed files
def MakeMultiNucFile(file_1, file_2):
    new_position = position2
    new_chromosome = chromosome
    new_strand = strand
    new_position1 = position1
    new_trinucleotide = trinucleotide
    new_mutation = mutation
    new_allele = allele

  #Put indel in middle
    multinucleotide = []
    lengths = []
    for pos in range(0, len(new_position)):
        times = 0
        if new_chromosome[pos] == 'chrX':
            sequence = sequenceX
            #if strand[pos] == '+':

            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
                
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1

        elif new_chromosome[pos] == 'chrY':
            sequence = sequenceY
            #if strand[pos] == '+':

            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
                
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        
        
        
        
        elif new_chromosome[pos] == 'chr1':
            sequence = sequence1
            #if strand[pos] == '+':
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
                
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr2':
            sequence = sequence2
            #if strand[pos] == '+':
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr3':
            sequence = sequence3
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr4':
            sequence = sequence4
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr5':
            sequence = sequence5
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr6':
            sequence = sequence6
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr7':
            sequence = sequence7
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr8':
            sequence = sequence8
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr9':
            sequence = sequence9
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr10':
            sequence = sequence10
            #if strand[pos] == '+':
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr11':
            sequence = sequence11
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr12':
            sequence = sequence12
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr13':
            sequence = sequence13
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr14':
            sequence = sequence14
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr15':
            sequence = sequence15
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr16':
            sequence = sequence16
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr17':
            sequence = sequence17
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr18':
            sequence = sequence18
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr19':
            sequence = sequence19
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr20':
            sequence = sequence20
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        elif new_chromosome[pos] == 'chr21':
            sequence = sequence21
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1

        elif new_chromosome[pos] == 'chr22':
            sequence = sequence22
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                lengths.append(GetHomopolymerLength(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                times = times + 1
            elif int(new_position[pos]) <= 27:
                multinucleotide_seq = ""
                for i in range(0, 27- int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'
                    multinucleotide_seq = multinucleotide_seq + sequence[0: int(new_position[pos]) + 27]
                if new_strand[pos] == '+':
                        multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
           
                lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                times = times + 1
        else:
            multinucleotide.append('NA')
            lengths.append('NA')
            times = times + 1
        if times > 1:
            print(times)
            
            
    
    for i in range(0, len(new_position)):
        #if lengths[i] != 'NA' and int(lengths[i]) < 4:
        if lengths[i] != 'NA':
       
                file_1.write(new_chromosome[i] + '\t' + str(int(new_position1[i]))  + '\t' + str(new_position[i]) + '\t' + new_trinucleotide[i] + '\t' +  new_mutation[i] + '\t' + new_strand[i] + '\t' + multinucleotide[i] + '\t' + str(GetHomopolymerLength(multinucleotide[i], new_position1[i], new_position[i], len(multinucleotide[i])/2 + 1)[0]) + '\t' + str(lengths[i]))
                file_1.write('\n')
                if new_strand[i] == '+':
                    if multinucleotide[i] != 'NA':
                        file_2.write(multinucleotide[i])
                        file_2.write('\n')

                if new_strand[i] == '-':
                    if multinucleotide[i] != 'NA':
                            file_2.write(multinucleotide[i][::-1])
                            file_2.write('\n')
            
        

    return multinucleotide, lengths

lengths = MakeMultiNucFile(file_1, file_2)[1]

def MutationsPerTumor(allele, tumor_number):
   
    substitution_number = 0
  
    for i in range(0, len(allele)):
        if len(allele[i]) == 1:
            substitution_number = substitution_number + 1
        
            
    print(substitution_number / tumor_number)
  
    return substitution_number

MutationsPerTumor(allele, 140)

def HomopolymerMutationsPerTumor(lengths, tumor_number, allele):
    nonhomopolymer = 0
    homopolymer = 0
    for i in range(0, len(lengths)):
        if lengths[i] != 'NA':
            if lengths[i] > 3:
                homopolymer = homopolymer + 1
                #if len(allele[i]) == 1:
                    #homopolymer = homopolymer + 1
            else:
                nonhomopolymer = nonhomopolymer + 1
    print(nonhomopolymer / tumor_number)
    print(homopolymer / tumor_number)

    return lengths


HomopolymerMutationsPerTumor(lengths, 140, allele)


homopolymer_file = open("hg_homopolymers", 'w+')

#Obtain the total number of bases in the human genome in sequences of each homopolymer length
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


homopolymer_dict = NormalizeHomopolymer(sequences, homopolymer_file)
new_homopolymer_dict = homopolymer_dict[0]
C_dict = homopolymer_dict[1]
T_dict = homopolymer_dict[2]



file_1.close()
file_2.close()
homopolymer_file.close()

file_3 = open('melanoma_sub_homopolymers', 'w+')

def PrintFrequencies(homopolymer_dict, C_dict, T_dict, lengths, allele, mutation, file, tumor_number ):
    homopolymer_counts = {}
    C_homopolymer_counts = {}
    T_homopolymer_counts = {}
    homopolymer_counts['7+'] = 0
    homopolymer_counts['21+'] = 0
    homopolymer_counts['16-20'] = 0
    homopolymer_counts['11-15'] = 0
    C_homopolymer_counts['7+'] = 0
    C_homopolymer_counts['21+'] = 0
    C_homopolymer_counts['16-20'] = 0
    C_homopolymer_counts['11-15'] = 0
    T_homopolymer_counts['7+'] = 0
    T_homopolymer_counts['21+'] = 0
    T_homopolymer_counts['16-20'] = 0
    T_homopolymer_counts['11-15'] = 0
    for l in range(0, len(lengths)):
        if allele[l] != 'NA' and lengths[l] != 'NA':
            if len(allele[l]) == 1:
                if lengths[l] not in homopolymer_counts.keys():
                    homopolymer_counts[lengths[l]] = 1
                    #homopolymer_dict[lengths[l]] = 0
                else:
                    homopolymer_counts[lengths[l]] = homopolymer_counts[lengths[l]] + 1
                if int(lengths[l]) >= 7:
                    homopolymer_counts['7+'] = homopolymer_counts['7+'] + 1
                if int(lengths[l]) >= 11:
                    if int(lengths[l]) >= 21:
                        homopolymer_counts['21+'] = homopolymer_counts['21+'] + 1
                    elif int(lengths[l]) >= 16:
                            homopolymer_counts['16-20'] = homopolymer_counts['16-20'] + 1
                    else:
                        homopolymer_counts['11-15'] = homopolymer_counts['11-15'] + 1
                
                if allele[l] == 'C' or allele[l] == 'G'or mutation[l] == 'C' or mutation[l] == 'G':
                    if lengths[l] not in C_homopolymer_counts.keys():
                        C_homopolymer_counts[lengths[l]] = 1
                        #homopolymer_dict[lengths[l]] = 0
                    else:
                        C_homopolymer_counts[lengths[l]] = C_homopolymer_counts[lengths[l]] + 1
    
                    if int(lengths[l]) >= 7:
                            C_homopolymer_counts['7+'] = C_homopolymer_counts['7+'] + 1
                    if int(lengths[l]) >= 11:
                        if int(lengths[l]) >= 21:
                            C_homopolymer_counts['21+'] = C_homopolymer_counts['21+'] + 1
                        elif int(lengths[l]) >= 16:
                            C_homopolymer_counts['16-20'] = C_homopolymer_counts['16-20'] + 1
                        else:
                            C_homopolymer_counts['11-15'] = C_homopolymer_counts['11-15'] + 1
                    
                if allele[l] == 'T' or allele[l] == 'A'or mutation[l] == 'T' or mutation[l] == 'A':
                    if lengths[l] not in T_homopolymer_counts.keys():
                        T_homopolymer_counts[lengths[l]] = 1
                        #homopolymer_dict[lengths[l]] = 0
                    else:
                        T_homopolymer_counts[lengths[l]] = T_homopolymer_counts[lengths[l]] + 1
    
                    if int(lengths[l]) >= 7:
                            T_homopolymer_counts['7+'] = T_homopolymer_counts['7+'] + 1
                    if int(lengths[l]) >= 11:
                        if int(lengths[l]) >= 21:
                            T_homopolymer_counts['21+'] = T_homopolymer_counts['21+'] + 1
                        elif int(lengths[l]) >= 16:
                            T_homopolymer_counts['16-20'] = T_homopolymer_counts['16-20'] + 1
                        else:
                            T_homopolymer_counts['11-15'] = T_homopolymer_counts['11-15'] + 1
    for key in homopolymer_counts.keys():
        if homopolymer_dict[key] != 0:
            homopolymer_counts[key] = homopolymer_counts[key]/homopolymer_dict[key]
            file.write(str(key) + ':' + str(homopolymer_counts[key]/tumor_number))
            file.write('\n')
    for key in C_homopolymer_counts.keys():
        if C_dict[key] != 0:
            C_homopolymer_counts[key] = C_homopolymer_counts[key]/C_dict[key]
            file.write(str(key) + ':' + str(C_homopolymer_counts[key]/tumor_number))
            file.write('\n')
    for key in T_homopolymer_counts.keys():
        if T_dict[key] != 0:
            T_homopolymer_counts[key] = T_homopolymer_counts[key]/T_dict[key]
            file.write(str(key) + ':' + str(T_homopolymer_counts[key]/ tumor_number))
            file.write('\n')


PrintFrequencies(homopolymer_dict, C_dict, T_dict, lengths, allele, mutation, file_3, 140 )
file_3.close()