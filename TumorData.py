

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


#Returns the start and end positions and the length of the homopolymer containing a deleted base
def GetHomopolymerLength(sequence, mut_pos1, mut_pos, mut_position):
    #mut_position = int(len(sequence)/2) + 1
    mut_pos1 = int(mut_pos1)
    mut_pos = int(mut_pos) 
    mut_position = int(mut_position)
    difference = mut_pos - mut_pos1 - 1
    mut_position2 = mut_position - difference
    print(sequence[mut_position - 1])
    
    if difference > 1:
        print(sequence[mut_position2 - 1:mut_position])
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
               
                if mut_position <= end_position  and mut_position >= end_position - (length - 1 ):
                     # or [(mut_position2 <= end_position)  and (mut_position2 >= end_position - (length - 1 ))]:
                     lengths[(end_position  + (mut_pos - mut_position) - (length - 1 ) , end_position  + (mut_pos - mut_position))] = length
                elif mut_position2 <= end_position  and mut_position2 >= end_position - (length - 1 ):
                    lengths[(end_position  + (mut_pos - mut_position) - (length - 1 ) , end_position  + (mut_pos - mut_position))] = length
    
        if mut_position <= end_position  and mut_position >= end_position - (length - 1 ):
        
    
            lengths[(end_position  + (mut_pos - mut_position) - (length - 1 ) , end_position  + (mut_pos - mut_position) )] = length
        elif mut_position2 <= end_position  and mut_position2 >= end_position - (length - 1 ):
                    lengths[(end_position  + (mut_pos - mut_position) - (length - 1 ) , end_position  + (mut_pos - mut_position))] = length   
    if len(lengths) == 0:
        return 1, (mut_pos, mut_pos)
    for key in lengths:
        if lengths[key] == max(lengths.values()):
            return lengths[key], key

#Returns the start and end positions and the length of the homopolymer containing a inserted base
def GetHomopolymerLengthIns(sequence, mut_pos1, mut_pos, mut_position):
    #mut_position = int(len(sequence)/2) + 1
    mut_position = mut_position 
    mut_pos1 = int(mut_pos1)
    mut_pos = int(mut_pos) 
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
                #if strand == '-':
                    #if mut_position <= end_position  and mut_position >= end_position - (length - 1 ) or (mut_position + 1) <= end_position  and (mut_position + 1) >= end_position - (length - 1 ):
                        #lengths[(end_position  + (mut_pos - mut_position) - (length - 1 ) , end_position  + (mut_pos - mut_position))] = length
                #else:
                if mut_position <= end_position + 1  and (mut_position) >= end_position - (length - 1):
                        lengths[(end_position  + (mut_pos - mut_position) - (length - 1 ) , end_position  + (mut_pos - mut_position))] = length
              
    
       
            if mut_position <= end_position + 1   and mut_position >= end_position - (length - 1):
                lengths[(end_position  + (mut_pos - mut_position) - (length - 1 ) , end_position  + (mut_pos - mut_position))] = length
         
    
           
    if len(lengths) == 0:
        return 1, (mut_pos, mut_pos)
    for key in lengths:
        if lengths[key] == max(lengths.values()):
            return lengths[key], key
      

#Finds and extracts the sequence of a chromosome from the human genome fasta file
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



acral = "Acral_simple_somatic_mutation_MELA_AU_DIPs.bed"
cutaneous = "Cutaneous_simple_somatic_mutation_MELA_AU_DIPs.bed"



def InputData(filename):
    with open(filename) as f:
        data = f.read()
    lines = [s.strip().split() for s in data.splitlines()]

    chromosome = []
    position1 = []
    position2 = []
    mutation_type = []
    allele = []
    mutation = []
    trinucleotide = []
    multiNuc = []
    strand = []
    multid = 0
    multii = 0
    deletions = 0
    insertions = 0

    for line in lines:
            chromosome.append(line[0])
            position1.append(line[1])
            position2.append(line[2])
            if 'deletion' in line[3]:
                if len(line[6]) > 1:
                    multid = multid + 1
                    multiNuc.append('Y')
                else:
                    deletions = deletions + 1
                    multiNuc.append('N')
                mutation_type.append('deletion')
        
            elif 'insertion' in line[3]:
                if len(line[7]) > 1:
                    multii = multii + 1
                    multiNuc.append('Y')
                else:
                    insertions = insertions + 1
                    multiNuc.append('N')
                mutation_type.append('insertion')
            
            else:
                mutation_type.append('NA')
            allele.append(line[6])
            mutation.append(line[7])
            
        

    for pos in range(0, len(position2)):
        if chromosome[pos] == 'chrX':
            sequence = sequenceX
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chrY':
            sequence = sequenceY
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr1':
            sequence = sequence1
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr2':
            sequence = sequence2
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr3':
            sequence = sequence3
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr4':
            sequence = sequence4
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr5':
            sequence = sequence5
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr6':
            sequence = sequence6
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr7':
            sequence = sequence7
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr8':
            sequence = sequence8
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr9':
            sequence = sequence9
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr10':
            sequence = sequence10
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr11':
            sequence = sequence11
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr12':
            sequence = sequence12
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr13':
            sequence = sequence13
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr14':
            sequence = sequence14
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr15':
            sequence = sequence15
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr16':
            sequence = sequence16
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr17':
            sequence = sequence17
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr18':
            sequence = sequence18
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr19':
            sequence = sequence19
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr20':
            sequence = sequence20
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr21':
            sequence = sequence21
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        elif chromosome[pos] == 'chr22':
            sequence = sequence22
            tri = sequence[int(position2[pos]) - 2: int(position2[pos]) + 1]
            if mutation[pos] == '-':
                if tri[1] == 'A' or tri[1] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
            else:
                if mutation[pos] == 'A' or mutation[pos] == 'G':
                    trinucleotide.append(GetReverseComplement(tri))
                    strand.append('-')
                else:
                    trinucleotide.append(tri)
                    strand.append('+')
        
        else:
            trinucleotide.append('NA')
            strand.append('+')
    #Counts the number of deletions and insertions of different lengths
    insertions = {}
    deletions = {}
    for m in range(0, len(mutation_type)):
        if mutation_type[m] == 'insertion':
            if len(mutation[m]) in insertions.keys():
                    insertions[len(mutation[m])] = insertions[len(mutation[m])] + 1
            else: 
                    insertions[len(mutation[m])] = 1
        elif mutation_type[m] == 'deletion':
            if len(allele[m]) in deletions.keys():
                deletions[len(allele[m])] = deletions[len(allele[m])] + 1
            else: 
                deletions[len(allele[m])] = 1
    

    #else:
        #trinucleotide.append('NA')
    return [chromosome, position1, position2, mutation_type, allele, mutation, trinucleotide, strand, multiNuc]


#InputData(acral)
#InputData(cutaneous)




def MakeMultiNucFile(file_1, file_2, file_3, file_4, file_5, file_6, file_7, file_8, file_9): #Finds the 27 bases on either side of a mutation and the length of the homopolymer containing the mutation
    new_data = InputData(file_1)
    new_chromosome = new_data[0]
    new_position1 = new_data[1]
    new_position = new_data[2]
    new_mutation_type = new_data[3]
    new_allele = new_data[4]
    new_mutation = new_data[5]
    new_trinucleotide = new_data[6]
    new_strand = new_data[7]
    new_multiNuc = new_data[8]
   

  #Put indel in middle
    multinucleotide = []
    lengths = []
    for pos in range(0, len(new_position)):
        if new_chromosome[pos] == 'chrX':
            sequence = sequenceX
            #if strand[pos] == '+':
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        
        
        elif new_chromosome[pos] == 'chrY':
            sequence = sequenceY
            #if strand[pos] == '+':
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])  
        
        
        elif new_chromosome[pos] == 'chr1':
            sequence = sequence1
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chr2':
            sequence = sequence2
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chr3':
            sequence = sequence3
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chr4':
            sequence = sequence4
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chr5':
            sequence = sequence5
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chr6':
            sequence = sequence6
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chr7':
            sequence = sequence7
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chr8':
            sequence = sequence8
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chr9':
            sequence = sequence9
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chr10':
            sequence = sequence10
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chr11':
            sequence = sequence11
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chr12':
            sequence = sequence12
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chr13':
            sequence = sequence13
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chr14':
            sequence = sequence14
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chr15':
            sequence = sequence15
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chr16':
            sequence = sequence16
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chr17':
            sequence = sequence17
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])

        elif new_chromosome[pos] == 'chr18':
            sequence = sequence18
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chr19':
            sequence = sequence19
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chr20':
            sequence = sequence20
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        elif new_chromosome[pos] == 'chr21':
            sequence = sequence21
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])

        elif new_chromosome[pos] == 'chr22':
            sequence = sequence22
            if int(new_position[pos]) > 27 and int(new_position[pos]) < len(sequence) -26:
                if new_strand[pos] == '+':
                    multinucleotide.append(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])
                else:
                    multinucleotide.append(GetReverseComplement(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27]))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[0])
                    if lengths[len(lengths) - 1] >= 4:
                        h_coord = GetHomopolymerLengthIns(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27], int(new_position1[pos]), int(new_position[pos]), int(len(sequence[int(new_position[pos]) - 28: int(new_position[pos]) + 27])/2) + 1)[1]
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+' #Insertions in C and T homopolymers are assumed to be on the '+' strand
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
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
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
            else:#elif int(new_position[pos]) >= len(sequence) -26:
                multinucleotide_seq = sequence[int(new_position[pos]) - 28: len(sequence)]
                for i in range(0, len(sequence) - 26 - int(new_position[pos]) + 1):
                    multinucleotide_seq = multinucleotide_seq + '*'

                if new_strand[pos] == '+':
                    multinucleotide.append(multinucleotide_seq)
                else:
                     multinucleotide.append(GetReverseComplement(multinucleotide_seq))
                if new_mutation_type[pos] == 'insertion':
                    lengths.append(GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
                    h_coord = GetHomopolymerLengthIns(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[1]
                    if lengths[len(lengths) - 1] >= 4:
                        if 'T' in sequence[h_coord[0] - 1: h_coord[1]] or 'C' in sequence[h_coord[0] - 1: h_coord[1]]:
                            new_strand[pos] == '+'
                        else:
                            new_strand[pos] = '-'
                elif new_mutation_type[pos] == 'deletion':
                    lengths.append(GetHomopolymerLength(multinucleotide_seq, int(new_position1[pos]), int(new_position[pos]), int(len(multinucleotide_seq)/2) + 1)[0])
        else:
            multinucleotide.append('NA')
            lengths.append('NA')
            
            
    
    for i in range(0, len(new_position)):
        #if lengths[i] != 'NA' and int(lengths[i]) < 4:
        if lengths[i] != 'NA' and (new_allele[i] == 'C' or new_allele[i] == 'G'):
        #if lengths[i] != 'NA' and (new_allele[i] == 'T' or new_allele[i] == 'A'):
        #if new_strand[i] == '-':
            #file_2.write(new_chromosome[i] + '\t' + str(int(new_position[i]) - 1) + '\t'  + new_position[i]  + '\t' + new_trinucleotide[i] + '\t' + new_mutation[i] + '\t' + new_strand[i]  + '\t'+ GetReverseComplement(multinucleotide[i]) + '\t' +str(GetHomopolymerLength(GetReverseComplement(multinucleotide[i]), new_position[i])[0]) + '\t' + str(lengths[i]))
            #file_2.write('\n')
        #if new_strand[i] == '+':
            if new_mutation_type[i] == "deletion" and new_multiNuc[i] == 'N':
                file_2.write(new_chromosome[i] + '\t' + str(int(new_position1[i]))  + '\t' + new_position[i] + '\t' + new_trinucleotide[i] + '\t' +  new_mutation[i] + '\t' + new_strand[i] + '\t' + multinucleotide[i] + '\t' + str(GetHomopolymerLength(multinucleotide[i], new_position1[i], new_position[i], len(multinucleotide[i])/2 + 1)[0]) + '\t' + str(lengths[i]))
                file_2.write('\n')
                if new_strand[i] == '+':
                    if multinucleotide[i] != 'NA':
                        file_3.write(multinucleotide[i])
                        file_3.write('\n')

                if new_strand[i] == '-':
                    if multinucleotide[i] != 'NA':
                            file_3.write(multinucleotide[i][::-1])
                            file_3.write('\n')


            elif new_mutation_type[i] == "deletion" and new_multiNuc[i] == 'Y':
                file_4.write(new_chromosome[i] + '\t' + str(int(new_position1[i]))  + '\t' + new_position[i] + '\t' + new_trinucleotide[i] + '\t' +  new_mutation[i] + '\t' + new_strand[i] + '\t' + multinucleotide[i] + '\t' + str(GetHomopolymerLength(multinucleotide[i], new_position1[i], new_position[i], len(multinucleotide[i])/2 + 1)[0]) + '\t' + str(lengths[i]))
                file_4.write('\n')
                if new_strand[i] == '+':
                    if multinucleotide[i] != 'NA':
                        file_5.write(multinucleotide[i])
                        file_5.write('\n')

                if new_strand[i] == '-':
                    if multinucleotide[i] != 'NA':
                            file_5.write(multinucleotide[i][::-1])
                            file_5.write('\n')
        #elif lengths[i] != 'NA' and (new_mutation[i] == 'T' or new_mutation[i] == 'A'):
        elif lengths[i] != 'NA' and (new_mutation[i] == 'C' or new_mutation[i] == 'G'):
    
            
            if new_mutation_type[i] == "insertion" and new_multiNuc[i] == 'N':
                file_6.write(new_chromosome[i] + '\t' + str(int(new_position1[i]))  + '\t' + new_position[i] + '\t' + new_trinucleotide[i] + '\t' +  new_mutation[i] + '\t' + new_strand[i] + '\t' + multinucleotide[i] + '\t' + str(GetHomopolymerLength(multinucleotide[i], new_position1[i], new_position[i], len(multinucleotide[i])/2 + 1)[0]) + '\t' + str(lengths[i]))
                file_6.write('\n')
                if new_strand[i] == '+':
                    if multinucleotide[i] != 'NA':
                        file_7.write(multinucleotide[i])
                        file_7.write('\n')
                if new_strand[i] == '-':
                    if multinucleotide[i] != 'NA':
                        file_7.write(multinucleotide[i][::-1])
                        file_7.write('\n')
            
            if new_mutation_type[i] == "insertion" and new_multiNuc[i] == 'Y':
                file_8.write(new_chromosome[i] + '\t' + str(int(new_position1[i]))  + '\t' + new_position[i] + '\t' + new_trinucleotide[i] + '\t' +  new_mutation[i] + '\t' + new_strand[i] + '\t' + multinucleotide[i] + '\t' + str(GetHomopolymerLength(multinucleotide[i], new_position1[i], new_position[i], len(multinucleotide[i])/2 + 1)[0]) + '\t' + str(lengths[i]))
                file_8.write('\n')
                if new_strand[i] == '+':
                    if multinucleotide[i] != 'NA':
                        file_9.write(multinucleotide[i])
                        file_9.write('\n')
                if new_strand[i] == '-':
                    if multinucleotide[i] != 'NA':
                        file_9.write(multinucleotide[i][::-1])
                        file_9.write('\n')

            
        

    return multinucleotide, new_chromosome, new_position, new_trinucleotide, lengths, new_allele, new_mutation, new_mutation_type, new_strand


f3 = open('acral_deletions_C.bed', 'w+')
f2 = open('acral_deletions.bed', 'w+')
f5 = open('acral_multiNucdeletions_C.bed', 'w+')
f4 = open('acral_multiNuc_deletions.bed', 'w+')

f13 = open('acral_insertions_C.bed', 'w+')
f12 = open('acral_insertions.bed', 'w+')
f15 = open('acral_multiNucinsertions_C.bed', 'w+')
f14 = open('acral_multiNuc_insertions.bed', 'w+')


f7 = open('cutaneous_deletions_C.bed', 'w+')
f6 = open('cutaneous_deletions.bed', 'w+')
f11 = open('cutaneous_multiNucdeletions_T.bed', 'w+')
f10 = open('cutaneous_multiNuc_deletions.bed', 'w+')
f17 = open('cutaneous_insertions_C.bed', 'w+')
f16 = open('cutaneous_insertions.bed', 'w+')
f19 = open('cutaneous_multiNucinsertions_C.bed', 'w+')
f18 = open('cutaneous_multiNuc_insertions.bed', 'w+')

acral_data = MakeMultiNucFile(acral, f2, f3, f4, f5, f12, f13, f14, f15)
chromosome = acral_data[1]
position = acral_data[2]
trinucleotide = acral_data[3]
lengths = acral_data[4]
allele = acral_data[5]
mutation = acral_data[6]
mutation_type = acral_data[7]
strand = acral_data[8]
f2.close()
f3.close()
f4.close()

cutaneous_data = MakeMultiNucFile(cutaneous, f6, f7, f10, f11, f16, f17, f18, f19)
chromosome2 = cutaneous_data[1]
position2 = cutaneous_data[2]
trinucleotide2 = cutaneous_data[3]
lengths2 = cutaneous_data[4]
allele2 = cutaneous_data[5]
mutation2 = cutaneous_data[6]
mutationtype_2 = cutaneous_data[7]
strand2 = cutaneous_data[8]






def MakeTriNucFile(chr_name, file1, file2, chromosome, position, trinucleotide, mutation, allele, strand, mutation_type, spectra, lengths ):
    #Makes bed files with single-base deletions
    file_data = []
    chromosome1_dict = {}
    chromosome1 =[]
    for x in range(0, len(position)):
        if trinucleotide[x].isalpha() and 'N' not in trinucleotide[x] and chromosome[x] == chr_name:
            if mutation_type[x] == 'deletion' and len(allele[x]) == 1:
                    if trinucleotide[x][1] == 'A' or trinucleotide[x][1] == 'G':
                        trinucleotide[x] = GetReverseComplement(trinucleotide[x])# len(trinucleotide[x]) == 3 and trinucleotide[x].isalpha():
                    if trinucleotide[x] not in spectra.keys():
                        spectra[trinucleotide[x]] = 1
                    else:
                        spectra[trinucleotide[x]]= spectra[trinucleotide[x]] + 1
                    if chromosome[x] == chr_name: 
                        file_data.append(( int(position[x]) - 1, int(position[x]), trinucleotide[x], mutation[x], strand[x], lengths[x]))
                        #chromosome1_dict[int(position1[x]), isolate[x]] = x
                


   
    new_file_data = sorted(file_data)
    for fd in range(0, len(new_file_data)):
        #if 'CC' in new_file_data[fd][2] or 'CT' in new_file_data[fd][2] or 'TC' in new_file_data[fd][2] or 'TT' in new_file_data[fd][2]:
        if int(new_file_data[fd][5]) >= 4: #Deletions in homopolymers
            file1.write(chr_name + '\t')  
            file1.write(str(new_file_data[fd][0]) + '\t' + str(new_file_data[fd][1]) + '\t' + new_file_data[fd][2] + '\t' + new_file_data[fd][3] + '\t' + new_file_data[fd][4])
            file1.write('\n')
            file2.write(chr_name + '\t')  
            file2.write(str(new_file_data[fd][0]) + '\t' + str(new_file_data[fd][1]) + '\t' + new_file_data[fd][2] + '\t' + new_file_data[fd][3] + '\t' + new_file_data[fd][4])
            file2.write('\n')
        else:
            file2.write(chr_name + '\t')  
            file2.write(str(new_file_data[fd][0]) + '\t' + str(new_file_data[fd][1]) + '\t' + new_file_data[fd][2] + '\t' + new_file_data[fd][3] + '\t' + new_file_data[fd][4])
            file2.write('\n')

    return spectra




def MakeInsertionFile(chr_name, file, file2, chromosome, position, trinucleotide, mutation, allele, strand, mutation_type, lengths ):
    #Makes bed files for single-base insertions
    file_data = []
    for x in range(0, len(position)):
        if mutation_type[x] == 'insertion' and len(mutation[x]) == 1:
                if trinucleotide[x].isalpha():# len(trinucleotide[x]) == 3 and trinucleotide[x].isalpha():
                    if chromosome[x] == chr_name: 
                        file_data.append(( int(position[x]) - 1, int(position[x]), trinucleotide[x], mutation[x], strand[x], lengths[x]))

   
    new_file_data = sorted(file_data)
    for fd in range(0, len(new_file_data)):
        #if new_file_data[fd][2][1] == 'A' or new_file_data[fd][2][1] == 'G':
        if int(new_file_data[fd][5] >= 4): #Insertions in homopolymers
            file.write(chr_name + '\t')  
            file.write(str(new_file_data[fd][0]) + '\t' + str(new_file_data[fd][1]) + '\t' + GetReverseComplement(new_file_data[fd][2][1:3]) + '\t' + new_file_data[fd][3] + '\t' + new_file_data[fd][4])
            file.write('\n')
            file2.write(chr_name + '\t')  
            file2.write(str(new_file_data[fd][0]) + '\t' + str(new_file_data[fd][1]) + '\t' + GetReverseComplement(new_file_data[fd][2][1:3]) + '\t' + new_file_data[fd][3] + '\t' + new_file_data[fd][4])
            file2.write('\n')
        else:
            file.write(chr_name + '\t')  
            file.write(str(new_file_data[fd][0]) + '\t' + str(new_file_data[fd][1]) + '\t' + new_file_data[fd][2][0:2] + '\t' + new_file_data[fd][3] + '\t' + new_file_data[fd][4])
            file.write('\n')



f5.close()
f6.close()
f7.close()
f10.close()
f11.close()
f12.close()
f13.close()
f14.close()
f15.close()
f16.close()
f17.close()
f18.close()
f19.close()

fa = open('Acral_singlebase_deletions_homopolymer.bed', 'w+')
fa2 = open('Acral_singlebase_deletions_all.bed', 'w+')
fb = open('Cutaneous_singlebase_deletions_homopolymer.bed', 'w+')
fb2 = open('Cutaneous_singlebase_deletions_all.bed', 'w+')
fc = open('Acral_singlebase_insertions.bed', 'w+')
fd = open('Cutaneous_singlebase_insertions.bed', 'w+')
fc2 = open('Acral_singlebase_insertions_homopolymer.bed', 'w+')
fd2 = open('Cutaneous_singlebase_insertions_homopolymer.bed', 'w+')

spectra = {}
spectra1 = MakeTriNucFile('chr1', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra, lengths)
spectra2 = MakeTriNucFile('chr2', fa, fa2, chromosome, position, trinucleotide, mutation, allele, strand, mutation_type, spectra1, lengths)
spectra3 = MakeTriNucFile('chr3', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra2, lengths)
spectra4 = MakeTriNucFile('chr4', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra3, lengths)
spectra5 = MakeTriNucFile('chr5', fa, fa2, chromosome, position, trinucleotide, mutation, allele, strand, mutation_type, spectra4, lengths)
spectra6 = MakeTriNucFile('chr6', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra5, lengths)
spectra7 = MakeTriNucFile('chr7', fa, fa2, chromosome, position, trinucleotide, mutation, allele, strand, mutation_type, spectra6, lengths)
spectra8 = MakeTriNucFile('chr8', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra7, lengths)
spectra9 = MakeTriNucFile('chr9', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra8, lengths)
spectra10 = MakeTriNucFile('chr10', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra9, lengths)
spectra11 = MakeTriNucFile('chr11', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra10, lengths)
spectra12 = MakeTriNucFile('chr12', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra11, lengths)
spectra13 = MakeTriNucFile('chr13', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra12, lengths)
spectra14 = MakeTriNucFile('chr14', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra13, lengths)
spectra15 = MakeTriNucFile('chr15', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra14, lengths)
spectra16 = MakeTriNucFile('chr16', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra15, lengths)
spectra17 = MakeTriNucFile('chr17', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra16, lengths)
spectra18 = MakeTriNucFile('chr18', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra17, lengths)
spectra19 = MakeTriNucFile('chr19', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra18, lengths)
spectra20 = MakeTriNucFile('chr20', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra19, lengths)
spectra21 = MakeTriNucFile('chr21', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra20, lengths)
spectra22 = MakeTriNucFile('chr22', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra21, lengths)
spectra23 = MakeTriNucFile('chrX', fa, fa2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, spectra21, lengths)
final_spectra = MakeTriNucFile('chrY', fa, fa2, chromosome, position, trinucleotide, mutation, allele, strand, mutation_type, spectra22, lengths)

spectrab = {}
spectra1b = MakeTriNucFile('chr1', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectrab, lengths2)
spectra2b = MakeTriNucFile('chr2', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra1b, lengths2)
spectra3b= MakeTriNucFile('chr3', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra2b, lengths2)
spectra4b = MakeTriNucFile('chr4', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra3b, lengths2)
spectra5b = MakeTriNucFile('chr5', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra4b, lengths2)
spectra6b = MakeTriNucFile('chr6', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra5b, lengths2)
spectra7b = MakeTriNucFile('chr7', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra6b, lengths2)
spectra8b = MakeTriNucFile('chr8', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra7b, lengths2)
spectra9b = MakeTriNucFile('chr9', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra8b, lengths2)
spectra10b = MakeTriNucFile('chr10', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra9b, lengths2)
spectra11b = MakeTriNucFile('chr11', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra10b, lengths2)
spectra12b = MakeTriNucFile('chr12', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra11b, lengths2)
spectra13b = MakeTriNucFile('chr13', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra12b, lengths2)
spectra14b = MakeTriNucFile('chr14', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra13b, lengths2)
spectra15b = MakeTriNucFile('chr15', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra14b, lengths2)
spectra16b = MakeTriNucFile('chr16', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra15b, lengths2)
spectra17b = MakeTriNucFile('chr17', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra16b, lengths2)
spectra18b = MakeTriNucFile('chr18', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra17b, lengths2)
spectra19b = MakeTriNucFile('chr19', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra18b, lengths2)
spectra20b = MakeTriNucFile('chr20', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra19b, lengths2)
spectra21b = MakeTriNucFile('chr21', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra20b, lengths2)
spectra22b = MakeTriNucFile('chr22', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra21b, lengths2)
spectra23b = MakeTriNucFile('chrX', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra21b, lengths2)
final_spectra2 = MakeTriNucFile('chrY', fb, fb2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, spectra22b, lengths2)

fa.close()
fa2.close()
fb.close()
fb2.close()

MakeInsertionFile('chr1', fc,fc2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, lengths)
MakeInsertionFile('chr2', fc, fc2, chromosome, position, trinucleotide, mutation, allele, strand, mutation_type, lengths)
MakeInsertionFile('chr3', fc, fc2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, lengths)
MakeInsertionFile('chr4', fc, fc2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, lengths)
MakeInsertionFile('chr5', fc, fc2, chromosome, position, trinucleotide, mutation, allele, strand, mutation_type, lengths)
MakeInsertionFile('chr6', fc, fc2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, lengths)
MakeInsertionFile('chr7', fc, fc2, chromosome, position, trinucleotide, mutation, allele, strand, mutation_type, lengths)
MakeInsertionFile('chr8', fc, fc2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, lengths)
MakeInsertionFile('chr9', fc, fc2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, lengths)
MakeInsertionFile('chr10', fc, fc2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, lengths)
MakeInsertionFile('chr11', fc, fc2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, lengths)
MakeInsertionFile('chr12', fc, fc2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, lengths)
MakeInsertionFile('chr13', fc, fc2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, lengths)
MakeInsertionFile('chr14', fc, fc2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, lengths)
MakeInsertionFile('chr15', fc, fc2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, lengths)
MakeInsertionFile('chr16', fc, fc2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, lengths)
MakeInsertionFile('chr17', fc, fc2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, lengths)
MakeInsertionFile('chr18', fc,fc2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, lengths)
MakeInsertionFile('chr19', fc, fc2,chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, lengths)
MakeInsertionFile('chr20', fc, fc2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, lengths)
MakeInsertionFile('chr21', fc, fc2, chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, lengths)
MakeInsertionFile('chr22', fc, fc2,  chromosome, position, trinucleotide, mutation,allele, strand, mutation_type, lengths)
MakeInsertionFile('chrX', fc, fc2, chromosome, position, trinucleotide, mutation, allele, strand, mutation_type, lengths)
MakeInsertionFile('chrY', fc, fc2, chromosome, position, trinucleotide, mutation, allele, strand, mutation_type, lengths)

MakeInsertionFile('chr1', fd, fd2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr2', fd, fd2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr3', fd, fd2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr4', fd, fd2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr5', fd, fd2,  chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr6', fd, fd2,  chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr7', fd, fd2,  chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr8', fd, fd2,  chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr9', fd, fd2,  chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr10', fd, fd2,  chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr11', fd, fd2,  chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr12', fd, fd2,  chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr13', fd, fd2,  chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr14', fd, fd2,  chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr15', fd, fd2,  chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr16', fd, fd2,  chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr17', fd, fd2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr18', fd, fd2, chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr19', fd, fd2,  chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr20', fd, fd2,  chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr21', fd, fd2,  chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chr22', fd, fd2,  chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chrX', fd, fd2,  chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)
MakeInsertionFile('chrY', fd, fd2,  chromosome2, position2, trinucleotide2, mutation2,allele2, strand2, mutationtype_2, lengths2)


fc.close()
fd.close()
fc2.close()
fd2.close()


#Calculates the number of mutations per tumor (140 for cutaneous, 35 for acral)
def MutationsPerTumor(mutation_type, allele, mutation, tumor_number):
   
    deletion_number = 0
    insertion_number = 0
  
    for i in range(0, len(mutation_type)):
        if mutation_type[i] == 'deletion':# and len(allele[i]) > 1:
            deletion_number = deletion_number + 1
        if mutation_type[i] == 'insertion:':# and len(mutation[i]) > 1:
            insertion_number = insertion_number + 1
            
    print(deletion_number / tumor_number)
    print(insertion_number / tumor_number)
  
    return deletion_number, insertion_number

def HomopolymerMutationsPerTumor(lengths, mutation_type, tumor_number, allele, mutation):
    nonhomopolymer = 0
    homopolymer = 0
    deletion_nonhomopolymer = 0
    deletion_homopolymer = 0
    insertion_nonhomopolymer = 0
    insertion_homopolymer = 0
    deletion_lengths = []
    insertion_lengths = []
    for i in range(0, len(lengths)):
        if lengths[i] != 'NA':
            if lengths[i] > 3:
                homopolymer = homopolymer + 1
                if mutation_type[i] == 'deletion' and len(allele[i]) == 1:
                    deletion_lengths.append(lengths[i])
                    deletion_homopolymer = deletion_homopolymer + 1
                if mutation_type[i] == 'insertion' and len(mutation[i]) == 1:
                    insertion_lengths.append(lengths[i])
                    insertion_homopolymer = insertion_homopolymer + 1
            else:
                nonhomopolymer = nonhomopolymer + 1
                if mutation_type[i] == 'deletion' and len(allele[i]) == 1:
                    deletion_lengths.append(lengths[i])
                    deletion_nonhomopolymer = deletion_nonhomopolymer + 1
                if mutation_type[i] == 'insertion' and len(mutation[i]) == 1:
                    insertion_lengths.append(lengths[i])
                    insertion_nonhomopolymer = insertion_nonhomopolymer + 1
    print(nonhomopolymer / tumor_number)
    print(homopolymer / tumor_number)
    print(deletion_nonhomopolymer/tumor_number)
    print(deletion_homopolymer/tumor_number)
    print(insertion_nonhomopolymer/tumor_number)
    print(insertion_homopolymer/tumor_number)
    print(deletion_lengths)
    print(insertion_lengths)
        
    return deletion_lengths, insertion_lengths


HomopolymerMutationsPerTumor(lengths, mutation_type, 35, allele, mutation)
HomopolymerMutationsPerTumor(lengths2, mutationtype_2, 140, allele2, mutation2)
MutationsPerTumor(mutation_type, allele, mutation, 35), 
MutationsPerTumor(mutationtype_2, allele2, mutation2, 140)


f8 = open('new_acral_polymer_lengths_single_deletions10+.txt', 'w+')
f9 = open('new_acral_polymer_lengths__single_insertions10+.txt', 'w+')
f11 = open('new_cutaneous_polymer_lengths_single_deletions10+.txt', 'w+')
f12 = open('new_cutaneous_polymer_lengths_single_insertions10+.txt', 'w+')
f13 = open('new_homopolymer_counts10+.txt', 'w+')





#Obtain the total number of bases in the human genome in sequences of each homopolymer length
def NormalizeHomopolymer(sequences, lengths, file):
    homopolymer_dict = {}
    C_dict = {}
    T_dict = {}
    homopolymer_counts = {}
    chromosome_sequences = {}
    chromosome_sequences['chrX'] = sequences[0]
    chromosome_sequences['chr1'] = sequences[1]
    chromosome_sequences['chr2'] = sequences[2]
    chromosome_sequences['chr3'] = sequences[3]
    chromosome_sequences['chr4'] = sequences[4]
    chromosome_sequences['chr5'] = sequences[5]
    chromosome_sequences['chr6'] = sequences[6]
    chromosome_sequences['chr7'] = sequences[7]
    chromosome_sequences['chr8'] = sequences[8]
    chromosome_sequences['chr9'] = sequences[9]
    chromosome_sequences['chr10'] = sequences[10]
    chromosome_sequences['chr11'] = sequences[11]
    chromosome_sequences['chr12'] = sequences[12]
    chromosome_sequences['chr13'] = sequences[13]
    chromosome_sequences['chr14'] = sequences[14]
    chromosome_sequences['chr15'] = sequences[15]
    chromosome_sequences['chr16'] = sequences[16]
    chromosome_sequences['chr17'] = sequences[17]
    chromosome_sequences['chr18'] = sequences[18]
    chromosome_sequences['chr19'] = sequences[19]
    chromosome_sequences['chr20'] = sequences[20]
    chromosome_sequences['chr21'] = sequences[21]
    chromosome_sequences['chr22'] = sequences[22]
    chromosome_sequences['chrY'] = sequences[23]



   
    homopolymer_dict['11-15'] = 0
    homopolymer_dict['16-20'] = 0
    homopolymer_dict['21+'] = 0
    homopolymer_dict['7+'] = 0
    homopolymer_dict['10+'] = 0
    C_dict['11-15'] = 0
    C_dict['16-20'] = 0
    C_dict['21+'] = 0
    C_dict['7+'] = 0
    C_dict['10+'] = 0
    T_dict['11-15'] = 0
    T_dict['16-20'] = 0
    T_dict['21+'] = 0
    T_dict['7+'] = 0
    T_dict['10+'] = 0

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
                
                if polymer_count in homopolymer_dict.keys():
                          
                            homopolymer_dict[polymer_count] = homopolymer_dict[polymer_count] + length
                               
                            if sequence[b] == 'C' or sequence[b] == 'G':
                                if polymer_count in C_dict.keys():
                                   
                                    C_dict[polymer_count] = C_dict[polymer_count] + length
                                else:
                                    C_dict[polymer_count] = length
                                    
                            if sequence[b] == 'T' or sequence[b] == 'A':
                                if polymer_count in T_dict.keys():
                                  
                                    T_dict[polymer_count] = T_dict[polymer_count] + length
                                else:
                                    T_dict[polymer_count] = length
                                 

    

                else:
                            homopolymer_dict[polymer_count] = length
                          
                            if sequence[b] == 'C' or sequence[b] == 'G':
                                C_dict[polymer_count] = length
                             
                            if sequence[b] == 'T' or sequence[b] == 'A':
                                T_dict[polymer_count] = length
                              
                if polymer_count >= 7:
                    homopolymer_dict['7+'] = homopolymer_dict['7+'] + length
                    if sequence[b] == 'C' or sequence[b] == 'G':
                            C_dict['7+'] = C_dict['7+'] + length
                    if sequence[b] == 'T' or sequence[b] == 'A':
                            T_dict['7+'] = T_dict['7+'] + length
                                 #for n in range (polymer_range[0], polymer_range[1] + 1):
                if polymer_count >= 10:
                    homopolymer_dict['10+'] = homopolymer_dict['10+'] + length
                    if sequence[b] == 'C' or sequence[b] == 'G':
                            C_dict['10+'] = C_dict['10+'] + length
                    if sequence[b] == 'T' or sequence[b] == 'A':
                            T_dict['10+'] = T_dict['10+'] + length
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




def PrintFrequencies(homopolymer_dict, C_dict, T_dict, lengths, mutation_type, allele, mutation, file, file2, tumor_number ):
    homopolymer_counts = {}
    deletion_counts = {}
    insertion_counts = {}
    C_homopolymer_counts = {}
    C_deletion_counts = {}
    C_insertion_counts = {}
    T_homopolymer_counts = {}
    T_deletion_counts = {}
    T_insertion_counts = {}
    homopolymer_counts['7+'] = 0
    deletion_counts['7+'] = 0
    insertion_counts['7+'] = 0
    homopolymer_counts['10+'] = 0
    deletion_counts['10+'] = 0
    insertion_counts['10+'] = 0
    homopolymer_counts['21+'] = 0
    deletion_counts['21+'] = 0
    insertion_counts['21+'] = 0
    homopolymer_counts['16-20'] = 0
    deletion_counts['16-20'] = 0
    insertion_counts['16-20'] = 0
    homopolymer_counts['11-15'] = 0
    deletion_counts['11-15'] = 0
    insertion_counts['11-15'] = 0
    C_homopolymer_counts['7+'] = 0
    C_deletion_counts['7+'] = 0
    C_insertion_counts['7+'] = 0
    C_homopolymer_counts['10+'] = 0
    C_deletion_counts['10+'] = 0
    C_insertion_counts['10+'] = 0
    C_homopolymer_counts['21+'] = 0
    C_deletion_counts['21+'] = 0
    C_insertion_counts['21+'] = 0
    C_homopolymer_counts['16-20'] = 0
    C_deletion_counts['16-20'] = 0
    C_insertion_counts['16-20'] = 0
    C_homopolymer_counts['11-15'] = 0
    C_deletion_counts['11-15'] = 0
    C_insertion_counts['11-15'] = 0
    T_homopolymer_counts['7+'] = 0
    T_deletion_counts['7+'] = 0
    T_insertion_counts['7+'] = 0
    T_homopolymer_counts['10+'] = 0
    T_deletion_counts['10+'] = 0
    T_insertion_counts['10+'] = 0
    T_homopolymer_counts['21+'] = 0
    T_deletion_counts['21+'] = 0
    T_insertion_counts['21+'] = 0
    T_homopolymer_counts['16-20'] = 0
    T_deletion_counts['16-20'] = 0
    T_insertion_counts['16-20'] = 0
    T_homopolymer_counts['11-15'] = 0
    T_deletion_counts['11-15'] = 0
    T_insertion_counts['11-15'] = 0
    for l in range(0, len(lengths)):
        if allele[l] != 'NA' and lengths[l] != 'NA':
            if (len(allele[l]) == 1 and mutation_type[l] == "deletion") or (len(mutation[l]) == 1 and mutation_type[l] == "insertion") :
                if lengths[l] not in homopolymer_counts.keys():
                    homopolymer_counts[lengths[l]] = 1
                    #homopolymer_dict[lengths[l]] = 0
                else:
                    homopolymer_counts[lengths[l]] = homopolymer_counts[lengths[l]] + 1
                if mutation_type[l] == "deletion":
                    if lengths[l] not in deletion_counts.keys():
                        deletion_counts[lengths[l]] = 1
                    else:
                        deletion_counts[lengths[l]] = deletion_counts[lengths[l]] + 1
                    if int(lengths[l]) >= 7:
                        homopolymer_counts['7+'] = homopolymer_counts['7+'] + 1
                        deletion_counts['7+'] = deletion_counts['7+'] + 1
                    if int(lengths[l]) >= 10:
                        homopolymer_counts['10+'] = homopolymer_counts['10+'] + 1
                        deletion_counts['10+'] = deletion_counts['10+'] + 1
                    if int(lengths[l]) >= 11:
                        if int(lengths[l]) >= 21:
                            homopolymer_counts['21+'] = homopolymer_counts['21+'] + 1
                            deletion_counts['21+'] = deletion_counts['21+'] + 1
                        elif int(lengths[l]) >= 16:
                            homopolymer_counts['16-20'] = homopolymer_counts['16-20'] + 1
                            deletion_counts['16-20'] = deletion_counts['16-20'] + 1
                        else:
                            homopolymer_counts['11-15'] = homopolymer_counts['11-15'] + 1
                            deletion_counts['11-15'] = deletion_counts['11-15'] + 1
                if mutation_type[l] == "insertion":
                    if lengths[l] not in insertion_counts.keys():
                        insertion_counts[lengths[l]] = 1
                    else:
                        insertion_counts[lengths[l]] = insertion_counts[lengths[l]] + 1
                    if int(lengths[l]) >= 7:
                        homopolymer_counts['7+'] = homopolymer_counts['7+'] + 1
                        insertion_counts['7+'] = insertion_counts['7+'] + 1
                    if int(lengths[l]) >= 10:
                        homopolymer_counts['10+'] = homopolymer_counts['10+'] + 1
                        insertion_counts['10+'] = insertion_counts['10+'] + 1
                    if int(lengths[l]) >= 11:
                        if int(lengths[l]) >= 21:
                            homopolymer_counts['21+'] = homopolymer_counts['21+'] + 1
                            insertion_counts['21+'] = insertion_counts['21+'] + 1
                        elif int(lengths[l]) >= 16:
                            homopolymer_counts['16-20'] = homopolymer_counts['16-20'] + 1
                            insertion_counts['16-20'] = insertion_counts['16-20'] + 1
                        else:
                            homopolymer_counts['11-15'] = homopolymer_counts['11-15'] + 1
                            insertion_counts['11-15'] = insertion_counts['11-15'] + 1
                
                if allele[l] == 'C' or allele[l] == 'G'or mutation[l] == 'C' or mutation[l] == 'G':
                    if lengths[l] not in C_homopolymer_counts.keys():
                        C_homopolymer_counts[lengths[l]] = 1
                        #homopolymer_dict[lengths[l]] = 0
                    else:
                        C_homopolymer_counts[lengths[l]] = C_homopolymer_counts[lengths[l]] + 1
                    if mutation_type[l] == "deletion":
                        if lengths[l] not in C_deletion_counts.keys():
                            C_deletion_counts[lengths[l]] = 1
                        else:
                            C_deletion_counts[lengths[l]] = C_deletion_counts[lengths[l]] + 1
                        if int(lengths[l]) >= 7:
                            C_homopolymer_counts['7+'] = C_homopolymer_counts['7+'] + 1
                            C_deletion_counts['7+'] = C_deletion_counts['7+'] + 1
                        if int(lengths[l]) >= 10:
                            C_homopolymer_counts['10+'] = C_homopolymer_counts['10+'] + 1
                            C_deletion_counts['10+'] = C_deletion_counts['10+'] + 1
                        if int(lengths[l]) >= 11:
                            if int(lengths[l]) >= 21:
                                C_homopolymer_counts['21+'] = C_homopolymer_counts['21+'] + 1
                                C_deletion_counts['21+'] = C_deletion_counts['21+'] + 1
                            elif int(lengths[l]) >= 16:
                                C_homopolymer_counts['16-20'] = C_homopolymer_counts['16-20'] + 1
                                C_deletion_counts['16-20'] = C_deletion_counts['16-20'] + 1
                            else:
                                C_homopolymer_counts['11-15'] = C_homopolymer_counts['11-15'] + 1
                                C_deletion_counts['11-15'] = C_deletion_counts['11-15'] + 1
                    if mutation_type[l] == "insertion":
                        if lengths[l] not in C_insertion_counts.keys():
                            C_insertion_counts[lengths[l]] = 1
                        else:
                            C_insertion_counts[lengths[l]] = C_insertion_counts[lengths[l]] + 1
                        if int(lengths[l]) >= 7:
                            C_homopolymer_counts['7+'] = C_homopolymer_counts['7+'] + 1
                            C_insertion_counts['7+'] = C_insertion_counts['7+'] + 1
                        if int(lengths[l]) >= 10:
                            C_homopolymer_counts['10+'] = C_homopolymer_counts['10+'] + 1
                            C_insertion_counts['10+'] = C_insertion_counts['10+'] + 1
                        if int(lengths[l]) >= 11:
                            if int(lengths[l]) >= 21:
                                C_homopolymer_counts['21+'] = C_homopolymer_counts['21+'] + 1
                                C_insertion_counts['21+'] = C_insertion_counts['21+'] + 1
                            elif int(lengths[l]) >= 16:
                                C_homopolymer_counts['16-20'] = C_homopolymer_counts['16-20'] + 1
                                C_insertion_counts['16-20'] = C_insertion_counts['16-20'] + 1
                            else:
                                C_homopolymer_counts['11-15'] = C_homopolymer_counts['11-15'] + 1
                                C_insertion_counts['11-15'] = C_insertion_counts['11-15'] + 1
                if allele[l] == 'T' or allele[l] == 'A'or mutation[l] == 'T' or mutation[l] == 'A':
                    if lengths[l] not in T_homopolymer_counts.keys():
                        T_homopolymer_counts[lengths[l]] = 1
                        #homopolymer_dict[lengths[l]] = 0
                    else:
                        T_homopolymer_counts[lengths[l]] = T_homopolymer_counts[lengths[l]] + 1
                    if mutation_type[l] == "deletion":
                        if lengths[l] not in T_deletion_counts.keys():
                            T_deletion_counts[lengths[l]] = 1
                        else:
                            T_deletion_counts[lengths[l]] = T_deletion_counts[lengths[l]] + 1
                        if int(lengths[l]) >= 7:
                            T_homopolymer_counts['7+'] = T_homopolymer_counts['7+'] + 1
                            T_deletion_counts['7+'] = T_deletion_counts['7+'] + 1
                        if int(lengths[l]) >= 10:
                            T_homopolymer_counts['10+'] = T_homopolymer_counts['10+'] + 1
                            T_deletion_counts['10+'] = T_deletion_counts['10+'] + 1
                        if int(lengths[l]) >= 11:
                            if int(lengths[l]) >= 21:
                                T_homopolymer_counts['21+'] = T_homopolymer_counts['21+'] + 1
                                T_deletion_counts['21+'] = T_deletion_counts['21+'] + 1
                            elif int(lengths[l]) >= 16:
                                T_homopolymer_counts['16-20'] = T_homopolymer_counts['16-20'] + 1
                                T_deletion_counts['16-20'] = T_deletion_counts['16-20'] + 1
                            else:
                                T_homopolymer_counts['11-15'] = T_homopolymer_counts['11-15'] + 1
                                T_deletion_counts['11-15'] = T_deletion_counts['11-15'] + 1
                    if mutation_type[l] == "insertion":
                        if lengths[l] not in T_insertion_counts.keys():
                            T_insertion_counts[lengths[l]] = 1
                        else:
                            T_insertion_counts[lengths[l]] = T_insertion_counts[lengths[l]] + 1
                        if int(lengths[l]) >= 7:
                            T_homopolymer_counts['7+'] = T_homopolymer_counts['7+'] + 1
                            T_insertion_counts['7+'] = T_insertion_counts['7+'] + 1
                        if int(lengths[l]) >= 10:
                            T_homopolymer_counts['10+'] = T_homopolymer_counts['10+'] + 1
                            T_insertion_counts['10+'] = T_insertion_counts['10+'] + 1
                        if int(lengths[l]) >= 11:
                            if int(lengths[l]) >= 21:
                                T_homopolymer_counts['21+'] = T_homopolymer_counts['21+'] + 1
                                T_insertion_counts['21+'] = T_insertion_counts['21+'] + 1
                            elif int(lengths[l]) >= 16:
                                T_homopolymer_counts['16-20'] = T_homopolymer_counts['16-20'] + 1
                                T_insertion_counts['16-20'] = T_insertion_counts['16-20'] + 1
                            else:
                                T_homopolymer_counts['11-15'] = T_homopolymer_counts['11-15'] + 1
                                T_insertion_counts['11-15'] = T_insertion_counts['11-15'] + 1
    for key in deletion_counts.keys():
        if homopolymer_dict[key] != 0:
            deletion_counts[key] = deletion_counts[key]/homopolymer_dict[key]
            file.write(str(key) + ':' + str(deletion_counts[key]/tumor_number))
            file.write('\n')
    for key in insertion_counts.keys():
        if homopolymer_dict[key] != 0:
            insertion_counts[key] = insertion_counts[key]/homopolymer_dict[key]
            file2.write(str(key) + ':' + str(insertion_counts[key]/tumor_number))
            file2.write('\n')
    for key in C_deletion_counts.keys():
        if C_dict[key] != 0:
            C_deletion_counts[key] = C_deletion_counts[key]/C_dict[key]
            file.write(str(key) + ':' + str(C_deletion_counts[key]/tumor_number))
            file.write('\n')
    for key in C_insertion_counts.keys():
        if C_dict[key] != 0:
            C_insertion_counts[key] = C_insertion_counts[key]/C_dict[key]
            file2.write(str(key) + ':' + str(C_insertion_counts[key]/ tumor_number))
            file2.write('\n')
    for key in T_deletion_counts.keys():
        if T_dict[key] != 0:
            T_deletion_counts[key] = T_deletion_counts[key]/T_dict[key]
            file.write(str(key) + ':' + str(T_deletion_counts[key]/ tumor_number))
            file.write('\n')
    for key in T_insertion_counts.keys():
        if T_dict[key] != 0:
            T_insertion_counts[key] = T_insertion_counts[key]/T_dict[key]
            file2.write(str(key) + ':' + str(T_insertion_counts[key] / tumor_number))
            file2.write('\n')


homopolymer_numbers = NormalizeHomopolymer(sequences, lengths, f13)       
PrintFrequencies(homopolymer_numbers[0], homopolymer_numbers[1], homopolymer_numbers[2], lengths, mutation_type, allele, mutation, f8, f9, 35)
PrintFrequencies(homopolymer_numbers[0], homopolymer_numbers[1], homopolymer_numbers[2], lengths2, mutationtype_2, allele2, mutation2, f11, f12, 140)
coordinates = homopolymer_numbers[3]


f8.close()
f9.close()
f11.close()
f12.close()
f13.close()

