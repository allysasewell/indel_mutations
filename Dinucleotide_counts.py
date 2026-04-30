
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


di_counts = {0:0, 1:0, 2:0, 3:0, 4:0, 5:0}

#Counts the number of bases in dipyrimidine and non-dipyrimidine sequence contexts
def DinucleotideContext(file, di_dict):
  
    ft = open(file)
    sequence = ft.read()
    sequence = str(sequence.strip().upper())
    sequence = ''.join(sequence.split('\n'))
    
    
    for i in range(1, len(sequence) - 1):
        if 'CC' in sequence[i-1:i+2] or 'CT' in sequence[i-1:i+2] or 'TC' in sequence[i-1:i+2] or 'TT' in sequence[i-1:i+2] or 'AA' in sequence[i-1:i+2] or  'AG' in sequence[i-1:i+2] or 'GA' in sequence[i-1:i+2] or 'GG' in sequence[i-1:i+2]:
            di_dict[0] = di_dict[0] + 1
            if sequence[i-1:i+2][1] == 'C' or sequence[i-1:i+2][1] == 'G':
                di_dict[2] = di_dict[2] + 1
            else:
                di_dict[4] = di_dict[4] + 1

        else:
            di_dict[1] = di_dict[1] + 1
            if sequence[i-1:i+2][1] == 'C' or sequence[i-1:i+2][1] == 'G':
                di_dict[3] = di_dict[3] + 1
            else:
                di_dict[5] = di_dict[5] + 1
    
       
        


                                         
    ft.close()
    return di_dict
di_dict1 = DinucleotideContext('chr1.txt',  di_counts)
di_dict2 = DinucleotideContext('chr2.txt',  di_dict1)
di_dict3 = DinucleotideContext('chr3.txt', di_dict2)
di_dict4 = DinucleotideContext('chr4.txt',  di_dict3)
di_dict5 = DinucleotideContext('chr5.txt', di_dict4)
di_dict6 = DinucleotideContext('chr6.txt',  di_dict5)
di_dict7 = DinucleotideContext('chr7.txt', di_dict6)
di_dict8 = DinucleotideContext('chr8.txt',  di_dict7)
di_dict9 = DinucleotideContext('chr9.txt',  di_dict8)
di_dict10 = DinucleotideContext('chr10.txt',  di_dict9)
di_dict11 = DinucleotideContext('chr11.txt', di_dict10)
di_dict12 = DinucleotideContext('chr12.txt',  di_dict11)
di_dict13 = DinucleotideContext('chr13.txt',  di_dict12)
di_dict14 = DinucleotideContext('chr14.txt',  di_dict13)
di_dict15 = DinucleotideContext('chr15.txt',  di_dict14)
dinuc_counts = DinucleotideContext('chr16.txt',  di_dict15)
print(dinuc_counts)

def CountTotals(trinuc_counts):
    dipyrimidine_count = 0
    nondipyrimidine_count = 0
    for key in trinuc_counts.keys():
        if 'TT' in key or 'TC' in key or 'CT' in key or 'CC' in key:
            dipyrimidine_count = dipyrimidine_count + int(trinuc_counts[key])
        else:
            nondipyrimidine_count = nondipyrimidine_count + int(trinuc_counts[key])
    return dipyrimidine_count, nondipyrimidine_count


#CountTotals(trinuc_counts)


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
print(len(sequenceX))
sequence1 = GetSequence('hg19.fa', 'chr1')
sequences.append(sequence1)
print(len(sequence1))
sequence2 = GetSequence('hg19.fa', 'chr2')
sequences.append(sequence2)
print(len(sequence2))
sequence3 = GetSequence('hg19.fa', 'chr3')
print(len(sequence3))
sequences.append(sequence3)
sequence4 = GetSequence('hg19.fa', 'chr4')
sequences.append(sequence4)
print(len(sequence4))
sequence5 = GetSequence('hg19.fa', 'chr5')
sequences.append(sequence5)
print(len(sequence5))
sequence6 = GetSequence('hg19.fa', 'chr6')
sequences.append(sequence6)
print(len(sequence6))
sequence7 = GetSequence('hg19.fa', 'chr7')
sequences.append(sequence7)
print(len(sequence7))
sequence8 = GetSequence('hg19.fa', 'chr8')
sequences.append(sequence8)
print(len(sequence8))
sequence9 = GetSequence('hg19.fa', 'chr9')
sequences.append(sequence9)
print(len(sequence9))
sequence10 = GetSequence('hg19.fa', 'chr10')
sequences.append(sequence10)
print(len(sequence10))
sequence11 = GetSequence('hg19.fa', 'chr11')
sequences.append(sequence11)
print(len(sequence11))
sequence12 = GetSequence('hg19.fa', 'chr12')
sequences.append(sequence12)
print(len(sequence12))
sequence13 = GetSequence('hg19.fa', 'chr13')
sequences.append(sequence13)
print(len(sequence13))
sequence14 = GetSequence('hg19.fa', 'chr14')
sequences.append(sequence14)
print(len(sequence14))
sequence15 = GetSequence('hg19.fa', 'chr15')
sequences.append(sequence15)
print(len(sequence15))
sequence16 = GetSequence('hg19.fa', 'chr16')
sequences.append(sequence16)
print(len(sequence16))
sequence17 = GetSequence('hg19.fa', 'chr17')
sequences.append(sequence17)
print(len(sequence17))
sequence18 = GetSequence('hg19.fa', 'chr18')
sequences.append(sequence18)
print(len(sequence18))
sequence19 = GetSequence('hg19.fa', 'chr19')
sequences.append(sequence19)
print(len(sequence19))
sequence20 = GetSequence('hg19.fa', 'chr20')
sequences.append(sequence20)
print(len(sequence20))
sequence21 = GetSequence('hg19.fa', 'chr21')
print(len(sequence21))
sequences.append(sequence21)
sequence22 = GetSequence('hg19.fa', 'chr22')
sequences.append(sequence22)
print(len(sequence22))
sequenceY = GetSequence('hg19.fa', 'chrY')
sequences.append(sequenceY)
print(len(sequenceY))

#Counts bases in dipyrimidine and non-dipyrimidine sequence contexts in the human genome
def DinucleotideContext2(sequences):
    dipyrimidine = 0
    C_dipyrimidine = 0
    T_dipyrimidine = 0
    C_nondipyrimidine = 0
    T_nondipyrimidine = 0
    non_dipyrimidine = 0
    for sequence in sequences:
        
       
    
        for i in range(1, len(sequence) - 1):
            if 'CC' in sequence[i-1:i+2] or 'CT' in sequence[i-1:i+2] or 'TC' in sequence[i-1:i+2] or 'TT' in sequence[i-1:i+2] or 'AA' in sequence[i-1:i+2] or  'AG' in sequence[i-1:i+2] or 'GA' in sequence[i-1:i+2] or 'GG' in sequence[i-1:i+2]:
                dipyrimidine = dipyrimidine + 1
                if sequence[i-1:i+2][1] == 'G' or sequence[i-1:i+2][1] == 'C':
                    C_dipyrimidine = C_dipyrimidine + 1
                else:
                    T_dipyrimidine = T_dipyrimidine + 1

            else:
                non_dipyrimidine = non_dipyrimidine + 1
                if sequence[i-1:i+2][1] == 'G' or sequence[i-1:i+2][1] == 'C':
                    C_nondipyrimidine = C_nondipyrimidine + 1
                else:
                    T_nondipyrimidine = T_nondipyrimidine + 1

        


                                         

    return dipyrimidine, non_dipyrimidine




dinuc_hg19 = DinucleotideContext2(sequences)
