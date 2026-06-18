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





#Takes file from perl scripts and finds total transcribed strand and non-transcribed strand counts
def ContextCounts(file1, file2):
    context = []
    TS = []
    TS_C = []
    TS_T = []
    NTS = []
    NTS_C_count = 0
    NTS_T_count = 0
    TS_C_count = 0
    TS_T_count = 0
    intergenic = []
    data = file1.read()
    lines = [s.strip().split() for s in data.splitlines()]
    for l in range (2, len(lines)):
        context.append(lines[l][0])
        TS.append(int(lines[l][3]))

        NTS.append(int(lines[l][4]))
        if lines[l][2][0] == 'C' or lines[l][2][0] == 'G':
            TS_C_count = TS_C_count + int(lines[l][3])
            NTS_C_count = NTS_C_count + int(lines[l][4])
        if lines[l][2][0] == 'T' or lines[l][2][0] == 'A':
            TS_T_count = TS_T_count + int(lines[l][3])
            NTS_T_count = NTS_T_count + int(lines[l][4])
        intergenic.append(int(lines[l][5]))
    dipyrimidine_count = 0
    other_count = 0
    genic_count = 0
    intergenic_count = 0
    transcribed_count = 0
    non_transcribed_count = 0
    for a in range(0, len(TS)):
        transcribed_count = transcribed_count + TS[a]
        non_transcribed_count = non_transcribed_count + NTS[a]
        genic_count = genic_count + TS[a] + NTS[a]
        intergenic_count = intergenic_count + intergenic[a]
        if 'TT' in context[a] or 'CC' in context[a] or 'CT' in context[a] or 'TC' in context[a]:
            dipyrimidine_count = dipyrimidine_count + TS[a] + NTS[a] + intergenic[a]
        else:
            other_count = other_count + TS[a] + NTS[a] + intergenic[a]
    file2.write('Transcribed Count: ' + str(transcribed_count))
    file2.write('\n')
    file2.write('Non-Transcribed Count: ' + str(non_transcribed_count))
    file2.write('\n')
    file2.write('Genic Count: ' + str(genic_count))
    file2.write('\n')
    file2.write('Intergenic Count: ' + str(intergenic_count))
    file2.write('\n')
    file2.write('Ratio: ' + str(genic_count/intergenic_count))
    file2.write('\n')
    file2.write('Dipyrimidine Count: ' + str(dipyrimidine_count))
    file2.write('\n')
    file2.write('Non-dipyrimidine Count: ' + str(other_count))
    file2.write('\n')
    if other_count != 0:
        file2.write('Ratio: ' + str(dipyrimidine_count/other_count))
    file2.close()
    return non_transcribed_count, transcribed_count




d = open("WT_singlebase_deletions_Trinuc_TrxStrandBins.txt")
d1 = open("rad16_singlebase_deletions_Trinuc_TrxStrandBins.txt")
d2 = open("rad26_singlebase_deletions_Trinuc_TrxStrandBins.txt")
d3 = open("rad30_singlebase_deletions_Trinuc_TrxStrandBins.txt")
e = open("WT_deletions_strandresults.txt", 'w+')
e1 = open("rad16_deletions_strandresults.txt", 'w+')
e2 = open("rad26_deletions_strandresults.txt", 'w+')
e3 = open("rad30_deletions_strandresults.txt", 'w+')


d4 = open("WT_singlebase_deletions_Trinuc_homopolymer_TrxStrandBins.txt")
d5 = open("rad16_singlebase_deletions_Trinuc_homopolymer_TrxStrandBins.txt")
d6 = open("rad26_singlebase_deletions_Trinuc_homopolymer_TrxStrandBins.txt")
d10 = open("rad30_singlebase_deletions_Trinuc_homopolymer_TrxStrandBins.txt")
d7 = open("WT_singlebase_insertions_Trinuc_homopolymer_TrxStrandBins.txt")
d8 = open("rad16_singlebase_insertions_Trinuc_homopolymer_TrxStrandBins.txt")
d9 = open("rad26_singlebase_insertions_Trinuc_homopolymer_TrxStrandBins.txt")
d11 = open("rad30_singlebase_insertions_Trinuc_homopolymer_TrxStrandBins.txt")
e4 = open("WT_homopolymer_singledeletions_strandresults.txt", 'w+')
e5 = open("rad16_homopolymer_singledeletions_strandresults.txt", 'w+')
e6 = open("rad26_homopolymer_singledeletions_strandresults.txt", 'w+')
e10 = open("rad30_homopolymer_singledeletions_strandresults.txt", 'w+')
e7 = open("WT_homopolymer_singleinsertions_strandresults.txt", 'w+')
e8 = open("rad16_homopolymer_singleinsertions_strandresults.txt", 'w+')
e9 = open("rad26_homopolymer_singleinsertions_strandresults.txt", 'w+')
e11 = open("rad30_homopolymer_singleinsertions_strandresults.txt", 'w+')


h = open("UVB_WT_singlebase_deletions_Trinuc_TrxStrandBins.txt")
h1 = open("UVB_rad16_singlebase_deletions_Trinuc_TrxStrandBins.txt")
h2 = open("UVB_rad26_singlebase_deletions_Trinuc_TrxStrandBins.txt")
h3 = open("UVB_rad30_singlebase_deletions_Trinuc_TrxStrandBins.txt")
i = open("UVB_WT_deletions_strandresults.txt", 'w+')
i1 = open("UVB_rad16_deletions_strandresults.txt", 'w+')
i2 = open("UVB_rad26_deletions_strandresults.txt", 'w+')
i3 = open("UVB_rad30_deletions_strandresults.txt", 'w+')



h4 = open("UVB_WT_singlebase_deletions_Trinuc_homopolymer_TrxStrandBins.txt")
h5 = open("UVB_rad16_singlebase_deletions_Trinuc_homopolymer_TrxStrandBins.txt")
h6 = open("UVB_rad26_singlebase_deletions_Trinuc_homopolymer_TrxStrandBins.txt")
h7 = open("UVB_WT_singlebase_insertions_Trinuc_homopolymer_TrxStrandBins.txt")
h8 = open("UVB_rad16_singlebase_insertions_Trinuc_homopolymer_TrxStrandBins.txt")
h9 = open("UVB_rad26_singlebase_insertions_Trinuc_homopolymer_TrxStrandBins.txt")
i4 = open("UVB_WT_homopolymer_singledeletions_strandresults.txt", 'w+')
i5 = open("UVB_rad16_homopolymer_singledeletions_strandresults.txt", 'w+')
i6 = open("UVB_rad26_homopolymer_singledeletions_strandresults.txt", 'w+')
i7 = open("UVB_WT_homopolymer_singleinsertions_strandresults.txt", 'w+')
i8 = open("UVB_rad16_homopolymer_singleinsertions_strandresults.txt", 'w+')
i9 = open("UVB_rad26_homopolymer_singleinsertions_strandresults.txt", 'w+')

j = open("WT_singlebase_deletions_TriNuc_nondipyrimidine_TrxStrandBins.txt")
j1 = open("WT_singlebase_deletions_TriNuc_dipyrimidine_TrxStrandBins.txt")
k = open("rad16_singlebase_deletions_TriNuc_nondipyrimidine_TrxStrandBins.txt")
k1 = open("rad16_singlebase_deletions_TriNuc_dipyrimidine_TrxStrandBins.txt")
m = open("rad26_singlebase_deletions_TriNuc_nondipyrimidine_TrxStrandBins.txt")
m1 = open("rad26_singlebase_deletions_TriNuc_dipyrimidine_TrxStrandBins.txt")
n = open("rad30_singlebase_deletions_TriNuc_nondipyrimidine_TrxStrandBins.txt")
n1 = open("rad30_singlebase_deletions_TriNuc_dipyrimidine_TrxStrandBins.txt")
j2 = open("WT_deletions_nondipyrimidine_strandresults.txt", 'w+')
j3 = open("WT_deletions_dipyrimidine_strandresults.txt", 'w+')
k2 = open("rad16_deletions_nondipyrimidine_strandresults.txt", 'w+')
k3 = open("rad16_deletions_dipyrimidine_strandresults.txt", 'w+')
m2 = open("rad26_deletions_nondipyrimidine_strandresults.txt", 'w+')
m3 = open("rad26_deletions_dipyrimidine_strandresults.txt", 'w+')
n2 = open("rad30_deletions_nondipyrimidine_strandresults.txt", 'w+')
n3 = open("rad30_deletions_dipyrimidine_strandresults.txt", 'w+')

p = open("UVB_WT_singlebase_deletions_Trinuc_nondipyrimidine_TrxStrandBins.txt")
p1 = open("UVB_WT_singlebase_deletions_Trinuc_dipyrimidine_TrxStrandBins.txt")
q = open("UVB_rad16_singlebase_deletions_Trinuc_nondipyrimidine_TrxStrandBins.txt")
q1 = open("UVB_rad16_singlebase_deletions_Trinuc_dipyrimidine_TrxStrandBins.txt")
r = open("UVB_rad26_singlebase_deletions_Trinuc_nondipyrimidine_TrxStrandBins.txt")
r1 = open("UVB_rad26_singlebase_deletions_Trinuc_dipyrimidine_TrxStrandBins.txt")
s = open("UVB_rad30_singlebase_deletions_Trinuc_nondipyrimidine_TrxStrandBins.txt")
s1 = open("UVB_rad30_singlebase_deletions_Trinuc_dipyrimidine_TrxStrandBins.txt")
p2 = open("UVB_WT_deletions_nondipyrimidine_strandresults.txt", 'w+')
p3 = open("UVB_WT_deletions_dipyrimidine_strandresults.txt", 'w+')
q2 = open("UVB_rad16_deletions_nondipyrimidine_strandresults.txt", 'w+')
q3 = open("UVB_rad16_deletions_dipyrimidine_strandresults.txt", 'w+')
r2 = open("UVB_rad26_deletions_nondipyrimidine_strandresults.txt", 'w+')
r3 = open("UVB_rad26_deletions_dipyrimidine_strandresults.txt", 'w+')
s2 = open("UVB_rad30_deletions_nondipyrimidine_strandresults.txt", 'w+')
s3 = open("UVB_rad30_deletions_dipyrimidine_strandresults.txt", 'w+')


t = open("WT_singlebase_insertions_TrxStrandBins.txt")
t1 = open("rad16_singlebase_insertions_TrxStrandBins.txt")
t2 = open("rad26_singlebase_insertions_TrxStrandBins.txt")
t3 = open("rad30_singlebase_insertions_TrxStrandBins.txt")
u = open("WT_insertions_strandresults.txt", 'w+')
u1 = open("rad16_insertions_strandresults.txt", 'w+')
u2 = open("rad26_insertions_strandresults.txt", 'w+')
u3 = open("rad30_insertions_strandresults.txt", 'w+')
v = open("UVB_WT_singlebase_insertions_TrxStrandBins.txt")
v1 = open("UVB_rad16_singlebase_insertions_TrxStrandBins.txt")
v2 = open("UVB_rad26_singlebase_insertions_TrxStrandBins.txt")
v3 = open("UVB_rad30_singlebase_insertions_TrxStrandBins.txt")
w = open("UVB_WT_insertions_strandresults.txt", 'w+')
w1 = open("UVB_rad16_insertions_strandresults.txt", 'w+')
w2 = open("UVB_rad26_insertions_strandresults.txt", 'w+')
w3 = open("UVB_rad30_insertions_strandresults.txt", 'w+')


z1 = open("WT_singlebase_deletionsv2_Trinuc_TrxStrandBins.txt")
z2 = open("rad30_singlebase_deletions_Trinuc_TrxStrandBins.txt")
z3 = open("WT_singlebase_deletionsv2_strandresults.txt", 'w+')
z4 = open("rad30_singlebase_deletions_strandresults.txt", 'w+')

z5 = open("WT_singlebase_deletionsv2_Trinuc_homopolymer_TrxStrandBins.txt")
z6 = open("rad30_singlebase_deletions_Trinuc_homopolymer_TrxStrandBins.txt")
z7 = open("WT_singlebase_deletions_homopolymer_strandresults.txt", 'w+')
z8 = open("rad30_singlebase_deletions_homopolymer_strandresults.txt", 'w+')

z9 = open("WT_singlebase_insertionsv2_Trinuc_homopolymer_TrxStrandBins.txt")
z10 = open("rad30_singlebase_insertions_Trinuc_homopolymer_TrxStrandBins.txt")
z11 = open("WT_singlebase_insertionsv2_homopolymer_strandresults.txt", 'w+')
z12 = open("rad30_singlebase_insertions_homopolymer_strandresults.txt", 'w+')




WT_counts = ContextCounts(d, e)
rad16_counts = ContextCounts(d1, e1)
rad26_counts = ContextCounts(d2, e2)
ContextCounts(d3, e3)
WT_HPdel_counts = ContextCounts(d4, e4)
rad16_HPdel_counts = ContextCounts(d5, e5)
rad26_HPdel_counts = ContextCounts(d6, e6)
WT_HPins_counts = ContextCounts(d7, e7)
rad16_HPins_counts = ContextCounts(d8, e8)
rad26_HPins_counts = ContextCounts(d9, e9)
ContextCounts(d10, e10)
ContextCounts(d11, e11)
WT_countsB = ContextCounts(h, i)
rad16_countsB = ContextCounts(h1, i1)
rad26_countsB = ContextCounts(h2, i2)
ContextCounts(h3, i3)
WT_HPdel_countsB = ContextCounts(h4, i4)
rad16_HPdel_countsB = ContextCounts(h5, i5)
rad26_HPdel_countsB = ContextCounts(h6, i6)
WT_HPins_countsB = ContextCounts(h7, i7)
rad16_HPins_countsB = ContextCounts(h8, i8)
rad26_HPins_countsB = ContextCounts(h9, i9)
WT_counts2 = ContextCounts(z1, z3)
rad30_counts = ContextCounts(z2, z4)
WT_HPdel_counts2 = ContextCounts(z5, z7)
rad30_HPdel_counts = ContextCounts(z6, z8)
WT_HPins_counts2 = ContextCounts(z9, z11)
rad30_HPins_counts = ContextCounts(z10, z12)

d.close()
d1.close()
d2.close()
d3.close()
d4.close()
d5.close()
d6.close()
d7.close()
d8.close()
d9.close()
e.close()
e1.close()
e2.close()
e3.close()
e4.close()
e5.close()
e6.close()
e7.close()
e8.close()
e9.close()
h.close()
h1.close()
h2.close()
h3.close()
i.close()
i1.close()
i2.close()
i3.close()
h4.close()
h5.close()
h6.close()
h7.close()
h8.close()
h9.close()
i4.close()
i5.close()
i6.close()
i7.close()
i8.close()
i9.close()
z1.close()
z2.close()
z3.close()
z4.close()
z5.close()
z6.close()
z7.close()
z8.close()
z9.close()
z10.close()
z11.close()
z12.close()


j.close()
j1.close()
j2.close()
j3.close()
k.close()
k1.close()
k2.close()
k3.close()
m.close()
m1.close()
m2.close()
m3.close()
n.close()
n1.close()
n2.close()
n3.close()
p.close()
p1.close()
p2.close()
p3.close()
q.close()
q1.close()
q2.close()
q3.close()
r.close()
r1.close()
r2.close()
r3.close()
s.close()
s1.close()
s2.close()
s3.close()

ContextCounts(t, u)
ContextCounts(t1, u1)
ContextCounts(t2, u2)
ContextCounts(t3, u3)
ContextCounts(v, w)
ContextCounts(v1, w1)
ContextCounts(v2, w2)
ContextCounts(v3, w3)

t.close()
t1.close()
t2.close()
t3.close()
u.close()
u1.close()
u2.close()
u3.close()
v.close()
v1.close()
v2.close()
v3.close()
w.close()
w1.close()
w2.close()
w3.close()


file = open("Iyer_saccer3_genecoord_noverlap_expand_strand_trinucfreq.txt")
#Finds the number of T and C bases on the NTS and TS
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
    return AT_NTS, AT_TS, CG_NTS, CG_TS

strand_numbers = FindTC(file)

NTS_mutations = strand_numbers[0] + strand_numbers[2]
TS_mutations = strand_numbers[1] + strand_numbers[3]

file.close()

with open('Iyer_saccer3_genecoord_noverlap.bed') as f2:
    data3 = f2.read()
lines3 = [s.strip().split() for s in data3.splitlines()]
chr3 = []
pos1_3 = []
pos2_3 = []
strand_data3 = []
for line in lines3:
    chr3.append(line[0])
    pos1_3.append(int(line[1]))
    pos2_3.append(int(line[2]))
    strand_data3.append(line[5])



TS_base = 0
NTS_base = 0
for x in range (0, len(chr3)):
    if chr3[x] == 'chrI':
        sequence = sequence1
    if chr3[x] == 'chrII':
        sequence = sequence2
    if chr3[x] == 'chrIII':
        sequence = sequence3
    if chr3[x] == 'chrIV':
        sequence = sequence4
    if chr3[x] == 'chrV':
        sequence = sequence5
    if chr3[x] == 'chrVI':
        sequence = sequence6
    if chr3[x] == 'chrVII':
        sequence = sequence7
    if chr3[x] == 'chrVIII':
        sequence = sequence8
    if chr3[x] == 'chrIX':
        sequence = sequence9
    if chr3[x] == 'chrX':
        sequence = sequence10
    if chr3[x] == 'chrXI':
        sequence = sequence11
    if chr3[x] == 'chrXII':
        sequence = sequence12
    if chr3[x] == 'chrXIII':
        sequence = sequence13
    if chr3[x] == 'chrXIV':
        sequence = sequence14
    if chr3[x] == 'chrXV':
        sequence = sequence15
    if chr3[x] == 'chrXVI':
        sequence = sequence16
    #for i in range(0, (pos2_3[x] - pos1_3[x]) + 1):
    b = pos1_3[x]
    while b <= (pos2_3[x]) - 1:
        #for b in range(0, len(sequence) - 1):
            polymer_range = {}
            length = 1
            if sequence[b] != 'N':# and b not in coordinates:
                length = 1
                i = 1
                polymer_range[0] = b
                while sequence[b] == sequence[b + i]:
                    if (b + i) <= (pos2_3[x]) - 1:
                        polymer_range[1] = b + i
                        length = length + 1
                        i = i + 1
                    else:
                        break
                   
              

                if length >= 4:
                
                    if strand_data3[x] == '+':
                                if sequence[b] == 'T' or sequence[b] == 'C':
                                        NTS_base = NTS_base + length
                                        
                                if sequence[b] == 'A' or sequence[b] == 'G':
                                        TS_base = TS_base + length
                    elif strand_data3[x] == '-':
                                if sequence[b] == 'A' or sequence[b] == 'G':
                                        NTS_base = NTS_base + length
                                if sequence[b] == 'T' or sequence[b] == 'C':
                                        TS_base = TS_base + length
                b = b + length
print(NTS_base)
print(TS_base)




#Caculates normalized percentages on the NTS and TS based on sequence context
def FindNormalizedRatio(NTS_count, TS_count, NTS_numbers, TS_numbers):
    mutation_ratio = TS_numbers/NTS_numbers
    base_ratio = NTS_count/TS_count
    ratio = mutation_ratio * base_ratio
    NTS_percentage = (ratio * 100)/ (ratio + 1)
    TS_percentage = 100/ (ratio + 1)
    return NTS_percentage, TS_percentage


FindNormalizedRatio(WT_counts[0], WT_counts[1], NTS_mutations, TS_mutations)
FindNormalizedRatio(rad16_counts[0], rad16_counts[1], NTS_mutations, TS_mutations)
FindNormalizedRatio(rad26_counts[0], rad26_counts[1], NTS_mutations, TS_mutations)

FindNormalizedRatio(WT_countsB[0], WT_countsB[1], NTS_mutations, TS_mutations)
FindNormalizedRatio(rad16_countsB[0], rad16_countsB[1], NTS_mutations, TS_mutations)
FindNormalizedRatio(rad26_countsB[0], rad26_countsB[1], NTS_mutations, TS_mutations)

FindNormalizedRatio(WT_HPdel_counts[0] + WT_HPins_counts[0], WT_HPdel_counts[1] + WT_HPins_counts[1], NTS_base, TS_base)
FindNormalizedRatio(rad16_HPdel_counts[0] + rad16_HPins_counts[0], rad16_HPdel_counts[1] + rad16_HPins_counts[1], NTS_base, TS_base)
FindNormalizedRatio(rad26_HPdel_counts[0] + rad26_HPins_counts[0], rad26_HPdel_counts[1] + rad26_HPins_counts[1], NTS_base, TS_base)


# FindNormalizedRatio(WT_HPins_counts[0], WT_HPins_counts[1], NTS_base, TS_base)
# FindNormalizedRatio(rad16_HPins_counts[0], rad16_HPins_counts[1], NTS_base, TS_base)
# FindNormalizedRatio(rad26_HPins_counts[0], rad26_HPins_counts[1], NTS_base, TS_base)

FindNormalizedRatio(WT_HPdel_countsB[0] + WT_HPins_countsB[0], WT_HPdel_countsB[1] + WT_HPins_countsB[1], NTS_base, TS_base)
FindNormalizedRatio(rad16_HPdel_countsB[0] + rad16_HPins_countsB[0], rad16_HPdel_countsB[1] + rad16_HPins_countsB[1], NTS_base, TS_base)
FindNormalizedRatio(rad26_HPdel_countsB[0] + rad26_HPins_countsB[0], rad26_HPdel_countsB[1] + rad26_HPins_countsB[1], NTS_base, TS_base)


# FindNormalizedRatio(WT_HPins_countsB[0], WT_HPins_countsB[1], NTS_base, TS_base)
# FindNormalizedRatio(rad16_HPins_countsB[0], rad16_HPins_countsB[1], NTS_base, TS_base)
# FindNormalizedRatio(rad26_HPins_countsB[0], rad26_HPins_countsB[1], NTS_base, TS_base)

FindNormalizedRatio(WT_counts2[0], WT_counts2[1], NTS_mutations, TS_mutations)
FindNormalizedRatio(rad30_counts[0], rad30_counts[1], NTS_mutations, TS_mutations)

FindNormalizedRatio(WT_HPdel_counts2[0] + WT_HPins_counts2[0], WT_HPdel_counts2[1] + WT_HPins_counts2[1], NTS_base, TS_base)
FindNormalizedRatio(rad30_HPdel_counts[0] + rad30_HPins_counts[0], rad30_HPdel_counts[1] + rad30_HPins_counts[1], NTS_base, TS_base)







