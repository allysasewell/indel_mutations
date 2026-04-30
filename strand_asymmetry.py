
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




d = open("WT_singlebase_deletions_TriNuc_TrxStrandBins.txt")
d1 = open("rad16_singlebase_deletions_TriNuc_TrxStrandBins.txt")
d2 = open("rad26_singlebase_deletions_TriNuc_TrxStrandBins.txt")
d3 = open("rad30_singlebase_deletions_TriNuc_TrxStrandBins.txt")
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


h = open("UVB_WT_singlebase_deletions_TriNuc_TrxStrandBins.txt")
h1 = open("UVB_rad16_singlebase_deletions_TriNuc_TrxStrandBins.txt")
h2 = open("UVB_rad26_singlebase_deletions_TriNuc_TrxStrandBins.txt")
h3 = open("UVB_rad30_singlebase_deletions_TriNuc_TrxStrandBins.txt")
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


z1 = open("WT_singlebase_substitutionsv2_Trinuc_TrxStrandBins.txt")
z2 = open("rad30_singlebase_substitutionsv2_Trinuc_TrxStrandBins.txt")
z3 = open("WT_singlebase_substitutionsv2_strandresults.txt", 'w+')
z4 = open("rad30_singlebase_substitutionsv2_strandresults.txt", 'w+')




ContextCounts(d, e)
ContextCounts(d1, e1)
ContextCounts(d2, e2)
ContextCounts(d3, e3)
ContextCounts(d4, e4)
ContextCounts(d5, e5)
ContextCounts(d6, e6)
ContextCounts(d7, e7)
ContextCounts(d8, e8)
ContextCounts(d9, e9)
ContextCounts(d10, e10)
ContextCounts(d11, e11)
ContextCounts(h, i)
ContextCounts(h1, i1)
ContextCounts(h2, i2)
ContextCounts(h3, i3)
ContextCounts(h4, i4)
ContextCounts(h5, i5)
ContextCounts(h6, i6)
ContextCounts(h7, i7)
ContextCounts(h8, i8)
ContextCounts(h9, i9)
ContextCounts(z1, z3)
ContextCounts(z2, z4)

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

#Caculates normalized percentages on the NTS and TS based on sequence context
def FindNormalizedRatio(NTS_count, TS_count, NTS_mutations, TS_mutations):
    mutation_ratio = NTS_mutations/TS_mutations
    base_ratio = TS_count/NTS_count
    ratio = mutation_ratio * base_ratio
    NTS_percentage = (ratio * 100)/ (ratio + 1)
    TS_percentage = 100/ (ratio + 1)
    return NTS_percentage, TS_percentage