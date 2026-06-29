import os
import numpy as np
import pandas as pd

workingDirectory = os.path.dirname(__file__)
allChrTxt = ["chr1.txt", "chr2.txt", "chr3.txt", "chr4.txt", "chr5.txt", "chr6.txt", "chr7.txt", "chr8.txt",
             "chr9.txt", "chr10.txt", "chr11.txt", "chr12.txt", "chr13.txt", "chr14.txt", "chr15.txt",
             "chr16.txt"]
roman = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]
output = "genomeHPcount.txt"

A_count = T_count = G_count = C_count = 0
total_A = total_T = total_G = total_C = 0
HPLengthsAll = np.zeros(45, dtype =int)
data = (HPLengthsAll, 0, A_count, T_count, G_count, C_count, total_A, total_T, total_G, total_C)


def countHomopoly(dna, lenDNA, data):

        TempHPLengthsAll = np.zeros(45, dtype = int)
        base = 0
        HPcount = 0
        
        A_count = data[2]
        T_count = data[3]
        G_count = data[4]
        C_count = data[5]
        total_A = data[6]
        total_T = data[7]
        total_G = data[8]
        total_C = data[9]

        while base < lenDNA:

                HPlength = 1            #start with one because minimum is one
                HPbaseType = dna[base]
                while base + 1 < lenDNA and dna[base] == dna[base + 1]: #iterate while not at the end / is the same nucleotide
                        HPlength += 1
                        base += 1 

                if HPlength > 3:
                      HPcount+= 1
                      if HPbaseType == "A":
                        A_count += 1
                        total_A += HPlength
                      elif HPbaseType == "T":
                        T_count += 1
                        total_T += HPlength
                      elif HPbaseType == "G":
                        G_count += 1
                        total_G += HPlength
                      elif HPbaseType == "C":
                        C_count += 1
                        total_C += HPlength
                else:
                     if HPbaseType == "A":
                        total_A += HPlength
                     elif HPbaseType == "T":
                        total_T += HPlength
                     elif HPbaseType == "G":
                        total_G += HPlength
                     elif HPbaseType == "C":
                        total_C += HPlength
                     

                TempHPLengthsAll[HPlength - 1] += 1
                


                base += 1


        finalHPLengthsAll = np.add(data[0], TempHPLengthsAll) #add to previous values
        finalHPCount = data[1] + HPcount


        return (finalHPLengthsAll, finalHPCount, A_count, T_count, G_count, C_count, total_A, total_T, total_G, total_C)


for x in range(len(allChrTxt)):
    with open(os.path.join(workingDirectory, allChrTxt[x]), "r") as file:
        dna = file.read().replace('\n', '')
        dnaLen = len(dna)

        data = countHomopoly(dna, dnaLen, data)


df = pd.DataFrame(data[0]) #make the array from the tuple into a dataframe

with open(os.path.join(workingDirectory, output), "w") as w:
    df_str = df.to_string(header = False, index = False)
    w.write(df_str)
    w.write("\n\n\nHP Count: ")
    w.write(str(data[1]))
    
    w.write("\n\nA Count:\t")
    w.write(str(data[2]))
    w.write("\nT Count:\t")
    w.write(str(data[3]))  
    w.write("\nG Count:\t")  
    w.write(str(data[4]))
    w.write("\nC Count:\t")
    w.write(str(data[5]))
    w.write("\n\nA NT Count:\t")    
    w.write(str(data[6])) 
    w.write("\nT NT Count:\t")
    w.write(str(data[7]))
    w.write("\nG NT Count:\t")
    w.write(str(data[8]))
    w.write("\nC NT Count:\t")
    w.write(str(data[9]))