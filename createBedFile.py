import os
import pandas as pd
workingDirectory = os.path.dirname(__file__)

UVYeast = {
        "Chromosome" : [],        
        "Start(0)" : [],
        "Start(1)" : [],
        "Original" : [],
        "Indel" : [],
        "Number" : [],
        "Isolate" : [],
        "Strand" : []
}

wtDF = pd.DataFrame(UVYeast)
rad16DF = pd.DataFrame(UVYeast)
rad26DF = pd.DataFrame(UVYeast)
rad30DF = pd.DataFrame(UVYeast)

def makeRow(line):
    newRow = { "Chromosome" : "chr" +  line[6],  
                "Start(0)" : str(int(line[8]) - 1), "Start(1)" : line[8], "Original" : line[12],
                "Indel" : line[13], "Number" : line[11],  "Isolate" : line[0]}
    newRow = findStrand(newRow)

    return newRow

def findStrand(row): #T, C or TT, TC, CC, CT on + strand. Otherwise on -
    if row["Number"] == "1":
        if row["Original"] != "-":
            if row["Original"] == "C" or row["Original"] == "T":
                row["Strand"] = "+"
            else:
                row["Strand"] = "-"
        else:
            if row["Indel"] == "C" or row["Indel"] == "T":
                row["Strand"] = "+"
            else:
                row["Strand"] = "-"
    else:
        row["Strand"] = "+"
        '''if "-" not in row["Original"]:
            if "CC" in row["Original"] or "TC" in row["Original"] or "CT" in row["Original"] or "TT" in row["Original"]:
                row["Strand"] = "+"
            else:
                 row["Strand"] = "-"
        elif  "CC" in row["Indel"] or "TC" in row["Indel"] or "CT" in row["Indel"] or "TT" in row["Indel"]:
            row["Strand"] = "+"
        else:
                row["Strand"] = "-"'''

    return row


with open(os.path.join(workingDirectory,"UV_yeast_Merged_SNP_DIPs_sorted_anz2_tandem.txt")) as f:
    data = f.read()
    lines = [s.strip().split() for s in data.splitlines()]

    for line in lines:
        if 'Mito' not in line[6] and 'ySR128' not in line[3] and 'dRad16_25' not in line[1]:
            if line[10] == "DIP" or line[10] == "Complex DIP":
                if line[4] == "WT":
                    #if strainName == "WT":
                        newRow = makeRow(line)

                        newRowDF = pd.DataFrame([newRow], columns = wtDF.columns)
                        wtDF = pd.concat([wtDF, newRowDF], ignore_index = True)

                elif line[4] == "rad16D/rad16D":
                        #if strainName == "rad16":
                        newRow = makeRow(line)

                        newRowDF = pd.DataFrame([newRow], columns = wtDF.columns)
                        rad16DF = pd.concat([rad16DF, newRowDF], ignore_index = True)
                elif line[4] == "rad26D/rad26D":
                        #if strainName == "rad26":
                        newRow = makeRow(line)

                        newRowDF = pd.DataFrame([newRow], columns = wtDF.columns)
                        rad26DF = pd.concat([rad26DF, newRowDF], ignore_index = True)
                elif line[4] == "rad30D/rad30D":
                        #if strainName == "rad30" :
                        newRow = makeRow(line)

                        newRowDF = pd.DataFrame([newRow], columns = wtDF.columns)
                        rad30DF = pd.concat([rad30DF, newRowDF], ignore_index = True)


rad16DF.to_csv("UVC/Combined_Clusters/rad16_Combined_all.bed", sep = '\t', header = False, index = False)
rad26DF.to_csv("UVC/Combined_Clusters/rad26_Combined_all.bed", sep = '\t', header = False, index = False)
rad30DF.to_csv("UVC/Combined_Clusters/rad30_Combined_all.bed", sep = '\t', header = False, index = False)
wtDF.to_csv("UVC/Combined_Clusters/WT_Combined_all.bed", sep = '\t', header = False, index = False)