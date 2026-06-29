import os
import pandas as pd
workingDirectory = os.path.dirname(__file__)

in_or_del= "del"
name = "WT"

BedFile = "UVC-BedForFASTA/UVC-" + name +"_"+ in_or_del+".bed"
output_P = "UVC/UVC-BedForFASTA/UVC-" + name +"_"+ in_or_del+"_+.bed"
output_M = "UVC/UVC-BedForFASTA/UVC-"+ name +"_"+ in_or_del+"_-.bed"

fasta = {
    "Chromosome" : [],
    "Start(0)" : [],
    "Start(1)" : [],
    "Isolate" : [],
    "HP Count" : [],
    "Strand" : []
}


def parseLine(line, dataframe):

        newRow = { "Chromosome" : line[0], "Start(0)" : line[1], "Start(1)" : line[2], "Isolate" : line[3],
                 "HP Count" : line[4], "Strand" : line[5]}
     
        newRowDf = pd.DataFrame([newRow], columns = dataframe.columns) #make the new row in dataframe format
        dataframe = pd.concat([dataframe, newRowDf], ignore_index = True) #combine the two dataframes
     
        return dataframe



with open(os.path.join(workingDirectory, BedFile), "r") as file:
                dfPlus = pd.DataFrame(fasta) #make dictionary into dataframe
                dfMinus = pd.DataFrame(fasta)

                for line in file: #iterate through each line
                        if line.strip(): #if line not empty:
                                line = line.split()
                                if line[5] == '+':
                                    dfPlus = parseLine(line, dfPlus)
                                else:
                                    dfMinus = parseLine(line, dfMinus)


dfPlus.to_csv(output_P, sep = '\t', header = False, index = False, float_format= "%.0f")
dfMinus.to_csv(output_M, sep = '\t', header = False, index = False, float_format= "%.0f")

