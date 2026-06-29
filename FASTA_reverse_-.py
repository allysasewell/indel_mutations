import os
import numpy as np
workingDirectory = os.path.dirname(__file__)

name = "rad16"
in_or_del = "ins"

#UVC
#FASTAFile = "UVC/Combined_Clusters/nonCluster_fasta/"+in_or_del+ " fasta/UVC-" + name + "_" + in_or_del+"_NONclusters.fasta"
#output = "UVC/Combined_Clusters/nonCluster_fasta/"+in_or_del+ " fasta/UVC-" + name + "_" + in_or_del+"_NONclusters_REVERSED.fasta"
#FASTAFile = "UVC/Combined_Clusters/cluster_fasta/"+in_or_del+ " fasta/UVC-" + name + "_" + in_or_del+"_clusters.fasta"
#output = "UVC/Combined_Clusters/cluster_fasta/"+in_or_del+ " fasta/UVC-" + name + "_" + in_or_del+"_clusters_REVERSED.fasta"

#UVB
#FASTAFile = "UVB/UVB-Clusters/non_cluster/"+in_or_del+ " fasta/UVB-" + name + "_" + in_or_del + "_NONclusters.fasta"
#output = "UVB/UVB-Clusters/non_cluster/"+in_or_del+ " fasta/UVB-" + name + "_" + in_or_del+"_NONclusters_REVERSED.fasta"
FASTAFile = "UVB/UVB-Clusters/"+in_or_del+ " fasta/UVB-" + name + "_" + in_or_del + "_clusters.fasta"
output = "UVB/UVB-Clusters/"+in_or_del+ " fasta/UVB-" + name + "_" + in_or_del+"_clusters_REVERSED.fasta"



fasta = []
counter = 0
with open(os.path.join(workingDirectory, FASTAFile), "r") as file:
    for line in file:
        if line.strip(): #if line not empty:
                line = line.split()
                
                if ">" in line[0]: #if title, add to file
                    fasta += line
                else:
                     if "-" in fasta[counter - 1]: #if a minus strand
                        line = line[0][::-1] #reverse to make reverse complement just complement
                        fasta.append(line)
                     else:
                          fasta += line
        counter += 1
             

                
with open(os.path.join(workingDirectory, output), "w") as file2:
     for item in fasta:
          file2.write(item + "\n")