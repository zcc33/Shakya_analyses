import os
import sys
import math
import random

os.chdir("/Users/yancyliao/Dropbox/bioinformatics/covariation_analysis")
pdb_file = "1AW2.pdb"
output_file = "proximate_positions"
dist_threshold = 4
sample_index = "ABDEGHJK"

#get the distance between two points in 3d space
def distance_calc(array1, array2):
    return math.sqrt((array1[0] - array2[0])**2 + (array1[1] - array2[1])**2 + (array1[2] - array2[2])**2)

#find the minimum r-group atom distance for a set of atom coordinates between two amino acids
def min_rgroup_dist(index1, index2, pos, coordinates):
    min_dist = 999999
    indices1 = [x for x,val in enumerate(pos) if val==index1]
    indices2 = [x for x,val in enumerate(pos) if val==index2]
    for k in indices1[4:]:
        for l in indices2[4:]:
            coordinates1 = coordinates[k]
            coordinates2 = coordinates[l]
            distance = distance_calc(coordinates1, coordinates2)
            if distance < min_dist:
                min_dist = distance
    return min_dist
    
#results{"A"} = [(pos1, pos2, distance)...]
def main():
    results = {}
    for i in sample_index:
        results[i] = []
    with open(pdb_file, "r") as pdb:
        for sample in sample_index:
            pdb.seek(0)
            coordinates = []
            pos = []
            atom = []
            for line in pdb:
                line=line.split()
                if line[0] == "ATOM" and line[4] == sample:
                    coordinates.append([float(line[6]), float(line[7]), float(line[8])])
                    pos.append(int(line[5]))
                    atom.append(line[2])

            unique_pos = []
            for i in range(len(pos)):
                if not pos[i] in unique_pos:
                    unique_pos.append(pos[i])
    
            for i in range(len(unique_pos)):
                for j in range(i+1, len(unique_pos)):
                    pos1 = unique_pos[i]
                    pos2 = unique_pos[j]
                    distance = min_rgroup_dist(pos1, pos2, pos, coordinates)
                    if distance < dist_threshold:
                        results[sample].append((pos1, pos2, distance))
   
    positions = {}
    for sample in sample_index:
        positions[sample] = []
        for j in results[sample]:
            positions[sample].append(j[0:2])

    common = positions[sample_index[0]]
    for sample in sample_index:
        common = sorted(list(set(common) & set(positions[sample])))
        
    
    
    file = open(output_file, "w")
    file.write("Samp \t Pos1 \t Pos2 \t Dist \n")
    for i in common:
        file.write("Common \t" + str(i[0]) + "\t" + str(i[1]) + "\t N/A \n")
    file.write("####### Total: " + str(len(common))+"\n")
        
    for i in sample_index:
        for j in range(len(results[i])):
            file.write(i + "\t" + str(results[i][j][0]) + "\t" + str(results[i][j][1]) + "\t" +str(results[i][j][2]) + "\n")
        file.write("###### Total: " + str(len(results[i]))+"\n")
    file.close()
    
        
if __name__ == "__main__":
    main()