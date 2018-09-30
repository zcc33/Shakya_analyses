#!/usr/bin/python

"""
This program takes as input the reference analyses files of a particular gene in the 
Shakya data set, as well as the mHMM summary file for that gene. Then it extracts out
all the mHMM positions that were deemed high in diversity score and F1 of precision
and accuracy. Using a weighted interval scheduling algorithm, this program weights each
position by an F1-score of the scaled diversity and F1-scores. The output is a list of
optimal non-overlapping positions, given that each position corresponds to an interval
of length 6.
"""
import os, sys, csv
from collections import Counter
genomes_folder = sys.path[0] + "/RplB/reference_analyses/"
csv_folder = sys.path[0] + "/RplB/reference_csvs/"
summary_file = sys.path[0] + "/RplB/TIGR01171_umdB_mHMM_summary"

"""
Gets a list of the feasible_positions according to the reference_analyses files
returns feasible_positions (list)
"""
def getFeasiblePositions():
    feasible_positions = []
    genome_file = os.listdir(genomes_folder)[1]
    with open(genomes_folder+ genome_file, 'r') as handle:
        for line in handle:
            if (len(line) < 2):
                break
            else:
                splitted = line.split()
                if(not(splitted[0][-1] == "x")):
                    feasible_positions.append(int(splitted[0]))
    handle.close()
    return feasible_positions

"""
Gets the 
"""
def getPositionScores(feasible_positions):
    position_scores = {}
    with open(summary_file, 'r') as handle:
        for line in handle:
            splitted = line.split()
            try:
                if(int(splitted[0]) in feasible_positions):
                    scaled_f1score = ((0.01*float(splitted[12])) - 0.99)/0.01
                    scaled_div = ((0.01*float(splitted[14]))-0.9)/0.1
                    score = 2*(scaled_f1score*scaled_div)/(scaled_f1score+scaled_div)
                    position_scores[int(splitted[0])] = score
            except ValueError:
                continue
    handle.close()
    return position_scores

"""
we
"""
def getOptimalPositions(feasible_positions, position_scores):
    p = []
    for i in range(len(feasible_positions)):
        if(feasible_positions[i] <= 6):
            p.append(0)
        else:
            best = 0
            for j in range(p[i-1], i):
                if feasible_positions[j] <= feasible_positions[i]-6:
                    best = j
            p.append(best)
    
    opt = []
    opt.append(position_scores[feasible_positions[0]])
    for i in range(1, len(feasible_positions)):
        opt.append(max(position_scores[feasible_positions[i]] + opt[p[i]], opt[i-1]))
        
    positions = []
    i = len(feasible_positions)-1
    while(i >= 0):
        if(position_scores[feasible_positions[i]] + opt[p[i]] >= opt[i-1]):
            positions.append(feasible_positions[i])
            i = p[i]
        else:
            i=i-1
            
    return list(reversed(positions))
    

def turn_to_csv(opt_pos):
    os.chdir(genomes_folder)
    print(os.listdir(genomes_folder))
    for filename in os.listdir(genomes_folder):
        if(filename[0] == '.' or filename.endswith('.csv')):
            continue
        
        with open(filename, 'r') as handle:
            depth = []
            for line in handle:
                if (len(line) < 2):
                    break
                    
                splitted = line.split()
                if(not(splitted[0][-1] == "x")):
                    if(int(splitted[0]) in opt_pos and splitted[1] != 'MISSING' and splitted[1] != 'unassigned'):
                        
                         if(len(splitted) ==2):
                            if(splitted[1] == "0"):
                                depth.append(splitted[1])
                         elif(len(splitted) >2):
                            depth.append(splitted[2])
                        
            counts = Counter(depth).items()
            
            combined = []
            for i in range(len(counts)):
                combined.append([counts[i][0], str(counts[i][1])])
            
            with open(csv_folder + filename + ".csv", 'wb') as myfile:
                writer = csv.writer(myfile, delimiter=',')
                for line in combined:
                    writer.writerow(line)
            
            myfile.close()
        handle.close()
    
if __name__ == "__main__":
    feasible_positions = getFeasiblePositions()
    position_scores = getPositionScores(feasible_positions)
    opt_pos=getOptimalPositions(feasible_positions, position_scores)
    turn_to_csv(opt_pos)
    