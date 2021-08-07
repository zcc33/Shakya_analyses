import os
import numpy as np
import scipy.stats
import random
import math
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

os.chdir("/Users/yancyliao/Dropbox/bioinformatics/covariation_analysis")
positions_file = "proximate_positions"
seq_file = "results"
pdb_align_file = "pdb_align"
trusted_cut = 130.6
hmm_length = 228

def extractAlign(text):
    seq_id = text[0].split()[1]
    
    if len(text) < 4:
        return "", seq_id
        
    seq_score = float(text[3].split()[2])
    if seq_score < trusted_cut:
        return "", seq_id
        
    hmm_from = int(text[3].split()[6])
    hmm_to = int(text[3].split()[7])
    
    hmm_alignment = ""
    hmm_seq = ""
    seq = ""
    
    for i in range(hmm_from-1):
        hmm_alignment += "-"
    
    for index in range(len(text)):
        line = text[index].split()
        
        if len(line) > 2 and line[1] == "domain" and line[2] == "2":
            break
            
        if len(line) > 0 and line[0] == seq_id:
            hmm_line = text[index-2].split()
            hmm_seq = hmm_seq + hmm_line[2]
            seq = seq + line[2]
    
    for i in range(len(seq)):
        if seq[i] == "-":
            hmm_alignment += "-"
        elif seq[i].isupper():
            hmm_alignment += seq[i]
    
    
    for i in range(hmm_length - hmm_to):
        hmm_alignment += "-"

    return hmm_alignment, seq_id

def get_align(seq_file):
    alignments = []
    indices = []
    ids = []
    with open(seq_file, "r") as seq_file:
        s = seq_file.read().splitlines()
        for index in range(len(s)):
            line = s[index].split()
            if len(line) > 0 and ">>" == line[0]:
                indices.append(index)
        indices.append(len(s)-1)
        
        for i in range(len(indices)-1):
            first_index = indices[i]
            second_index = indices[i+1]
            new_alignment, new_id = extractAlign(s[first_index:second_index])
            if not new_alignment == "":
                alignments.append(new_alignment)
                ids.append(new_id)
            
    for i in range(len(alignments)):
        if not len(alignments[i]) == hmm_length:
            print "Warning, alignment length not correct."
            print i
    return alignments, indices, ids

def get_pos():
    position1 = []
    position2 = []
    with open(positions_file, "r") as pos:
        for line in pos:
            if line.split()[0] == "Common":
                position1.append(int(line.split()[1]))
                position2.append(int(line.split()[2]))
                
                
    #translate the positions from pdb files to the ones that are valid for the hmms
    with open(pdb_align_file, "r") as pdb_file:
        pdb_file = pdb_file.readlines()
        alignment = pdb_file[2].split()[1] + pdb_file[7].split()[1]
        alignment = alignment.replace('\n', '')
        valid_pos = []
        count = 0
        for i in range(len(alignment)):
            if alignment[i].isupper():
                count+=1
                valid_pos.append(count)
            elif alignment[i]=="-":
                count+=1
            else:
                valid_pos.append(-1)
    
    for i in range(len(position1)):
        pos1_val = position1[i]
        pos2_val = position2[i]
        if valid_pos[pos1_val-1] == -1 or valid_pos[pos2_val-1] == -1:
            position1[i] = -1
            position2[i] = -1
        else:
            position1[i], position2[i] = valid_pos[pos1_val-1], valid_pos[pos2_val-1]
    
    return position1, position2

def make_contingency(vec1, vec2):
    count1 = 0
    count2 = 0
    count3 = 0
    count4 = 0
    for i in range(len(vec1)):
        if vec1[i] == 1 and vec2[i] == 1:
            count1 +=1
        elif vec1[i] == 1 and vec2[i] == 2:
            count2 +=1
        elif vec1[i] == 2 and vec2[i] == 1:
            count3+=1
        else:
            count4+=1
    if count1==0 or count2==0 or count3==0 or count4 ==0:
        return -1
    else:
        return np.array([[count1, count2], [count3, count4]])

def chi_sq(vec1, vec2):
    count1 = 0
    count2 = 0
    count3 = 0
    count4 = 0
    for i in range(len(vec1)):
        if vec1[i] == 1 and vec2[i] == 1:
            count1 +=1
        elif vec1[i] == 1 and vec2[i] == 2:
            count2 +=1
        elif vec1[i] == 2 and vec2[i] == 1:
            count3+=1
        else:
            count4+=1
    if (count1==0 and count3==0) or (count1 ==0 and count2 == 0) or (count2 == 0 and count4 == 0) or (count3 ==0 and count4 == 0):
        return -1
    else:
        obs = np.array([[count1, count2], [count3, count4]])
        return scipy.stats.chi2_contingency(obs, lambda_="log-likelihood")[1]
    
def get_col(alignments, pos1, pos2):
    vec1 = []
    vec2 = []
    for i in range(len(alignments)):
        vec1.append(alignments[i][pos1-1])
        vec2.append(alignments[i][pos2-1])
    pruned_vec1, pruned_vec2 = [], []
    for i in range(len(vec1)):
        if not (vec1[i] == '-' or vec2[i] == '-'):
            pruned_vec1.append(vec1[i])
            pruned_vec2.append(vec2[i])
    return pruned_vec1, pruned_vec2
    
def get_col2(alignments, pos1, pos2):
    vec1 = []
    vec2 = []
    for i in range(len(alignments)):
        if not (alignments[i][pos1-1] == '-' or alignments[i][pos2-1] == '-'):
            vec1.append(alignments[i][pos1-1])
            vec2.append(alignments[i][pos2-1])
    return vec1, vec2
    
def convert_binary(vec1, vec2, cat1):
    bin1 = []
    bin2 = []
    for i in range(len(vec1)):
        if vec1[i] in cat1:
            bin1.append(1)
        else:
            bin1.append(2)
        if vec2[i] in cat1:
            bin2.append(1)
        else:
            bin2.append(2)
    return bin1, bin2
    
#values is an array of chi square p-values
def summarize(values):
    count = 0
    
    thresh1=0.01
    thresh2=10**-20
    thresh3=10**-75
    count1 = 0
    count2=0
    count3=0
    
    for i in range(len(values)):
        if values[i] == -1:
            continue
        else:
            count+=1
            
        if values[i] < thresh1:
            count1+=1
        if values[i] < thresh2:
            count2+=1
        if values[i] < thresh3:
            count3+=1
    print count, count1/float(count), count2/float(count), count3/float(count)
            

def adj_rand_pos(alignments, position1, position2):
    #these are the positions next to each other
    position1a = []
    position2a = []
    for i in range(len(alignments[0])-1):
        position1a.append(i)
        position2a.append(i+1)
    
    #these are the positions randomly chosen
    position1b = []
    position2b = []
    for i in range(len(alignments[0])-1):
        position1b.append(random.randint(0,len(alignments[0])))
        position2b.append(random.randint(0,len(alignments[0])))
    
    #make sure there is no overlap with the proximate positions
    bad_pos_a = []
    bad_pos_b = []
    curiosity = []
    for i in range(len(position1a)):
        for j in range(len(position1)):
            if position1a[i] == position1[j] and position2a[i] == position2[j]:
                bad_pos_a.append(i)
                curiosity.append(j)
            if position1b[i] == position1[j] and position2b[i] == position2[j]:
                bad_pos_b.append(i)
    
    for i in range(len(bad_pos_a)):
        index = bad_pos_a[i]
        position1a = position1a[:index] + position1a[index+1:]
        position2a = position2a[:index] + position2a[index+1:]
    for i in range(len(bad_pos_b)):
        index = bad_pos_b[i]
        position1b = position1b[:index] + position1b[index+1:]
        position2b = position2b[:index] + position2b[index+1:]
    return position1a, position2a, position1b, position2b

def calc_chisq(alignments, position1, position2, cat1):
    chisq_values = []
    for j in range(len(position1)):
        pos1 = position1[j]
        pos2 =position2[j]
        
        if pos1 == -1:
            continue
        
        vec1 = []
        vec2 = []
        for i in range(len(alignments)):
            if (not alignments[i][pos1-1] == "-") and (not alignments[i][pos2-1] == "-"):
                if alignments[i][pos1-1] in cat1:
                    vec1.append(1)
                else:
                    vec1.append(2)
                if alignments[i][pos2-1] in cat1:
                    vec2.append(1)
                else:
                    vec2.append(2)
    
        chisq_values.append(chi_sq(vec1, vec2))
        
        if chi_sq(vec1, vec2) == 0:
            print str(position1[j]) + " " + str(position2[j]) + "\n"
        
    return chisq_values


def rand_cat():
    master = "ARNDCQEGHILKMFPSTWYV"
    cat = ""
    for i in range(len(master)):
        if random.random() < 0.5:
            cat = cat + master[i]
    return cat

def get_bottom(alignments, position1, position2, chisq_values, cat1, output):
    file = open(output, "a")
    d = chisq_values
    for i in range(len(d)):
        if d[i] == -1:
            d[i] = 1
    newchi = np.array(d)
    bottom5indices = newchi.argsort()[:5]
    pos1_clean = [i for i in position1 if not i == -1]
    pos2_clean = [i for i in position2 if not i == -1]
    file.write("\n" + cat1 + "\n")
    for i in bottom5indices:
        c,d = get_col2(alignments, pos1_clean[i], pos2_clean[i])
        c,d = convert_binary(c,d,cat1)
        file.write(str(make_contingency(c,d)))
        file.write("\n")
        file.write(str(chisq_values[i]) + "\n")
        file.write(str(pos1_clean[i]) + " " + str(pos2_clean[i]) + "\n")
    file.close()

def log_chisq(x):
    d = []
    for i in x:
        if not i ==-1:
            d.append(i)
    x=d
    for i in range(len(x)):
        if not x[i] == 0:
            x[i] = math.log10(x[i])
    return x



def main2():
    alignments, indices, ids = get_align(seq_file)
    position1, position2 = get_pos()
    position1adj, position2adj, position1rand, position2rand = adj_rand_pos(alignments, position1, position2)
    
    #not large, not tiny, polar, nonpolar, large, medium, negative, positive
    catlist = ["ACDGILMNPSTV", "DEFHIKLMNPQRTVWY", "CDEHKNQRST", "AFILMPVW", "EFHKQRWY", "DILMNPTV", "DE", "KRH"]
    
    for i in range(10):
        catlist.append(rand_cat())
    
    output = "chisq_results.txt"
    file = open(output, "w")
    file.write("starting" + "\n")
    file.close()
    
    for cat in catlist:
        chisq_values = calc_chisq(alignments, position1,position2, cat)
        get_bottom(alignments, position1, position2, chisq_values, cat, output)


def main1():
    alignments, indices, ids = get_align(seq_file)
    position1, position2 = get_pos()
    position1adj, position2adj, position1rand, position2rand = adj_rand_pos(alignments, position1, position2)
    num_cat = 250
    catlist = []
    for i in range(num_cat):
        catlist.append(rand_cat())
    
    min_values = []
    for i in range(len(catlist)):
        print i
        cat = catlist[i]
        chisq_values = calc_chisq(alignments, position1,position2, cat)
        min_val = 99
        for i in range(len(chisq_values)):
            if chisq_values[i] < min_val and not chisq_values[i] ==-1:
                min_val = chisq_values[i]
        min_values.append(min_val)
    
    plot_vals = log_chisq(min_values)
    for i in range(len(plot_vals)):
        if plot_vals[i] == 0:
            plot_vals[i] = -350
        
    n, bins, patches = plt.hist(plot_vals, bins = [-350, -300, -250, -200, -150, -100, -50, 0], rwidth=0.9)
    plt.xlabel("Log p-values (base 10)")
    plt.ylabel("Counts")
    plt.title("Lowest chi-square pvalues among proximate positions across "+ str(num_cat)+ " random categories")
    
    ax = plt.gca()
    ax.set_xticks([-350, -300, -250, -200, -150, -100, -50, -0])
    ax.set_xticklabels(["TRUE ZERO", -300, -250, -200, -150, -100, -50, -0])
    
    plt.show()
    print min_values
    for i in range(len(chisq_values)):
        if chisq_values[i] == 0:
            print position1[i]
            print "\n"
        

def main3():
    alignments, indices, ids = get_align(seq_file)
    position1, position2 = get_pos()
    position1adj, position2adj, position1rand, position2rand = adj_rand_pos(alignments, position1, position2)

    position1 = []
    position2 = []
    for i in range(len(alignments[0])-1):
        for j in range(i+1, len(alignments[0])):
            position1.append(i+1)
            position2.append(j+1)
        print i+1

    catlist = ["ADFIKMPQTW"]
    #not large, not tiny, polar, nonpolar, large, medium, negative, positive, bad
    #catlist = ["ACDGILMNPSTV", "DEFHIKLMNPQRTVWY", "CDEHKNQRST", "AFILMPVW", "EFHKQRWY", "DILMNPTV", "DE", "KRH", "ADFIKMPQTW"]
    
    #catlist = []
    #for i in range(1):
    #    catlist.append(rand_cat())
    
    cat = catlist[0]
    chisq_values = calc_chisq(alignments, position1,position2, cat)
    
    plot_vals = log_chisq(chisq_values)
    for i in range(len(plot_vals)):
        if plot_vals[i] == 0:
            plot_vals[i] = -350
            
    n, bins, patches = plt.hist(plot_vals, bins = [-350, -300, -250, -200, -150, -100, -50, 0], rwidth=0.9)
    plt.xlabel("Log p-values (base 10)")
    plt.ylabel("Counts")
    plt.title("Chi-square p-values for the category "+ catlist[0]+ " across "+ str(len(position1))+ " proximate positions")
    
    ax = plt.gca()
    ax.set_xticks([-350, -300, -250, -200, -150, -100, -50, -0])
    ax.set_xticklabels(["TRUE ZERO", -300, -250, -200, -150, -100, -50, -0])
    
    plt.show()
    print chisq_values
    print n
    print alignments[0]

def main():
    alignments, indices, ids = get_align(seq_file)
    position1, position2 = get_pos()
    position1adj, position2adj, position1rand, position2rand = adj_rand_pos(alignments, position1, position2)
    num_cat = 250
    catlist = ["ACDGILMNPSTV", "DEFHIKLMNPQRTVWY", "CDEHKNQRST", "AFILMPVW", "EFHKQRWY", "DILMNPTV", "DE", "KRH"]
    
    #catlist=[]
    #for i in range(num_cat):
    #    catlist.append(rand_cat())
    
    num_zeros = []
    for i in range(len(catlist)):
        print i
        cat = catlist[i]
        chisq_values = calc_chisq(alignments, position1,position2, cat)
        count = 0
        for i in range(len(chisq_values)):
            if chisq_values[i] ==0:
                count+=1
        num_zeros.append(count)
    
    n, bins, patches = plt.hist(num_zeros, rwidth=0.9)
    plt.xlabel("Number of zeros")
    plt.ylabel("Counts of categories")
    plt.title("Number of zeros among the p-values for 8 meaningful categories")
    
    plt.show()
    print num_zeros

#chisq_valuesadj =calc_chisq(position1adj, position2adj)
#chisq_valuesrand = calc_chisq(position1rand, position2rand)

#summarize(chisq_values)
#summarize(chisq_valuesadj)
#summarize(chisq_valuesrand)

def get_low(x):
    min_val = 99
    for i in range(len(x)):
        if x[i] < min_val:
            min_val = x[i]
    return min_val

def get_zeros(x):
    count = 0
    for i in range(len(x)):
        if x[i] == 0:
            count+=1
    return count


if __name__ == "__main__": main()
