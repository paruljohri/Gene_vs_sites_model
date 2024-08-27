#This is to calculate LD stats on my own.
#D will be calculated so that it will always be p(11)-p(1)p(1). "1" can represent derived/minor allele.
#How to run:
#python statistics_LD_general_reps.py -fileExt ".ms" -sample_size 100 -regionLen 1000 -input_folder ${infolder}/${gamma}/${dom} -output_folder ${outfolder} -output_prefix ${fitness_model}_${gamma}_${dom}
from __future__ import print_function
import sys
#import pandas
import math
import argparse
import os

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-fileExt', dest = 'fileExt', action='store', nargs = 1, type = str, help = 'example: _masked.ms')
parser.add_argument('-sampleSize', dest = 'sampleSize', action='store', nargs = 1, default = 100, type = int, help = 'number of haploid genomes')#500 bp for small, 5000bp for big
parser.add_argument('-regionLen', dest = 'regionLen', action='store', nargs = 1, type = int, help = 'length in bp of region simulated')#Length of coding region simulated
parser.add_argument('-input_folder', dest = 'input_folder', action='store', nargs = 1, type = str, help = 'full path to folder with .ms files')
parser.add_argument('-output_folder', dest = 'output_folder', action='store', nargs = 1, type = str, help = 'full path to folder where you want to write the output')
parser.add_argument('-output_prefix', dest = 'output_prefix', action='store', nargs = 1, type = str, help = 'full path to output file')
args = parser.parse_args()
file_ext = args.fileExt[0]
chr_len =  args.regionLen[0]
sample_size = args.sampleSize[0]
#step_size = args.stepSize[0]/float(chr_len)
infolder = args.input_folder[0]
outfolder = args.output_folder[0]
prefix = args.output_prefix[0]

def read_ms_file(f_MS):
    l_POS = [] #list of float positions of SNPs
    d_MS = {} #integer position (0,1,2..) -> alleles
    for line in f_ms:
        line1 = line.strip('\n')
        if "positions" in line1:
            line2 = line1.split()
            i = 0
            for x in line2:
                if "position" not in x:
                    l_POS.append(x)
                    d_MS[str(i)] = ""
                    i += 1
        elif "//" not in line and "segsites" not in line:
            #print (d_MS)
            i = 0
            while i < len(line1):
                d_MS[str(i)] = d_MS[str(i)] + line1[i]
                i = i + 1
    print ("number of positions: " + str(len(l_POS)))
    print ("size of matrix: " + str(len(d_MS)))
    print(d_MS)
    return(l_POS, d_MS)

def convert_derived_to_minor(s_COL):
    if s_COL.count("1") <= int(sample_size/2):#derived allele is minor
        s_COL_minor = s_COL
    else:#derived allele is major
        s_COL_minor = ""
        for x in s_COL:
            if x == "0":
                s_COL_minor += "1"
            elif x == "1":
                s_COL_minor += "0"
    return(s_COL_minor)

def get_allele_freq(COL):
    s_AC = 0
    s_tot = 0
    for x in COL:
        if x == "0":
            s_tot += 1
        elif x == "1":
            s_tot += 1
            s_AC += 1
    #quick check:
    if s_tot != sample_size:
        print ("ERROR! Column size is not equal to sample size!")
    return(float(s_AC)/float(s_tot))

def calculate_D(COL1, COL2):
    #we are calculating p(11)-p(1)p(1)
    n = len(COL1)
    #calculate p(11)
    p11 = 0
    i = 0
    while i < n:
        if COL1[i] == "1" and COL2[i] == "1":
            p11 += 1
        i += 1
    #calcualte p(1) at COL1 and COL2
    p_1_locus1 = COL1.count("1")
    p_1_locus2 = COL2.count("1")
    #calculate D:
    STAT_D = (p11/float(n)) - (p_1_locus1/float(n))*(p_1_locus2/float(n))
    return(STAT_D)

def calculate_correlation_coefficient(COL1, COL2):
    q1 = get_allele_freq(COL1)
    q2 = get_allele_freq(COL2)
    s_D = calculate_D(COL1, COL2)
    denominator = float(q1)*(1.0-q1)*float(q2)*(1.0-q2)
    if denominator == 0.0:
        s_sigmaD = "NA"
        s_sigma1D = "NA"
    else:
        s_sigmaD = s_D/math.sqrt(denominator)
        s_sigma1D = s_D/denominator
    return(s_sigmaD, s_sigma1D)

def calculate_distance(POSN1, POSN2):
    return(abs(round(float(POSN1)*(chr_len-1)) - round(float(POSN2)*(chr_len-1))))

def get_pairwise_stats(l_POS, d_MS):
    l_D_derived, l_D_minor, l_DIST, l_sigmaD_derived, l_sigmaD_minor, l_sigma1D_derived, l_sigma1D_minor = [], [], [], [], [], [], []
    l_AF1, l_AF2 = [], []
    num_of_SNPs = len(l_POS)
    SNP1 = 0
    while SNP1 < num_of_SNPs:
        SNP2 = 0
        while SNP2 < num_of_SNPs:
            if SNP2 > SNP1: #to make sure that the same pair is not used twice
                l_D_derived.append(calculate_D(d_MS[str(SNP1)], d_MS[str(SNP2)]))
                l_D_minor.append(calculate_D(convert_derived_to_minor(d_MS[str(SNP1)]), convert_derived_to_minor(d_MS[str(SNP2)])))
                l_DIST.append(calculate_distance(l_POS[SNP1], l_POS[SNP2]))
                sigma_stats_derived = calculate_correlation_coefficient(d_MS[str(SNP1)], d_MS[str(SNP2)])
                l_sigmaD_derived.append(sigma_stats_derived[0])
                l_sigma1D_derived.append(sigma_stats_derived[1])
                sigma_stats_minor = calculate_correlation_coefficient(convert_derived_to_minor(d_MS[str(SNP1)]), convert_derived_to_minor(d_MS[str(SNP2)]))
                l_sigmaD_minor.append(sigma_stats_minor[0])
                l_sigma1D_minor.append(sigma_stats_minor[1])
                l_AF1.append(get_allele_freq(d_MS[str(SNP1)]))
                l_AF2.append(get_allele_freq(d_MS[str(SNP2)]))
            SNP2 += 1
        SNP1 += 1
    return(l_DIST, l_D_derived, l_D_minor, l_sigmaD_derived, l_sigmaD_minor, l_sigma1D_derived, l_sigma1D_minor, l_AF1, l_AF2)

result =  open(outfolder + "/" + prefix + ".D", 'w+')
result.write("filename" + '\t' + "distance" + '\t' + "D_derived" + '\t' + "D_minor" + '\t' + "sigmaD_derived" + '\t' + "sigmaD_minor" + '\t' + "sigma1D_derived" + '\t' + "sigma1D_minor" + '\t' + "af1" + '\t' + "af2" + '\n')

#go through all simulation replicates and read data into pylibseq format
#addin the option of ignoring some files if they don't exist
os.system("ls " + infolder + "/*" + file_ext + " > " + outfolder + "/" + prefix + ".list")


f_list = open(outfolder + "/" + prefix + ".list", 'r')
numsim = 1
s_absent = 0
for Aline in f_list:
    Aline1 = Aline.strip('\n')
    f_name = Aline1#.split(".")[0]
    f_name_small = Aline1.split("/").pop()
    print ("Reading file:" + Aline1)
    #try:
    if numsim > 0:
        #read ms files and store data in a matrix
        f_ms = open(f_name, 'r')
        t_data = read_ms_file(f_ms)
        f_ms.close()
        l_positions = t_data[0]
        d_ms = t_data[1]

	#calculate pairwise D
        t_LD = get_pairwise_stats(l_positions, d_ms)
        l_distances = t_LD[0]
        l_stat_D_derived = t_LD[1]
        l_stat_D_minor = t_LD[2]
        l_stat_sigmaD_derived = t_LD[3]
        l_stat_sigmaD_minor = t_LD[4]
        l_stat_sigma1D_derived = t_LD[5]
        l_stat_sigma1D_minor = t_LD[6]
        l_af1 = t_LD[7]
        l_af2 = t_LD[8]

        #write down all values:
        i = 0
        for x in l_stat_D_derived:
            result.write(f_name_small + '\t' + str(l_distances[i]) + '\t' + str(x) + '\t' + str(l_stat_D_minor[i]) + '\t' + str(l_stat_sigmaD_derived[i]) + '\t' + str(l_stat_sigmaD_minor[i]) + '\t' + str(l_stat_sigma1D_derived[i]) + '\t' + str(l_stat_sigma1D_minor[i]) + '\t' + str(l_af1[i]) + '\t' + str(l_af2[i]) + '\n')
            i += 1
    #except:    
    else:
        s_absent = s_absent + 1
        print ("This file does not exist or cannot be read or is empty")
	
    numsim = numsim + 1

result.close()
print ("Number of files not read:" + '\t' + str(s_absent))
print ("Finished")






