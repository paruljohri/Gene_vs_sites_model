#This is to get the mean allele frequency and the number of segregating sites from a .ms file:
#Note that this will include all sites, not just SNPs
#How to run:
##python get_af_from_ms_general_reps.py -input_folder -output_folder -extension .ms -output_prefix -num_indv 100 -num_sites 1000
##python get_af_from_ms_general_reps.py -input_folder -output_folder /work/users/p/j/pjohri/FitnessNote/statistics/thousand_site/theta0_005/noEpistasis -extension .ms -output_prefix additive_site_gamma2_mean_h_0_5 -num_indv 100 -num_sites 1000

import sys
import argparse
import os

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-input_folder', dest = 'input_folder', action='store', nargs = 1, type = str, help = 'full path to folder with .ms files')
parser.add_argument('-output_folder', dest = 'output_folder', action='store', nargs = 1, type = str, help = 'full path to folder where you want to write the output')
parser.add_argument('-extension', dest = 'extension', action='store', nargs = 1, type = str, help = '*extension files will be read')
parser.add_argument('-output_prefix', dest = 'output_prefix', action='store', nargs = 1, type = str, help = 'name of output file without extension')
parser.add_argument('-num_indv', dest = 'num_indv', action='store', nargs = 1, type = int, help = 'number  of individuals')
parser.add_argument('-num_sites', dest = 'num_sites', action='store', nargs = 1, type = int, help = 'number  of sites in total')

#read input parameters
args = parser.parse_args()
in_folder = args.input_folder[0]
out_folder = args.output_folder[0]
s_ext = args.extension[0]
prefix = args.output_prefix[0]
num_indv = args.num_indv[0]
num_sites = args.num_sites[0]
print (out_folder)

#reading in genotypes from ms
def get_af_from_ms(f_MS):
    d_AF = {}
    linecount = 0
    for line in f_MS:
        line1 = line.strip('\n')
        if "//" not in line1 and "segsites" not in line1 and "positions" not in line1:
            linecount += 1
            if linecount <= int(num_indv):
                col = 1
                for x in line1:
                    try:
                        d_AF[col] = int(d_AF[col]) + int(x)
                    except:
                        d_AF[col] = int(x)
                    col += 1
    return(d_AF)

#Open output file
result = open(out_folder + "/" + prefix + ".af", 'w+')
result.write("filename" + '\t' + "mean_allele_freq" + '\t' + "S" + '\n')

#Make a list of all .ms files:
os.system("ls " + in_folder + "/*" + s_ext + " > " + out_folder + "/" + prefix + ".list")

f_list = open(out_folder + "/" + prefix + ".list", 'r')
for Aline in f_list:
    Aline1 = Aline.strip('\n')
    f_name = Aline1.split("/").pop()
    print ("Reading file:" + Aline1)

    #reading in genotypes from ms
    f_ms = open(in_folder + "/" + f_name, 'r')
    d_af = get_af_from_ms(f_ms)
    f_ms.close()
    
    #get mean allele frequency and number of polymorphic sites:
    sum_af = 0.0
    s_seg_sites = 0
    for posn in d_af.keys():
        s_q = float(d_af[posn])/float(num_indv)
        sum_af = sum_af + s_q
        if s_q > 0.0 and s_q < 1.0:
            s_seg_sites += 1
    
    #Write the full result:
    result.write(f_name + '\t' + str(sum_af/float(num_sites)) + '\t' + str(s_seg_sites) + '\n')

f_list.close()
result.close()
print ("done")

