#This is to summarize the stats from each output file into one table:

import sys

num_sites="thousand_site"
input_folder="/work/users/p/j/pjohri/FitnessNote/simulations/" + num_sites #/additive_site/gamma20/h_0_2"
output_folder="/work/users/p/j/pjohri/FitnessNote/statistics/" + num_sites
fitness_model="additive_site" #multiplicative/additive_site/additive_gene
theta="theta0_01"
epistatic="noEpistasis"
l_gamma=["gamma1000_mean"] #gamma2/gamma20/gamma2_mean/gamma20_mean/gamma100_mean/gamma1000_mean
l_dom=["h_0_0", "h_0_2", "h_0_5"] #h_0_0/h_0_2/h_0_5

for gamma in l_gamma:
    for dom in l_dom:
        print (gamma + '\t' + dom)
        result = open(output_folder + "/" + theta + "/" + epistatic + "/" + fitness_model + "_" + gamma + "_" + dom + ".txt", 'w+')
        result.write("fitness_var" + '\t' + "genetic_load" + '\t' + "inbreeding_load" + '\t' + "af_site1" + '\t' + "af_site2" + '\t' + "af_site_center" + '\n')
        
        repnum=1
        while repnum <= 1000:
            s_fitness_var, s_genetic_load, s_inbreeding_load, s_af1, s_af2, s_af5 = "NA", "NA", "NA", "NA", "NA", "NA"
            f_out = open(input_folder + "/" + fitness_model + "/" + theta + "/" + epistatic + "/" + gamma + "/" + dom + "/output" + str(repnum) + ".txt", 'r')
            for line in f_out:
                line1 = line.strip('\n')
                line2 = line1.split(": ")
                if "variance in fitness" in line:
                    s_fitness_var = line2[1]
                elif "genetic load" in line:
                    s_genetic_load = line2[1]
                elif "inbreeding load" in line:
                    s_inbreeding_load = line2[1]
                elif "allele frequency at site1" in line:
                    s_af1 = line2[1]
                elif "allele frequency at site2" in line:
                    s_af2 = line2[1]
                elif "allele frequency at site" in line:
                    s_af_center = line2[1]
            f_out.close()
            result.write(s_fitness_var + '\t' + s_genetic_load + '\t' + s_inbreeding_load + '\t' + s_af1 + '\t' + s_af2 + '\t' + s_af_center + '\n')
            repnum += 1
        
        result.close()
print("done")



