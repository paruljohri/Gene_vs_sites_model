#To plot the distribution of loads etc for the fitness note:
#With epistasis

num_sites <- "thousand_site"
theta <- "theta0_005"
v_gammas <- c("gamma100_mean_lowrec") #keep this a single value
v_epsilon <- c("epsilon_0_00", "epsilon_0_02", "epsilon_0_04", "epsilon_0_08", "epsilon_0_16", "epsilon_0_32")
v_dom <- c("h_0_0", "h_0_2", "h_0_5")
s_folder <- "/work/users/p/j/pjohri/FitnessNote"
v_data <- c()
for (s_gamma in v_gammas) {
		print (s_gamma)
		for (s_epsilon in v_epsilon){
		print (s_epsilon)
			for (s_dom in v_dom){
				t_add_gene <- read.table(paste(s_folder, "/statistics/", num_sites, "/", theta, "/Epistasis/additive_gene_", s_gamma, "_", s_epsilon, "_", s_dom, ".txt", sep=""), h=T)
				t_add_site <- read.table(paste(s_folder, "/statistics/", num_sites, "/", theta, "/Epistasis/additive_site_", s_gamma, "_", s_epsilon, "_", s_dom, ".txt", sep=""), h=T)
				#genetic load:
				v_data <- rbind(v_data, c(s_gamma, s_epsilon, s_dom, "genetic_load", "mean", mean(t_add_gene$genetic_load), mean(t_add_site$genetic_load)))
				v_data <- rbind(v_data, c(s_gamma, s_epsilon, s_dom, "genetic_load", "SE", sd(t_add_gene$genetic_load)/sqrt(length(t_add_gene$genetic_load)), sd(t_add_site$genetic_load)/sqrt(length(t_add_site$genetic_load))))
				#inbreeding load:
				v_data <- rbind(v_data, c(s_gamma, s_epsilon, s_dom, "inbreeding_load", "mean", mean(t_add_gene$inbreeding_load), mean(t_add_site$inbreeding_load)))
				v_data <- rbind(v_data, c(s_gamma, s_epsilon, s_dom, "inbreeding_load", "SE", sd(t_add_gene$inbreeding_load)/sqrt(length(t_add_gene$inbreeding_load)), sd(t_add_site$inbreeding_load)/sqrt(length(t_add_site$inbreeding_load))))
				#inbreeding load sample:
				v_data <- rbind(v_data, c(s_gamma, s_epsilon, s_dom, "inbreeding_load_sample", "mean", mean(t_add_gene$inbreeding_load_sample), mean(t_add_site$inbreeding_load_sample)))
				v_data <- rbind(v_data, c(s_gamma, s_epsilon, s_dom, "inbreeding_load_sample", "SE", sd(t_add_gene$inbreeding_load_sample)/sqrt(length(t_add_gene$inbreeding_load_sample)), sd(t_add_site$inbreeding_load_sample)/sqrt(length(t_add_site$inbreeding_load_sample))))
				#variance in fitness:
				v_data <- rbind(v_data, c(s_gamma, s_epsilon, s_dom, "variance_in_fitness", "mean", mean(t_add_gene$fitness_var), mean(t_add_site$fitness_var)))
				v_data <- rbind(v_data, c(s_gamma, s_epsilon, s_dom, "variance_in_fitness", "SE", sd(t_add_gene$fitness_var)/sqrt(length(t_add_gene$fitness_var)), sd(t_add_site$fitness_var)/sqrt(length(t_add_site$fitness_var))))
				#allele frequency at a site in the center:
				v_data <- rbind(v_data, c(s_gamma, s_epsilon, s_dom, "allele_frequency", "mean", mean(t_add_gene$af_site_center), mean(t_add_site$af_site_center)))
				v_data <- rbind(v_data, c(s_gamma, s_epsilon, s_dom, "allele_frequency", "SE", sd(t_add_gene$af_site_center)/sqrt(length(t_add_gene$af_site_center)), sd(t_add_site$af_site_center)/sqrt(length(t_add_site$af_site_center))))
			}
		}
}

#name columns:
colnames(v_data) <- c("gamma", "epsilon", "dominance", "summary", "type", "additive_gene", "additive_site")
write.table(v_data, file = paste(s_folder, "/Tables/", num_sites, "/", theta, "/Epistasis/summary_load_", s_gamma, ".txt", sep=""), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)







