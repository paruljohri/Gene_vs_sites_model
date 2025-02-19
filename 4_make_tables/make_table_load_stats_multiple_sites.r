#To plot the distribution of loads etc for hte fitness note:

s_folder <- "/work/users/p/j/pjohri/FitnessNote"
num_sites <- "thousand_site"
theta <- "theta0_005"
v_gammas <- c("gamma2", "gamma20", "gamma2_mean", "gamma20_mean", "gamma100_mean", "gamma1000_mean")
v_dom <- c("h_0_0", "h_0_2", "h_0_5")
v_data <- c()
for (s_gamma in v_gammas) {
		print (s_gamma)
		for (s_dom in v_dom){
			t_add_gene <- read.table(paste(s_folder, "/statistics/", num_sites, "/", theta, "/noEpistasis/additive_gene_", s_gamma, "_", s_dom, ".txt", sep=""), h=T)
			t_add_site <- read.table(paste(s_folder, "/statistics/", num_sites, "/", theta, "/noEpistasis/additive_site_", s_gamma, "_", s_dom, ".txt", sep=""), h=T)
			t_af_gene <- read.table(paste(s_folder, "/statistics/", num_sites, "/", theta, "/noEpistasis/additive_gene_", s_gamma, "_", s_dom, ".af", sep=""), h=T)
			t_af_site <- read.table(paste(s_folder, "/statistics/", num_sites, "/", theta, "/noEpistasis/additive_site_", s_gamma, "_", s_dom, ".af", sep=""), h=T)
			#genetic load:
			v_data <- rbind(v_data, c(s_gamma, s_dom, "genetic_load", "mean", mean(t_add_gene$genetic_load), mean(t_add_site$genetic_load)))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "genetic_load", "SE", sd(t_add_gene$genetic_load)/sqrt(length(t_add_gene$genetic_load)), sd(t_add_site$genetic_load)/sqrt(length(t_add_site$genetic_load))))
			#inbreeding load:
			v_data <- rbind(v_data, c(s_gamma, s_dom, "inbreeding_load", "mean", mean(t_add_gene$inbreeding_load), mean(t_add_site$inbreeding_load)))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "inbreeding_load", "SE", sd(t_add_gene$inbreeding_load)/sqrt(length(t_add_gene$inbreeding_load)), sd(t_add_site$inbreeding_load)/sqrt(length(t_add_site$inbreeding_load))))
			#varience in fitness:
			v_data <- rbind(v_data, c(s_gamma, s_dom, "variance_in_fitness", "mean", mean(t_add_gene$fitness_var), mean(t_add_site$fitness_var)))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "variance_in_fitness", "SE", sd(t_add_gene$fitness_var)/sqrt(length(t_add_gene$fitness_var)), sd(t_add_site$fitness_var)/sqrt(length(t_add_site$fitness_var))))
			#allele frequency at a site in the center:
			v_data <- rbind(v_data, c(s_gamma, s_dom, "allele_frequency_at_center", "mean", mean(t_add_gene$af_site_center), mean(t_add_site$af_site_center)))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "allele_frequency_at_center", "SE", sd(t_add_gene$af_site_center)/sqrt(length(t_add_gene$af_site_center)), sd(t_add_site$af_site_center)/sqrt(length(t_add_site$af_site_center))))
			#mean allele frequency across all sites for a sample of 100 genomes:
			v_data <- rbind(v_data, c(s_gamma, s_dom, "allele_frequency_from_sample", "mean", mean(t_af_gene$mean_allele_freq), mean(t_af_site$mean_allele_freq)))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "allele_frequency_from_sample", "SE", sd(t_af_gene$mean_allele_freq)/sqrt(length(t_af_gene$mean_allele_freq)), sd(t_af_site$mean_allele_freq)/sqrt(length(t_af_site$mean_allele_freq))))
		}
}

#name columns:
colnames(v_data) <- c("gamma", "dominance", "summary", "type", "additive_gene", "additive_site")
write.table(v_data, file = paste(s_folder, "/Tables/", num_sites, "/", theta, "/noEpistasis/summary_load.txt", sep=""), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)







