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
			t_add_gene <- read.table(paste(s_folder, "/statistics/", num_sites, "/", theta, "/additive_gene_", s_gamma, "_", s_dom, ".txt", sep=""), h=T)
			t_add_site <- read.table(paste(s_folder, "/statistics/", num_sites, "/", theta, "/additive_site_", s_gamma, "_", s_dom, ".txt", sep=""), h=T)
			#genetic load:
			v_data <- rbind(v_data, c(s_gamma, s_dom, "genetic load", "mean", mean(t_add_gene$genetic_load), mean(t_add_site$genetic_load)))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "genetic load", "SE", sd(t_add_gene$genetic_load)/sqrt(length(t_add_gene$genetic_load)), sd(t_add_site$genetic_load)/sqrt(length(t_add_site$genetic_load))))
			#inbreeding load:
			v_data <- rbind(v_data, c(s_gamma, s_dom, "inbreeding load", "mean", mean(t_add_gene$inbreeding_load), mean(t_add_site$inbreeding_load)))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "inbreeding load", "SE", sd(t_add_gene$inbreeding_load)/sqrt(length(t_add_gene$inbreeding_load)), sd(t_add_site$inbreeding_load)/sqrt(length(t_add_site$inbreeding_load))))
			#varience in fitness:
			v_data <- rbind(v_data, c(s_gamma, s_dom, "variance in fitness", "mean", mean(t_add_gene$fitness_var), mean(t_add_site$fitness_var)))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "variance in fitness", "SE", sd(t_add_gene$fitness_var)/sqrt(length(t_add_gene$fitness_var)), sd(t_add_site$fitness_var)/sqrt(length(t_add_site$fitness_var))))
			#allele frequency at a site in the center:
			v_data <- rbind(v_data, c(s_gamma, s_dom, "allele frequency", "mean", mean(t_add_gene$af_site_center), mean(t_add_site$af_site_center)))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "allele frequency", "SE", sd(t_add_gene$af_site_center)/sqrt(length(t_add_gene$af_site_center)), sd(t_add_site$af_site_center)/sqrt(length(t_add_site$af_site_center))))
		}
}

#name columns:
colnames(v_data) <- c("gamma", "dominance", "summary", "type", "additive_gene", "additive_site")
write.table(v_data, file = paste(s_folder, "/Tables/", num_sites, "/", theta, "/summary_load.txt", sep=""), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)







