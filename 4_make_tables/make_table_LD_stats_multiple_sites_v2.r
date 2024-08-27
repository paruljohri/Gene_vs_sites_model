#To plot the distribution of LD stats for the fitness note:
#Here we are trying to get SEs for all stats across replicates
#setwd("OneDrive - University of North Carolina at Chapel Hill/Work/Projects/Fitness_note_Brian/")

num_sites <- "thousand_site"
theta <- "theta0_005"
v_gammas <- c("gamma100") #c("gamma2", "gamma20", "gamma100", "gamma2_mean", "gamma20_mean", "gamma100_mean", "gamma1000_mean")
v_dom <- c("h_0_0", "h_0_2", "h_0_5")
dist_min <- 1 #1 or 800 or 500
dist_max <- 100 #25 or 100 or 1000 or 600
s_folder <- "/work/users/p/j/pjohri/FitnessNote"
v_data <- c()
for (s_gamma in v_gammas) {
		print (s_gamma)
		for (s_dom in v_dom){
			print (s_dom)
			t_add_gene <- read.table(paste(s_folder, "/statistics/", num_sites, "/", theta, "/noEpistasis/additive_gene_", s_gamma, "_", s_dom, ".D", sep=""), h=T)
			t_add_site <- read.table(paste(s_folder, "/statistics/", num_sites, "/", theta, "/noEpistasis/additive_site_", s_gamma, "_", s_dom, ".D", sep=""), h=T)
			
			t_add_gene$product <- t_add_gene$af1*(1.0-t_add_gene$af1)*t_add_gene$af2*(1.0-t_add_gene$af2)
			t_add_site$product <- t_add_site$af1*(1.0-t_add_site$af1)*t_add_site$af2*(1.0-t_add_site$af2)
			
			#derived allele:
			
			#filtering by distance and for non-monomorphic sites
			v_D_gene <- t_add_gene$D_derived[which(t_add_gene$sigmaD_derived!="NA" & t_add_gene$distance > dist_min & t_add_gene$distance <= dist_max)]
			v_product_gene <- t_add_gene$product[which(t_add_gene$sigmaD_derived!="NA" & t_add_gene$distance > dist_min & t_add_gene$distance <= dist_max)]
			v_filename_gene <- t_add_gene$filename[which(t_add_gene$sigmaD_derived!="NA" & t_add_gene$distance > dist_min & t_add_gene$distance <= dist_max)]
			v_D_site <- t_add_site$D_derived[which(t_add_site$sigmaD_derived!="NA" & t_add_site$distance > dist_min & t_add_site$distance <= dist_max)]
			v_product_site <- t_add_site$product[which(t_add_site$sigmaD_derived!="NA" & t_add_site$distance > dist_min & t_add_site$distance <= dist_max)]
			v_filename_site <- t_add_site$filename[which(t_add_site$sigmaD_derived!="NA" & t_add_site$distance > dist_min & t_add_site$distance <= dist_max)]
			
			#getting means across replicates
			t_D_reps_gene <- aggregate(v_D_gene~v_filename_gene, FUN=mean)
			t_D_reps_site <- aggregate(v_D_site~v_filename_site, FUN=mean)
			t_product_reps_gene <- aggregate(v_product_gene~v_filename_gene, FUN=mean)
			t_product_reps_site <- aggregate(v_product_site~v_filename_site, FUN=mean)
			
			v_D_mean_gene <- t_D_reps_gene$v_D_gene
			v_D_mean_site <- t_D_reps_site$v_D_site
			v_product_mean_gene <- t_product_reps_gene$v_product_gene
			v_product_mean_site <- t_product_reps_site$v_product_site
			v_sigmaD_mean_gene <- t_D_reps_gene$v_D_gene/sqrt(v_product_mean_gene)
			v_sigmaD_mean_site <- t_D_reps_site$v_D_site/sqrt(v_product_mean_site)
			
			v_data <- rbind(v_data, c(s_gamma, s_dom, "D_derived", "mean", mean(v_D_mean_gene), mean(v_D_mean_site)))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "D_derived", "SE", sd(v_D_mean_gene)/sqrt(length(v_D_mean_gene)), sd(v_D_mean_site)/sqrt(length(v_D_mean_site))))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "sigmaD_derived", "mean", mean(v_sigmaD_mean_gene), mean(v_sigmaD_mean_site)))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "sigmaD_derived", "SE", sd(v_sigmaD_mean_gene)/sqrt(length(v_sigmaD_mean_gene)), sd(v_sigmaD_mean_site)/sqrt(length(v_sigmaD_mean_site))))
			
			#minor allele:
			#filtering by distance and for non-monomorphic sites
			v_D_gene <- t_add_gene$D_minor[which(t_add_gene$sigmaD_minor!="NA" & t_add_gene$distance > dist_min & t_add_gene$distance <= dist_max)]
			v_product_gene <- t_add_gene$product[which(t_add_gene$sigmaD_minor!="NA" & t_add_gene$distance > dist_min & t_add_gene$distance <= dist_max)]
			v_filename_gene <- t_add_gene$filename[which(t_add_gene$sigmaD_minor!="NA" & t_add_gene$distance > dist_min & t_add_gene$distance <= dist_max)]
			v_D_site <- t_add_site$D_minor[which(t_add_site$sigmaD_minor!="NA" & t_add_site$distance > dist_min & t_add_site$distance <= dist_max)]
			v_product_site <- t_add_site$product[which(t_add_site$sigmaD_minor!="NA" & t_add_site$distance > dist_min & t_add_site$distance <= dist_max)]
			v_filename_site <- t_add_site$filename[which(t_add_site$sigmaD_minor!="NA" & t_add_site$distance > dist_min & t_add_site$distance <= dist_max)]
			
			#getting means across replicates
			t_D_reps_gene <- aggregate(v_D_gene~v_filename_gene, FUN=mean)
			t_D_reps_site <- aggregate(v_D_site~v_filename_site, FUN=mean)
			t_product_reps_gene <- aggregate(v_product_gene~v_filename_gene, FUN=mean)
			t_product_reps_site <- aggregate(v_product_site~v_filename_site, FUN=mean)
			
			v_D_mean_gene <- t_D_reps_gene$v_D_gene
			v_D_mean_site <- t_D_reps_site$v_D_site
			v_product_mean_gene <- t_product_reps_gene$v_product_gene
			v_product_mean_site <- t_product_reps_site$v_product_site
			v_sigmaD_mean_gene <- t_D_reps_gene$v_D_gene/sqrt(v_product_mean_gene)
			v_sigmaD_mean_site <- t_D_reps_site$v_D_site/sqrt(v_product_mean_site)
			
			v_data <- rbind(v_data, c(s_gamma, s_dom, "D_minor", "mean", mean(v_D_mean_gene), mean(v_D_mean_site)))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "D_minor", "SE", sd(v_D_mean_gene)/sqrt(length(v_D_mean_gene)), sd(v_D_mean_site)/sqrt(length(v_D_mean_site))))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "sigmaD_minor", "mean", mean(v_sigmaD_mean_gene), mean(v_sigmaD_mean_site)))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "sigmaD_minor", "SE", sd(v_sigmaD_mean_gene)/sqrt(length(v_sigmaD_mean_gene)), sd(v_sigmaD_mean_site)/sqrt(length(v_sigmaD_mean_site))))
			
		}
}

#name columns:
colnames(v_data) <- c("gamma", "dominance", "summary", "type", "additive_gene", "additive_site")
write.table(v_data, file = paste(s_folder, "/Tables/", num_sites, "/", theta, "/noEpistasis/summary_LD_v2_", dist_min, "_", dist_max, ".txt", sep=""), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)







