#To plot the distribution of LD stats for the fitness note:
#Here we are trying to get SEs for all stats across replicates
#setwd("OneDrive - University of North Carolina at Chapel Hill/Work/Projects/Fitness_note_Brian/")

num_sites <- "thousand_site"
theta <- "theta0_005"
v_gammas <- c("gamma2_mean", "gamma20_mean", "gamma100_mean", "gamma1000_mean")
v_dom <- c("h_0_0", "h_0_2", "h_0_5")
dist_min <- 1 #1 or 800 or 500
dist_max <- 100 #25 or 100 or 1000 or 600
s_folder <- "/work/users/p/j/pjohri/FitnessNote"
v_data <- c()
for (s_gamma in v_gammas) {
		print (s_gamma)
		for (s_dom in v_dom){
			print (s_dom)
			t_multi <- read.table(paste(s_folder, "/statistics/", num_sites, "/", theta, "/noEpistasis/multiplicative_", s_gamma, "_", s_dom, ".D", sep=""), h=T)
			
			t_multi$product <- t_multi$af1*(1.0-t_multi$af1)*t_multi$af2*(1.0-t_multi$af2)
			
			#derived allele:
			
			#filtering by distance and for non-monomorphic sites
			v_D_multi <- t_multi$D_derived[which(t_multi$sigmaD_derived!="NA" & t_multi$distance > dist_min & t_multi$distance <= dist_max)]
			v_product_multi <- t_multi$product[which(t_multi$sigmaD_derived!="NA" & t_multi$distance > dist_min & t_multi$distance <= dist_max)]
			v_filename_multi <- t_multi$filename[which(t_multi$sigmaD_derived!="NA" & t_multi$distance > dist_min & t_multi$distance <= dist_max)]
			
			#getting means across replicates
			t_D_reps_multi <- aggregate(v_D_multi~v_filename_multi, FUN=mean)
			t_product_reps_multi <- aggregate(v_product_multi~v_filename_multi, FUN=mean)
			
			v_D_mean_multi <- t_D_reps_multi$v_D_multi
			v_product_mean_multi <- t_product_reps_multi$v_product_multi
			v_sigmaD_mean_multi <- t_D_reps_multi$v_D_multi/sqrt(v_product_mean_multi)
			
			v_data <- rbind(v_data, c(s_gamma, s_dom, "D_derived", "mean", mean(v_D_mean_multi)))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "D_derived", "SE", sd(v_D_mean_multi)/sqrt(length(v_D_mean_multi))))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "sigmaD_derived", "mean", mean(v_sigmaD_mean_multi)))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "sigmaD_derived", "SE", sd(v_sigmaD_mean_multi)/sqrt(length(v_sigmaD_mean_multi))))
			
			#minor allele:
			#filtering by distance and for non-monomorphic sites
			v_D_multi <- t_multi$D_minor[which(t_multi$sigmaD_minor!="NA" & t_multi$distance > dist_min & t_multi$distance <= dist_max)]
			v_product_multi <- t_multi$product[which(t_multi$sigmaD_minor!="NA" & t_multi$distance > dist_min & t_multi$distance <= dist_max)]
			v_filename_multi <- t_multi$filename[which(t_multi$sigmaD_minor!="NA" & t_multi$distance > dist_min & t_multi$distance <= dist_max)]
			
			#getting means across replicates
			t_D_reps_multi <- aggregate(v_D_multi~v_filename_multi, FUN=mean)
			t_product_reps_multi <- aggregate(v_product_multi~v_filename_multi, FUN=mean)
			
			v_D_mean_multi <- t_D_reps_multi$v_D_multi
			v_product_mean_multi <- t_product_reps_multi$v_product_multi
			v_sigmaD_mean_multi <- t_D_reps_multi$v_D_multi/sqrt(v_product_mean_multi)
			
			v_data <- rbind(v_data, c(s_gamma, s_dom, "D_minor", "mean", mean(v_D_mean_multi)))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "D_minor", "SE", sd(v_D_mean_multi)/sqrt(length(v_D_mean_multi))))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "sigmaD_minor", "mean", mean(v_sigmaD_mean_multi)))
			v_data <- rbind(v_data, c(s_gamma, s_dom, "sigmaD_minor", "SE", sd(v_sigmaD_mean_multi)/sqrt(length(v_sigmaD_mean_multi))))
			
		}
}

#name columns:
colnames(v_data) <- c("gamma", "dominance", "summary", "type", "multiplicative")
write.table(v_data, file = paste(s_folder, "/Tables/", num_sites, "/", theta, "/noEpistasis/summary_LD_v2_multi_", dist_min, "_", dist_max, ".txt", sep=""), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)







