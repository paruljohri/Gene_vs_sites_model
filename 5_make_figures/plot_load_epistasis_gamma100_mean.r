#to make a plot for epistasis:

setwd("OneDrive - University of North Carolina at Chapel Hill/Work/Projects/Fitness_note_Brian/")

options(scipen=999)

t <- read.table("Tables/thousand_site/theta0_01/Epistasis/summary_load.txt", h=T)

#A function to add arrows on the chart
error.bar <- function(x, y, upper, lower=upper, length=0.05){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length)
}

epsilon_values <- c("0.00", "0.02", "0.04", "0.08", "0.16", "0.32")
v_colors <- c("turquoise4" , "cyan2")

#One figure with genetic load and inbreeding load
par(mfcol=c(3,2))
par(mar=c(5,5,2,1))

#using basic R:
#genetic load; h=0.0
t_sub_mean <- t[which(t$summary=="genetic_load" & t$dominance=="h_0_0" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="genetic_load" & t$dominance=="h_0_0" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , ylab=expression(italic("L")), names.arg=epsilon_values, xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(0,0.012))
error.bar(my_barplot, data_mean, data_se)
title(main = expression(paste(italic("h"), "=0.0")), cex.main=2, adj=0)

#genetic load; h=0.2
t_sub_mean <- t[which(t$summary=="genetic_load" & t$dominance=="h_0_2" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="genetic_load" & t$dominance=="h_0_2" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , ylab=expression(italic("L")), names.arg=epsilon_values, xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(0,0.012))
error.bar(my_barplot, data_mean, data_se)
title(main = expression(paste(italic("h"), "=0.2")), cex.main=2, adj=0)

#genetic load; h=0.5
t_sub_mean <- t[which(t$summary=="genetic_load" & t$dominance=="h_0_5" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="genetic_load" & t$dominance=="h_0_5" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , ylab=expression(italic("L")), names.arg=epsilon_values, xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(0,0.012))
error.bar(my_barplot, data_mean, data_se)
title(main = expression(paste(italic("h"), "=0.5")), cex.main=2, adj=0)

#inbreeding load; h=0.0
t_sub_mean <- t[which(t$summary=="inbreeding_load" & t$dominance=="h_0_0" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="inbreeding_load" & t$dominance=="h_0_0" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , ylab=expression(italic("B")), names.arg=epsilon_values, xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(0,0.07))
error.bar(my_barplot, data_mean, data_se)
legend("topright", legend=c("gene", "sites"), col=v_colors, fill=v_colors, cex=1.5, box.lty=0)
title(main = expression(paste(italic("h"), "=0.0")), cex.main=2, adj=0)

#inbreeding load; h=0.2
t_sub_mean <- t[which(t$summary=="inbreeding_load" & t$dominance=="h_0_2" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="inbreeding_load" & t$dominance=="h_0_2" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , ylab=expression(italic("B")), names.arg=epsilon_values, xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(0,0.015))
error.bar(my_barplot, data_mean, data_se)
title(main = expression(paste(italic("h"), "=0.2")), cex.main=2, adj=0)

#inbreeding load; h=0.5
t_sub_mean <- t[which(t$summary=="inbreeding_load" & t$dominance=="h_0_5" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="inbreeding_load" & t$dominance=="h_0_5" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , ylab=expression(italic("B")), names.arg=epsilon_values, xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(min(t_sub_mean$additive_gene,t_sub_mean$additive_site),max(t_sub_mean$additive_gene,t_sub_mean$additive_site)))
error.bar(my_barplot, data_mean, data_se)
title(main = expression(paste(italic("h"), "=0.5")), cex.main=2, adj=0)


#One figure with allele frequency and variance in fitness
par(mfcol=c(3,2))
par(mar=c(5,5,2,1))

#using basic R:
#variance in fitness:
y_min <- min(c(t$additive_gene[which(t$summary=="variance_in_fitness" & t$type=="mean")]-t$additive_gene[which(t$summary=="variance_in_fitness" & t$type=="SE")], t$additive_site[which(t$summary=="variance_in_fitness" & t$type=="mean")]-t$additive_site[which(t$summary=="variance_in_fitness" & t$type=="SE")]))
y_max <- max(c(t$additive_gene[which(t$summary=="variance_in_fitness" & t$type=="mean")]+t$additive_gene[which(t$summary=="variance_in_fitness" & t$type=="SE")], t$additive_site[which(t$summary=="variance_in_fitness" & t$type=="mean")]+t$additive_site[which(t$summary=="variance_in_fitness" & t$type=="SE")]))

#variance in fitness; h=0.0
t_sub_mean <- t[which(t$summary=="variance_in_fitness" & t$dominance=="h_0_0" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="variance_in_fitness" & t$dominance=="h_0_0" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , ylab=expression(italic("Var")), names.arg=epsilon_values, xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(0, y_max))
error.bar(my_barplot, data_mean, data_se)
title(main = expression(paste(italic("h"), "=0.0")), cex.main=2, adj=0)

#variance in fitness; h=0.2
t_sub_mean <- t[which(t$summary=="variance_in_fitness" & t$dominance=="h_0_2" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="variance_in_fitness" & t$dominance=="h_0_2" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , ylab=expression(italic("Var")), names.arg=epsilon_values, xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(0, y_max))
error.bar(my_barplot, data_mean, data_se)
title(main = expression(paste(italic("h"), "=0.2")), cex.main=2, adj=0)

#variance in fitness; h=0.5
t_sub_mean <- t[which(t$summary=="variance_in_fitness" & t$dominance=="h_0_5" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="variance_in_fitness" & t$dominance=="h_0_5" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , ylab=expression(italic("Var")), names.arg=epsilon_values, xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(0, y_max))
error.bar(my_barplot, data_mean, data_se)
title(main = expression(paste(italic("h"), "=0.5")), cex.main=2, adj=0)

#allele frequency:
y_min <- min(c(t$additive_gene[which(t$summary=="allele_frequency" & t$type=="mean")]-t$additive_gene[which(t$summary=="allele_frequency" & t$type=="SE")], t$additive_site[which(t$summary=="allele_frequency" & t$type=="mean")]-t$additive_site[which(t$summary=="allele_frequency" & t$type=="SE")]))
y_max <- max(c(t$additive_gene[which(t$summary=="allele_frequency" & t$type=="mean")]+t$additive_gene[which(t$summary=="allele_frequency" & t$type=="SE")], t$additive_site[which(t$summary=="allele_frequency" & t$type=="mean")]+t$additive_site[which(t$summary=="allele_frequency" & t$type=="SE")]))

#allele frequency; h=0.0
t_sub_mean <- t[which(t$summary=="allele_frequency" & t$dominance=="h_0_0" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="allele_frequency" & t$dominance=="h_0_0" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , ylab=expression(italic("q")), names.arg=epsilon_values, xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(0, y_max))
error.bar(my_barplot, data_mean, data_se)
title(main = expression(paste(italic("h"), "=0.0")), cex.main=2, adj=0)
legend("topright", legend=c("gene", "sites"), col=v_colors, fill=v_colors, cex=1.5, box.lty=0, horiz=T)

#allele frequency; h=0.2
t_sub_mean <- t[which(t$summary=="allele_frequency" & t$dominance=="h_0_2" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="allele_frequency" & t$dominance=="h_0_2" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , ylab=expression(italic("q")), names.arg=epsilon_values, xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(0, y_max))
error.bar(my_barplot, data_mean, data_se)
title(main = expression(paste(italic("h"), "=0.2")), cex.main=2, adj=0)

#allele frequency; h=0.5
t_sub_mean <- t[which(t$summary=="allele_frequency" & t$dominance=="h_0_5" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="allele_frequency" & t$dominance=="h_0_5" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , ylab=expression(italic("q")), names.arg=epsilon_values, xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(0, y_max))
error.bar(my_barplot, data_mean, data_se)
title(main = expression(paste(italic("h"), "=0.5")), cex.main=2, adj=0)
