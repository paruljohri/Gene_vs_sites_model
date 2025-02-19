#to make a plot for epistasis:

setwd("Work/Projects/Gene_vs_sites_model/")

options(scipen=999)

#A function to add arrows on the chart
error.bar <- function(x, y, upper, lower=upper, length=0.05){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length)
}

t <- read.table("/Users/paruljohri/Work/Projects/Gene_vs_sites_model/Tables/thousand_site/theta0_005/noEpistasis/summary_LD_v2_1_100.txt", h=T)

y_min_derived <- min(c(t$additive_gene[which(t$summary=="sigmaD_derived" & t$type=="mean")]-t$additive_gene[which(t$summary=="sigmaD_derived" & t$type=="SE")], t$additive_site[which(t$summary=="sigmaD_derived" & t$type=="mean")]-t$additive_site[which(t$summary=="sigmaD_derived" & t$type=="SE")]))
y_max_derived <- max(c(t$additive_gene[which(t$summary=="sigmaD_derived" & t$type=="mean")]+t$additive_gene[which(t$summary=="sigmaD_derived" & t$type=="SE")], t$additive_site[which(t$summary=="sigmaD_derived" & t$type=="mean")]+t$additive_site[which(t$summary=="sigmaD_derived" & t$type=="SE")]))

y_min_minor <- min(c(t$additive_gene[which(t$summary=="sigmaD_minor" & t$type=="mean")]-t$additive_gene[which(t$summary=="sigmaD_minor" & t$type=="SE")], t$additive_site[which(t$summary=="sigmaD_minor" & t$type=="mean")]-t$additive_site[which(t$summary=="sigmaD_minor" & t$type=="SE")]))
y_max_minor <- max(c(t$additive_gene[which(t$summary=="sigmaD_minor" & t$type=="mean")]+t$additive_gene[which(t$summary=="sigmaD_minor" & t$type=="SE")], t$additive_site[which(t$summary=="sigmaD_minor" & t$type=="mean")]+t$additive_site[which(t$summary=="sigmaD_minor" & t$type=="SE")]))

#dom_values <- c(expression(paste(italic("h"), "=0.0")), expression(paste(italic("h"), "=0.2")), expression(paste(italic("h"), "=0.5")))
dom_values <- c("0.0", "0.2", "0.5")
v_colors <- c("firebrick2" , "gray39")

#One figure with sigmaD
par(mfrow=c(2,4))
par(mar=c(5,5,2,1))

#using basic R:
#sigmaD_derived; gamma2_mean
t_sub_mean <- t[which(t$summary=="sigmaD_derived" & t$gamma=="gamma2_mean" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="sigmaD_derived" & t$gamma=="gamma2_mean" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , ylab=expression(paste(italic(sigma)["D"], " (selected)")), names.arg=dom_values, cex.lab=1.8, cex.axis=1.5, cex.names=1.7, ylim=c(y_min_derived, y_max_derived))
error.bar(my_barplot, data_mean, data_se)
title(main = expression(paste(bar(italic(gamma)), "=2")), cex.main=2, adj=0)

#sigmaD_derived; gamma20_mean
t_sub_mean <- t[which(t$summary=="sigmaD_derived" & t$gamma=="gamma20_mean" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="sigmaD_derived" & t$gamma=="gamma20_mean" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , names.arg=dom_values, cex.lab=1.8, cex.axis=1.5, cex.names=1.7, ylim=c(y_min_derived, y_max_derived))
error.bar(my_barplot, data_mean, data_se)
title(main = expression(paste(bar(italic(gamma)), "=20")), cex.main=2, adj=0)

#sigmaD_derived; gamma100_mean
t_sub_mean <- t[which(t$summary=="sigmaD_derived" & t$gamma=="gamma100_mean" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="sigmaD_derived" & t$gamma=="gamma100_mean" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , names.arg=dom_values, cex.lab=1.8, cex.axis=1.5, cex.names=1.7, ylim=c(y_min_derived, y_max_derived))
error.bar(my_barplot, data_mean, data_se)
title(main = expression(paste(bar(italic(gamma)), "=100")), cex.main=2, adj=0)

#sigmaD_derived; gamma1000_mean
t_sub_mean <- t[which(t$summary=="sigmaD_derived" & t$gamma=="gamma1000_mean" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="sigmaD_derived" & t$gamma=="gamma1000_mean" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , names.arg=dom_values, cex.lab=1.8, cex.axis=1.5, cex.names=1.7, ylim=c(y_min_derived, y_max_derived))
error.bar(my_barplot, data_mean, data_se)
title(main = expression(paste(bar(italic(gamma)), "=1000")), cex.main=2, adj=0)



#sigmaD_minor; gamma2_mean
t_sub_mean <- t[which(t$summary=="sigmaD_minor" & t$gamma=="gamma2_mean" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="sigmaD_minor" & t$gamma=="gamma2_mean" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , ylab=expression(paste(italic(sigma)["D"], " (minor)")), names.arg=dom_values, xlab=expression(paste(italic("h"))), cex.lab=1.8, cex.axis=1.5, cex.names=1.7, ylim=c(y_min_minor, y_max_minor))
error.bar(my_barplot, data_mean, data_se)


#sigmaD_minor; gamma20_mean
t_sub_mean <- t[which(t$summary=="sigmaD_minor" & t$gamma=="gamma20_mean" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="sigmaD_minor" & t$gamma=="gamma20_mean" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , names.arg=dom_values, xlab=expression(paste(italic("h"))), cex.lab=1.8, cex.axis=1.5, cex.names=1.7, ylim=c(y_min_minor, y_max_minor))
error.bar(my_barplot, data_mean, data_se)


#sigmaD_minor; gamma100_mean
t_sub_mean <- t[which(t$summary=="sigmaD_minor" & t$gamma=="gamma100_mean" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="sigmaD_minor" & t$gamma=="gamma100_mean" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , names.arg=dom_values, xlab=expression(paste(italic("h"))), cex.lab=1.8, cex.axis=1.5, cex.names=1.7, ylim=c(y_min_minor, y_max_minor))
error.bar(my_barplot, data_mean, data_se)


#sigmaD_minor; gamma1000_mean
t_sub_mean <- t[which(t$summary=="sigmaD_minor" & t$gamma=="gamma1000_mean" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="sigmaD_minor" & t$gamma=="gamma1000_mean" & t$type=="SE"),]
data_mean <- rbind(t_sub_mean$additive_gene,t_sub_mean$additive_site)
data_se <- rbind(t_sub_se$additive_gene,t_sub_se$additive_site)
my_barplot <- barplot(data_mean , beside=T , legend.text=T,col=v_colors , names.arg=dom_values, xlab=expression(paste(italic("h"))), cex.lab=1.8, cex.axis=1.5, cex.names=1.7, ylim=c(y_min_minor, y_max_minor))
error.bar(my_barplot, data_mean, data_se)
legend("topright", legend=c("gene", "sites"), col=v_colors, fill=v_colors, cex=1.5, box.lty=0)


>> save as 8.5 x 5.0 (landscape)