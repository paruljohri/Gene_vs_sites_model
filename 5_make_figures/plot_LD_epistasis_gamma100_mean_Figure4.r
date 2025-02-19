#to make a plot of LD stats with epistasis:

setwd("OneDrive - University of North Carolina at Chapel Hill/Work/Projects/Fitness_note_Brian/")
options(scipen=999)

t <- read.table("Tables/thousand_site/theta0_005/Epistasis/summary_LD_v2_gamma100_mean_lowrec_1_100.txt", h=T)
    
#A function to add arrows on the chart
error.bar <- function(x, y, upper, lower=upper, length=0.05, color){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, col=color)
}

epsilon_values <- c(0.00, 0.02, 0.04, 0.08, 0.16, 0.32)
v_colors <- c("firebrick2" , "gray39")

#One figure with D (derived) and D (minor)
par(mfcol=c(3,2))
par(mar=c(5,5,2,1))

#D (derived):
#getting limits to y-axis for D_derived
y_min <- min(c(t$additive_gene[which(t$summary=="D_derived" & t$type=="mean")]-t$additive_gene[which(t$summary=="D_derived" & t$type=="SE")], t$additive_site[which(t$summary=="D_derived" & t$type=="mean")]-t$additive_site[which(t$summary=="D_derived" & t$type=="SE")]))
y_max <- max(c(t$additive_gene[which(t$summary=="D_derived" & t$type=="mean")]+t$additive_gene[which(t$summary=="D_derived" & t$type=="SE")], t$additive_site[which(t$summary=="D_derived" & t$type=="mean")]+t$additive_site[which(t$summary=="D_derived" & t$type=="SE")]))


#h=0.0
t_sub_mean <- t[which(t$summary=="D_derived" & t$dominance=="h_0_0" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="D_derived" & t$dominance=="h_0_0" & t$type=="SE"),]
my_plot <- plot(x=epsilon_values, y=t_sub_mean$additive_gene, type="b", col=v_colors[1] , ylab=expression(paste(italic('D'), " (selected)")), xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(y_min,y_max))
points(x=epsilon_values, y=t_sub_mean$additive_site, col=v_colors[2], type="b")
error.bar(x=epsilon_values, y=t_sub_mean$additive_gene, upper=t_sub_se$additive_gene, color=v_colors[1])
error.bar(x=epsilon_values, y=t_sub_mean$additive_site, upper=t_sub_se$additive_site, color=v_colors[2])
title(main = expression(paste(italic("h"), "=0.0")), cex.main=2, adj=0)

#h=0.2
t_sub_mean <- t[which(t$summary=="D_derived" & t$dominance=="h_0_2" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="D_derived" & t$dominance=="h_0_2" & t$type=="SE"),]
my_plot <- plot(x=epsilon_values, y=t_sub_mean$additive_gene, type="b", col=v_colors[1] , ylab=expression(paste(italic('D'), " (selected)")), xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(y_min,y_max))
points(x=epsilon_values, y=t_sub_mean$additive_site, col=v_colors[2], type="b")
error.bar(x=epsilon_values, y=t_sub_mean$additive_gene, upper=t_sub_se$additive_gene, color=v_colors[1])
error.bar(x=epsilon_values, y=t_sub_mean$additive_site, upper=t_sub_se$additive_site, color=v_colors[2])
title(main = expression(paste(italic("h"), "=0.2")), cex.main=2, adj=0)

#h=0.5
t_sub_mean <- t[which(t$summary=="D_derived" & t$dominance=="h_0_5" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="D_derived" & t$dominance=="h_0_5" & t$type=="SE"),]
my_plot <- plot(x=epsilon_values, y=t_sub_mean$additive_gene, type="b", col=v_colors[1] , ylab=expression(paste(italic('D'), " (selected)")), xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(y_min,y_max))
points(x=epsilon_values, y=t_sub_mean$additive_site, col=v_colors[2], type="b")
error.bar(x=epsilon_values, y=t_sub_mean$additive_gene, upper=t_sub_se$additive_gene, color=v_colors[1])
error.bar(x=epsilon_values, y=t_sub_mean$additive_site, upper=t_sub_se$additive_site, color=v_colors[2])
title(main = expression(paste(italic("h"), "=0.5")), cex.main=2, adj=0)

#D (minor):
#getting limits to y-axis for D_minor
y_min <- min(c(t$additive_gene[which(t$summary=="D_minor" & t$type=="mean")]-t$additive_gene[which(t$summary=="D_minor" & t$type=="SE")], t$additive_site[which(t$summary=="D_minor" & t$type=="mean")]-t$additive_site[which(t$summary=="D_minor" & t$type=="SE")]))
y_max <- max(c(t$additive_gene[which(t$summary=="D_minor" & t$type=="mean")]+t$additive_gene[which(t$summary=="D_minor" & t$type=="SE")], t$additive_site[which(t$summary=="D_minor" & t$type=="mean")]+t$additive_site[which(t$summary=="D_minor" & t$type=="SE")]))

#h=0.0
t_sub_mean <- t[which(t$summary=="D_minor" & t$dominance=="h_0_0" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="D_minor" & t$dominance=="h_0_0" & t$type=="SE"),]
my_plot <- plot(x=epsilon_values, y=t_sub_mean$additive_gene, type="b", col=v_colors[1] , ylab=expression(paste(italic('D'), " (minor)")), xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(y_min,y_max))
points(x=epsilon_values, y=t_sub_mean$additive_site, col=v_colors[2], type="b")
error.bar(x=epsilon_values, y=t_sub_mean$additive_gene, upper=t_sub_se$additive_gene, color=v_colors[1])
error.bar(x=epsilon_values, y=t_sub_mean$additive_site, upper=t_sub_se$additive_site, color=v_colors[2])
legend("topright", legend=c("gene", "sites"), lty=1, col=v_colors, cex=1.5, box.lty=0, horiz=T)
title(main = expression(paste(italic("h"), "=0.0")), cex.main=2, adj=0)

#h=0.2
t_sub_mean <- t[which(t$summary=="D_minor" & t$dominance=="h_0_2" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="D_minor" & t$dominance=="h_0_2" & t$type=="SE"),]
my_plot <- plot(x=epsilon_values, y=t_sub_mean$additive_gene, type="b", col=v_colors[1] , ylab=expression(paste(italic('D'), " (minor)")), xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(y_min,y_max))
points(x=epsilon_values, y=t_sub_mean$additive_site, col=v_colors[2], type="b")
error.bar(x=epsilon_values, y=t_sub_mean$additive_gene, upper=t_sub_se$additive_gene, color=v_colors[1])
error.bar(x=epsilon_values, y=t_sub_mean$additive_site, upper=t_sub_se$additive_site, color=v_colors[2])
title(main = expression(paste(italic("h"), "=0.2")), cex.main=2, adj=0)

#h=0.5
t_sub_mean <- t[which(t$summary=="D_minor" & t$dominance=="h_0_5" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="D_minor" & t$dominance=="h_0_5" & t$type=="SE"),]
my_plot <- plot(x=epsilon_values, y=t_sub_mean$additive_gene, type="b", col=v_colors[1] , ylab=expression(paste(italic('D'), " (minor)")), xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(y_min,y_max))
points(x=epsilon_values, y=t_sub_mean$additive_site, col=v_colors[2], type="b")
error.bar(x=epsilon_values, y=t_sub_mean$additive_gene, upper=t_sub_se$additive_gene, color=v_colors[1])
error.bar(x=epsilon_values, y=t_sub_mean$additive_site, upper=t_sub_se$additive_site, color=v_colors[2])
title(main = expression(paste(italic("h"), "=0.5")), cex.main=2, adj=0)







#Second figure with sigma (derived) and sigma (minor)
par(mfcol=c(3,2))
par(mar=c(5,5,2,1))

#sigma (derived):
#getting limits to y-axis for sigmaD_derived
y_min <- min(c(t$additive_gene[which(t$summary=="sigmaD_derived" & t$type=="mean")]-t$additive_gene[which(t$summary=="sigmaD_derived" & t$type=="SE")], t$additive_site[which(t$summary=="sigmaD_derived" & t$type=="mean")]-t$additive_site[which(t$summary=="sigmaD_derived" & t$type=="SE")]))
y_max <- max(c(t$additive_gene[which(t$summary=="sigmaD_derived" & t$type=="mean")]+t$additive_gene[which(t$summary=="sigmaD_derived" & t$type=="SE")], t$additive_site[which(t$summary=="sigmaD_derived" & t$type=="mean")]+t$additive_site[which(t$summary=="sigmaD_derived" & t$type=="SE")]))


#h=0.0
t_sub_mean <- t[which(t$summary=="sigmaD_derived" & t$dominance=="h_0_0" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="sigmaD_derived" & t$dominance=="h_0_0" & t$type=="SE"),]
my_plot <- plot(x=epsilon_values, y=t_sub_mean$additive_gene, type="b", col=v_colors[1] , ylab=expression(paste(italic(sigma)[d], " (selected)")), xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(y_min,y_max))
points(x=epsilon_values, y=t_sub_mean$additive_site, col=v_colors[2], type="b")
error.bar(x=epsilon_values, y=t_sub_mean$additive_gene, upper=t_sub_se$additive_gene, color=v_colors[1])
error.bar(x=epsilon_values, y=t_sub_mean$additive_site, upper=t_sub_se$additive_site, color=v_colors[2])
title(main = expression(paste(italic("h"), "=0.0")), cex.main=2, adj=0)

#h=0.2
t_sub_mean <- t[which(t$summary=="sigmaD_derived" & t$dominance=="h_0_2" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="sigmaD_derived" & t$dominance=="h_0_2" & t$type=="SE"),]
my_plot <- plot(x=epsilon_values, y=t_sub_mean$additive_gene, type="b", col=v_colors[1] , ylab=expression(paste(italic(sigma)[d], " (selected)")), xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(y_min,y_max))
points(x=epsilon_values, y=t_sub_mean$additive_site, col=v_colors[2], type="b")
error.bar(x=epsilon_values, y=t_sub_mean$additive_gene, upper=t_sub_se$additive_gene, color=v_colors[1])
error.bar(x=epsilon_values, y=t_sub_mean$additive_site, upper=t_sub_se$additive_site, color=v_colors[2])
title(main = expression(paste(italic("h"), "=0.2")), cex.main=2, adj=0)

#h=0.5
t_sub_mean <- t[which(t$summary=="sigmaD_derived" & t$dominance=="h_0_5" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="sigmaD_derived" & t$dominance=="h_0_5" & t$type=="SE"),]
my_plot <- plot(x=epsilon_values, y=t_sub_mean$additive_gene, type="b", col=v_colors[1] , ylab=expression(paste(italic(sigma)[d], " (selected)")), xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(y_min,y_max))
points(x=epsilon_values, y=t_sub_mean$additive_site, col=v_colors[2], type="b")
error.bar(x=epsilon_values, y=t_sub_mean$additive_gene, upper=t_sub_se$additive_gene, color=v_colors[1])
error.bar(x=epsilon_values, y=t_sub_mean$additive_site, upper=t_sub_se$additive_site, color=v_colors[2])
title(main = expression(paste(italic("h"), "=0.5")), cex.main=2, adj=0)


#sigma (minor):
#getting limits to y-axis for sigma_minor
y_min <- min(c(t$additive_gene[which(t$summary=="sigmaD_minor" & t$type=="mean")]-t$additive_gene[which(t$summary=="sigmaD_minor" & t$type=="SE")], t$additive_site[which(t$summary=="sigmaD_minor" & t$type=="mean")]-t$additive_site[which(t$summary=="sigmaD_minor" & t$type=="SE")]))
y_max <- max(c(t$additive_gene[which(t$summary=="sigmaD_minor" & t$type=="mean")]+t$additive_gene[which(t$summary=="sigmaD_minor" & t$type=="SE")], t$additive_site[which(t$summary=="sigmaD_minor" & t$type=="mean")]+t$additive_site[which(t$summary=="sigmaD_minor" & t$type=="SE")]))

#h=0.0
t_sub_mean <- t[which(t$summary=="sigmaD_minor" & t$dominance=="h_0_0" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="sigmaD_minor" & t$dominance=="h_0_0" & t$type=="SE"),]
my_plot <- plot(x=epsilon_values, y=t_sub_mean$additive_gene, type="b", col=v_colors[1] , ylab=expression(paste(italic(sigma)[d], " (minor)")), xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(y_min,y_max))
points(x=epsilon_values, y=t_sub_mean$additive_site, col=v_colors[2], type="b")
error.bar(x=epsilon_values, y=t_sub_mean$additive_gene, upper=t_sub_se$additive_gene, color=v_colors[1])
error.bar(x=epsilon_values, y=t_sub_mean$additive_site, upper=t_sub_se$additive_site, color=v_colors[2])
title(main = expression(paste(italic("h"), "=0.0")), cex.main=2, adj=0)
legend("topright", legend=c("gene", "sites"), lty=1, col=v_colors, cex=1.5, box.lty=0, horiz=T) #horiz=T

#h=0.2
t_sub_mean <- t[which(t$summary=="sigmaD_minor" & t$dominance=="h_0_2" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="sigmaD_minor" & t$dominance=="h_0_2" & t$type=="SE"),]
my_plot <- plot(x=epsilon_values, y=t_sub_mean$additive_gene, type="b", col=v_colors[1] , ylab=expression(paste(italic(sigma)[d], " (minor)")), xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(y_min,y_max))
points(x=epsilon_values, y=t_sub_mean$additive_site, col=v_colors[2], type="b")
error.bar(x=epsilon_values, y=t_sub_mean$additive_gene, upper=t_sub_se$additive_gene, color=v_colors[1])
error.bar(x=epsilon_values, y=t_sub_mean$additive_site, upper=t_sub_se$additive_site, color=v_colors[2])
title(main = expression(paste(italic("h"), "=0.2")), cex.main=2, adj=0)

#h=0.5
t_sub_mean <- t[which(t$summary=="sigmaD_minor" & t$dominance=="h_0_5" & t$type=="mean"),]
t_sub_se <- t[which(t$summary=="sigmaD_minor" & t$dominance=="h_0_5" & t$type=="SE"),]
my_plot <- plot(x=epsilon_values, y=t_sub_mean$additive_gene, type="b", col=v_colors[1] , ylab=expression(paste(italic(sigma)[d], " (minor)")), xlab="epistasis coefficient", cex.lab=2, cex.axis=1.5, ylim=c(y_min,y_max))
points(x=epsilon_values, y=t_sub_mean$additive_site, col=v_colors[2], type="b")
error.bar(x=epsilon_values, y=t_sub_mean$additive_gene, upper=t_sub_se$additive_gene, color=v_colors[1])
error.bar(x=epsilon_values, y=t_sub_mean$additive_site, upper=t_sub_se$additive_site, color=v_colors[2])
title(main = expression(paste(italic("h"), "=0.5")), cex.main=2, adj=0)
