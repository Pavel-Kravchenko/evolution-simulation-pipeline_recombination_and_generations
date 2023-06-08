args <- commandArgs(trailingOnly = T)

name_plot <- toString(args[4])
mask <- toString(args[5])
name <- toString(args[6])

"%+%" <- function(...){
  paste0(...)
}

name_plot0.0 <- name_plot %+% "/LD.txt"
name_plot0.1 <- name_plot %+% "/r2.txt"

LD_txt <- read.csv(name_plot0.0, sep = '\t', 
                 stringsAsFactors = F)
LD_txt <- LD_txt[order(LD_txt$Len, decreasing = F),] 
x <- LD_txt$Len
y <- LD_txt$LD

r2_txt <- read.csv(name_plot0.1, sep = '\t', 
                 stringsAsFactors = F)
r2_txt <- r2_txt[order(r2_txt$Len, decreasing = F),] 
a <- r2_txt$Len
b <- r2_txt$r2

wind <- as.integer(round(length(x)/50, 0))
if (wind == 0) {wind <- 2}
name_plot1 <- name_plot %+% "/LD_plot_" %+% name %+% "_" %+% mask %+% "_" %+% wind  %+% ".png"
name_plot2 <- name_plot %+% "/r2_plot_" %+% name %+% "_" %+% mask %+% "_" %+% wind  %+% ".png"
name_plot3 <- name_plot %+% "/LD_hist_" %+% name %+% "_" %+% mask %+% "_" %+% wind  %+% ".png"
name_plot4 <- name_plot %+% "/Len_hist_" %+% name %+% "_" %+% mask %+% "_" %+% wind  %+% ".png"
name_plot5 <- name_plot %+% "/Mix_" %+% name %+% "_" %+% mask %+% "_" %+% wind  %+% ".png"

#----------------------------------------------------------------
png(file= name_plot1,width=13000, height=8000, res=800)
plot(x, y, col=grey(.7),  main = "LD plot of " %+% name %+% ". " %+% "Mask = " %+% mask %+% ". Window = " %+% wind, xlab ="Distance (nuc.)",
 ylab = "LD", ylim=c(0, max(y)), cex.main=1.85, cex.lab=2.2, cex.axis=2.2)
grid()
f <- rep(1/wind, wind)
y_lag <- filter(y, f, sides=1)
lines(x, y_lag, col="red",  lwd=10)
dev.off()

# png(file= name_plot2,width=13000, height=8000, res=800)
# plot(a, b, col=grey(.7),  main = "r**2 plot of " %+% name %+% ". " %+% "Mask = " %+% mask %+% ". Window = " %+% wind, 
#      xlab ="Distance", ylab = "r**2", cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
# grid()
# q <- rep(1/wind, wind)
# b_lag <- filter(b, q, sides=1)
# lines(a, b_lag, col="blue")
# #abline(v=1, col="orange", lty=2)
# dev.off()

# png(file= name_plot3, width=5000, height=3000, res=800)
# hist(y, breaks = "Sturges", main = "LD histogram of " %+% name %+% ". " %+% "Mask = " %+% mask, xlab ="LD", ylab = "Frequency", pch=16, col = "orange")
# dev.off()
# 
# png(file= name_plot4,width=5000, height=3000, res=1000)
# hist(x, breaks = "Sturges", main = "Len histogram of " %+% name %+% ". " %+% "Mask = " %+% mask, xlab ="Len", ylab = "Frequency", pch=16, col = "green")
# dev.off()
# 
# png(file= name_plot5, width=5000, height=3000, res=1000)
# par(mfrow=c(2,2))
# plot(x, y, col=grey(.7),  main = "LD plot", xlab ="Distance", ylab = "LD", ylim=c(0, max(y)))
# 
# f <- rep(1/wind, wind)
# y_lag <- filter(y, f, sides=1)
# lines(x, y_lag, col="red")
# 
# plot(a, b, col=grey(.7),  main = "r**2 plot", xlab ="Distance", ylab = "r**2")
# 
# q <- rep(1/wind, wind)
# b_lag <- filter(b, q, sides=1)
# lines(a, b_lag, col="blue")
# 
# hist(y, breaks = "Sturges", main = "LD histogram", xlab ="LD", ylab = "Frequency", pch=16, col = "orange")
# hist(x, breaks = "Sturges", main = "Len histogram", xlab ="Len", ylab = "Frequency", pch=16, col = "green") 
# dev.off()
# 
# #-------------------------------------------------------------------------
# df1 <- data.frame(x, y)
# border <- as.integer(round(length(x)/5, 0))
# if (border > 1) {df1 <- df1[1:border,]}
# 
# #if (length(df1$x) > 50000) {df1 <- df1[sample(nrow(df1), 50000), ]}
# sink("Kendall_test.txt")
# cor.test(as.integer(df1$y), as.integer(df1$x), method = "kendall")
# sink()
