################################################################################
##	Project: 
##    ``Heterogeneous Treatment Effects in Regression Discontinuity Designs''
##	  by Calonico, Cattaneo, Farrell, Palomba, and Titiunik
##
##	Purpose: 
##    Replicate and extend the heterogeneity analysis of AER 2022, 112(2): 442–493
##    by Akhtari, Moreira, and Trucco, specifically heterogeneity by income
##    as reported in the original Table A.21.
##    https://doi.org/10.3886/E150323V1
##
##########################################################################


rm(list=ls())
install.packages("rdhte")
library(rdhte)
library(Hmisc)
sessionInfo()

##########################################################################
## Read in data and set up variables

data <- read.csv("AMT_2022_AER.csv")

# Set up main variables
Y <- data$expthisschl_lessthan2_DPB
X <- -1*data$pX
cluster = data$COD_MUNICIPIO
W <- data$medianincome2000
trunc <- t1_1 <- quantile(W, na.rm=TRUE, p=(0.99))
W[W>=trunc] <- trunc
W <- W/100
W2 <- W^2

##########################################################################
## Discretize income by different quantiles

# Median
quant_med <- quantile(W[data$tag==1], prob = seq(0, 1, length = 3), type = 5, na.rm=TRUE)
W_med <- findInterval(W, quant_med, rightmost.closed = TRUE)

# Quartiles
quant_qrt <- quantile(W[data$tag==1], prob = seq(0, 1, length = 5), type = 5, na.rm=TRUE)
W_qrt <- findInterval(W, quant_qrt, rightmost.closed = TRUE)

# Deciles
quant_dec <- quantile(W[data$tag==1], prob = seq(0, 1, length = 11), type = 5, na.rm=TRUE)
W_dec <- findInterval(W, quant_dec, rightmost.closed = TRUE)



##########################################################################
## Estimation and inference

# Standard RD treatment effect (no heterogeneity)
rd_ate <- rdhte(y = Y, x = X, cluster = cluster)

# HTE - Binary
rd_med <- rdhte(y = Y, x = X, covs.hte = factor(W_med), cluster = cluster) 

# HTE - Quartiles
rd_qrt <- rdhte(y = Y, x=X, covs.hte = factor(W_qrt), cluster = cluster) 

# HTE - Deciles
rd_dec <- rdhte(y = Y, x = X, covs.hte = factor(W_dec), cluster = cluster)

# HTE - Continuous, Linear in Income
rd_linear <- rdhte(y = Y, x = X, covs.hte = W, cluster = cluster)

# HTE - Continuous, Quadratic in Income
rd_quad <- rdhte(y = Y, x = X, covs.hte = cbind(W,W2), cluster = cluster)




##########################################################################
## Export results for Table 1


ci_1 <- formatC(cbind(rd_ate$ci_rb[,1], rd_ate$ci_rb[,2]),format = "f", digits = 3)
ci_2 <- formatC(cbind(rd_med$ci_rb[,1], rd_med$ci_rb[,2]),format = "f", digits = 3)
ci_3 <- formatC(cbind(rd_qrt$ci_rb[,1], rd_qrt$ci_rb[,2]),format = "f", digits = 3)
ci_4 <- formatC(cbind(rd_dec$ci_rb[,1], rd_dec$ci_rb[,2]),format = "f", digits = 3)
ci_5 <- formatC(cbind(rd_linear$ci_rb[,1], rd_linear$ci_rb[,2]),format = "f", digits = 3)
ci_6 <- formatC(cbind(rd_quad$ci_rb[,1], rd_quad$ci_rb[,2]),format = "f", digits = 3)

table1_est <- formatC(rd_ate$Estimate, format = "f", digits = 3)
table2_est <- formatC(rd_med$Estimate, format = "f", digits = 3)
table3_est <- formatC(rd_qrt$Estimate, format = "f", digits = 3)
table4_est <- formatC(rd_dec$Estimate, format = "f", digits = 3)
table5_est <- formatC(rd_linear$Estimate, format = "f", digits = 3)
table6_est <- formatC(rd_quad$Estimate, format = "f", digits = 3)

table1_ci <- paste("[",ci_1[,1]," ; ",ci_1[,2],"]",sep="")
table2_ci <- paste("[",ci_2[,1]," ; ",ci_2[,2],"]",sep="")
table3_ci <- paste("[",ci_3[,1]," ; ",ci_3[,2],"]",sep="")
table4_ci <- paste("[",ci_4[,1]," ; ",ci_4[,2],"]",sep="")
table5_ci <- paste("[",ci_5[,1]," ; ",ci_5[,2],"]",sep="")
table6_ci <- paste("[",ci_6[,1]," ; ",ci_6[,2],"]",sep="")

table1_pv <- formatC(rd_ate$pv, format = "f", digits = 3)
table2_pv <- formatC(rd_med$pv,   format = "f", digits = 3)
table3_pv <- formatC(rd_qrt$pv,   format = "f", digits = 3)
table4_pv <- formatC(rd_dec$pv,   format = "f", digits = 3)
table5_pv <- formatC(rd_linear$pv,   format = "f", digits = 3)
table6_pv <- formatC(rd_quad$pv,   format = "f", digits = 3)

table1_N <- formatC(rd_ate$Nh,  format="f", big.mark=",", digits=0)
table2_N <- formatC(rd_med$Nh,    format="f", big.mark=",", digits=0)
table3_N <- formatC(rd_qrt$Nh,    format="f", big.mark=",", digits=0)
table4_N <- formatC(rd_dec$Nh,    format="f", big.mark=",", digits=0)
table5_N <- formatC(rd_linear$Nh,    format="f", big.mark=",", digits=0)
table6_N <- formatC(rd_quad$Nh,    format="f", big.mark=",", digits=0)

table1_h <- formatC(rd_ate$h, format="f", big.mark=",", digits=3)
table2_h <- formatC(rd_med$h,    format="f", big.mark=",", digits=3)
table3_h <- formatC(rd_qrt$h,    format="f", big.mark=",", digits=3)
table4_h <- formatC(rd_dec$h,    format="f", big.mark=",", digits=3)
table5_h <- formatC(rd_linear$h,    format="f", big.mark=",", digits=3)
table6_h <- formatC(rd_quad$h,    format="f", big.mark=",", digits=3)

t1 <- cbind(table1_est, table1_ci, table1_pv, table1_N, table1_h)
t2 <- cbind(table2_est, table2_ci, table2_pv, table2_N, table2_h)
t3 <- cbind(table3_est, table3_ci, table3_pv, table3_N, table3_h)
t4 <- cbind(table4_est, table4_ci, table4_pv, table4_N, table4_h)
t5 <- cbind(table5_est, table5_ci, table5_pv, table5_N, table5_h)
t6 <- cbind(table6_est, table6_ci, table6_pv, table6_N, table6_h)

TT <- rbind(t1, NA, t2,NA, t3,NA, t4,NA, t5,NA, t6)

colnames(TT) <- c("Estimate",  "95\\% Robust CI","$p$-value", "N", "$h$")
rownames(TT) <- c("$\\tau$",
                 "\\emph{Binary}", "\\qquad Below median", "\\qquad Above median", 
                 "\\emph{By Quartile}:", paste("\\qquad",1:4,sep=""), 
                 "\\emph{By Decile}:", paste("\\qquad",1:10,sep=""), 
                 "\\emph{Linear}:", "$\\theta(0)$", "$\\xi(0)$", 
                 "\\emph{Quadratic}:", "$\\theta(0)$", "$\\xi_1(0)$", "$\\xi_2(0)$"
                 )


table1_tex <- latex(TT, file = paste("CCFPT_2025_HTERD_Table1_RAW_OUTPUT",".txt",sep=""), 
                   landscape=FALSE, outer.size='scriptsize', col.just=rep('c',ncol(t1)), 
                   n.rgroup=c(1,19,3,4), 
                   rgroup = c("Overall Average Effect", "Discretize Income", "Linear", "Quadratic"),
                   center='none', title='', table.env=FALSE)







##########################################################################
## Figure 1 - with linear & quadratic fits

# Median
pdf(file = "CCFPT_2025_HTERD_Figure1a.pdf",  width = 4, height = 6)
  curve(rd_linear$Estimate[1] + rd_linear$Estimate[2]*x, from=min(W, na.rm=TRUE), to=max(W, na.rm=TRUE),  xlab="Income", ylab="Estimate", yaxt = "n", ylim=c(-0.75, 0.75), cex.axis=0.85)
  W.grid <- seq(from=min(W, na.rm=TRUE), to=max(W, na.rm=TRUE), length.out=500)
  lines(W.grid, rd_quad$Estimate[1] + rd_quad$Estimate[2]*W.grid + rd_quad$Estimate[3]*W.grid^2, lty="dashed", lwd=1.5)
  axis(2, at=c(-0.75, 0.75, 0, round(rd_ate$Estimate[1],3), round(rd_linear$Estimate[1],3)), labels = FALSE)
  axis(2, at=c(round(rd_ate$Estimate[1],3), round(rd_linear$Estimate[1],3)), cex.axis=0.85)
  axis(2, at=c(-0.75, 0.75), cex.axis=0.85)
  axis(2, at=c(0.0), cex.axis=0.85)
  abline(h=rd_ate$Estimate[1], lty=2)
  ci_l <- c(rd_med$ci_rb[,1])
  ci_r <- c(rd_med$ci_rb[,2])
  est <- rd_med$Estimate
  a <- quant_med
  x <- a[-length(a)] + diff(a)/2
  points(x = x, y=est, pch=20)
  arrows(x0=x, y0=ci_l, x1=x, y1=ci_r, code=3, angle=90, length=0.025)
dev.off()


# Quartiles
pdf(file = "CCFPT_2025_HTERD_Figure1b.pdf",  width = 4, height = 6)
  curve(rd_linear$Estimate[1] + rd_linear$Estimate[2]*x, from=min(W, na.rm=TRUE), to=max(W, na.rm=TRUE), xlab="Income", ylab=NA, yaxt = "n", ylim=c(-0.75, 0.75), cex.axis=0.85)
  W.grid <- seq(from=min(W, na.rm=TRUE), to=max(W, na.rm=TRUE), length.out=500)
  lines(W.grid, rd_quad$Estimate[1] + rd_quad$Estimate[2]*W.grid + rd_quad$Estimate[3]*W.grid^2, lty="dashed", lwd=1.5)
  axis(1, labels = FALSE, tick = FALSE)
  abline(h=rd_ate$Estimate[1], lty=2)
  ci_l <- c(rd_qrt$ci_rb[,1])
  ci_r <- c(rd_qrt$ci_rb[,2])
  est <- rd_qrt$Estimate
  a <- quant_qrt
  x <- a[-length(a)] + diff(a)/2
  points(x = x, y=est, pch=20)
  arrows(x0=x, y0=ci_l, x1=x, y1=ci_r, code=3, angle=90, length=0.025)
dev.off()


pdf(file = "CCFPT_2025_HTERD_Figure1c.pdf",  width = 4, height = 6)
  curve(rd_linear$Estimate[1] + rd_linear$Estimate[2]*x, from=min(W, na.rm=TRUE), to=max(W, na.rm=TRUE), xlab="Income", ylab=NA, yaxt = "n", ylim=c(-0.75, 0.75), cex.axis=0.85)
  W.grid <- seq(from=min(W, na.rm=TRUE), to=max(W, na.rm=TRUE), length.out=500)
  lines(W.grid, rd_quad$Estimate[1] + rd_quad$Estimate[2]*W.grid + rd_quad$Estimate[3]*W.grid^2, lty="dashed", lwd=1.5)
  abline(h=rd_ate$Estimate[1], lty=2)
  ci_l <- c(rd_dec$ci_rb[,1])
  ci_r <- c(rd_dec$ci_rb[,2])
  est <- rd_dec$Estimate
  a <- quant_dec
  x <- a[-length(a)] + diff(a)/2
  points(x = x, y=est, pch=20)
  arrows(x0=x, y0=ci_l, x1=x, y1=ci_r, code=3, angle=90, length=0.025)
dev.off()

  
  
  

