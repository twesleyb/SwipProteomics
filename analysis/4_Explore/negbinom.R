#!/usr/bin/env Rscript

# title: "AML_course_libs.R"
# author: tyler w a bradshaw
# reference: http://sherrytowers.com/2018/04/10/negative-binomial-likelihood-fits-for-overdispersed-count-data/

# true model for log of y is linear in x
a = 4.00
b = 0.01
x = seq(0,100,0.01)
logy = a+b*x
ypred = exp(logy)

# randomly generate Poisson distributed
# and Negative Binomially distributed
# simulated data about the true model
set.seed(314300)

yobs = rpois(length(ypred),ypred)
alpha = 0.2
yobs_nb = my_rnbinom(length(ypred),ypred,alpha)

# plot the data
#ylim = c(min(c(yobs,yobs_nb)),max(c(yobs,yobs_nb)))
#plot(x,yobs,cex=1.0,xlab="x",ylab="y",ylim=ylim)
#points(x,yobs_nb,cex=1.2,col=2)
#points(x,yobs,cex=1.0)
#lines(x,ypred,col=3,lwd=5)
#legend("topleft",legend=c("Poisson distributed data","Negative Binomial distributed data","True model"),col=c(1,2,3),lwd=5,cex=0.9,bty="n")
library(ggplot2)
