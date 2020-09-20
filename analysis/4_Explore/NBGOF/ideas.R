# https://stats.stackexchange.com/questions/205403/fitting-negative-binomial-distribution-to-large-count-data/368656#368656

# read the file containing count data
data <- read.csv("data.txt", header=FALSE)

# plot the histogram
hist(data[[1]], prob=TRUE, breaks=145)

# load library
library(fitdistrplus)

# fit the negative binomial distribution
fit <- fitdist(data[[1]], "nbinom")

# get the fitted densities. mu and size from fit.
fitD <- dnbinom(0:145, size=25.05688, mu=31.56127)

# add fitted line (blue) to histogram
lines(fitD, lwd="3", col="blue")

# Goodness of fit with the chi squared test
# get the frequency table
t <- table(data[[1]])

# convert to dataframe
df <- as.data.frame(t)

# get frequencies
observed_freq <- df$Freq

# perform the chi-squared test
chisq.test(observed_freq, p=fitD)


# also try:
#https://stats.idre.ucla.edu/r/dae/negative-binomial-regression/
