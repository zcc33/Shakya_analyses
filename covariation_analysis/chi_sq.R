setwd("/Users/yancyliao/Dropbox/bioinformatics/covariation_analysis")
dfs <- read.csv("output.csv")
X <- rbinom(200, 1, 0.5)
Y <- rbinom(200, 1, 0.5)
for(i in 1:40000){
  chisq.test(X, Y, correct=FALSE)
}

