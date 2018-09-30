#makes statistical comparisons for these four genes, according to the google doc
#prints them out in text files
library(lme4)
library(fitdistrplus)

#these are the gene directories - there are 4 of them, hard-coded
#within each directory, there are ONLY the csv files obtained from reference_analyses
Frr <- "/Users/yancyliao/Dropbox/bioinformatics/Shakya_analyses/Frr"
RplB <- "/Users/yancyliao/Dropbox/bioinformatics/Shakya_analyses/RplB"
RplK <- "/Users/yancyliao/Dropbox/bioinformatics/Shakya_analyses/RplK"
RpsT <- "/Users/yancyliao/Dropbox/bioinformatics/Shakya_analyses/RpsT"

#get the file names for each gene
setwd(Frr)
files1 <- list.files(full.names=F, recursive=FALSE)
setwd(RplB)
files2 <- list.files(full.names=F, recursive=FALSE)
setwd(RplK)
files3 <- list.files(full.names=F, recursive=FALSE)
setwd(RpsT)
files4 <- list.files(full.names=F, recursive=FALSE)

#set up all the vectors to store the data
names <- ANOVA_p_values <- c()
AIC_glm <- AIC_lme <- AIC_lm <- is_glm_lowest <- c()
indices1 <- indices2 <- indices3 <-indices4 <- c()
anova_p <- c()
r_nbinom1 <- r_nbinom2 <- r_nbinom3 <- r_nbinom4 <- c()
glm_ste_fixed<- glm_est_fixed<-rand_std<-c()
random1<-random2<-random3<-random4<-c()
theta <-c()
se_nbinom <- rse_nbinom <-c()
mean <- std <- c()

anova_p1 <- anova_p2 <- anova_p3 <- anova_p4 <- anova_p5 <- anova_p6 <- c()

#only look at genomes for which each gene has data
for(i in 1:length(files1)){
  if(!grepl('all_together', files1[i]) & 
     (grepl("Acido", files1[i]) | #changed
     grepl("Akker", files1[i]) | #changed
     grepl("Bordetella", files1[i]) |
     grepl("Burkholder", files1[i]) | #missing, originally?
     grepl("Chloroflexus", files1[i]) |
     grepl("Deinococc", files1[i]) |
     grepl("Dictyo", files1[i]) | #typo
     grepl("Entero", files1[i]) |
     grepl("Fusob", files1[i]) | #changed
     grepl("Gemma", files1[i]) | #changed
     grepl("Geoba", files1[i]) |
     grepl("Herpeto", files1[i]) |
     grepl("Hydrogen", files1[i]) |
     grepl("Lepto", files1[i]) |
     grepl("Nitroso", files1[i]) |
     grepl("Nostoc", files1[i]) |
     grepl("Persephon", files1[i]) |
     grepl("Porphyro", files1[i]) | 
     grepl("Rhodo", files1[i]) |
     grepl("Clostri", files1[i]) | #not there.....
     grepl("Shewan", files1[i]) | #Frr which strain? not the all_together one right?
     grepl("Thermoanaer", files1[i]) |
     grepl("Thermus", files1[i]) |
     grepl("Trepon", files1[i]) |
     grepl("Woline", files1[i]) | #changed
     grepl("Zymom", files1[i])) & #changed
     (files1[i] %in% files2) & 
    (files1[i] %in% files3) & 
    (files1[i] %in% files4)) {
    
    #store the file name for genome
    names <- c(names, files1[i])
    
    #store the file-name indices of genome in case they're needed
    ind2 <- match(files1[i], files2)
    ind3 <- match(files1[i], files3)
    ind4 <- match(files1[i], files4)
    indices1 <- c(indices1, i)
    indices2<-c(indices2, ind2)
    indices3<-c(indices3, ind3)
    indices4<-c(indices4, ind4)
    
    #get the coverage data for each gene for genome
    setwd(Frr)
    mydata1 <- read.csv(files1[i], header = FALSE)
    setwd(RplB)
    mydata2 <- read.csv(files2[ind2], header = FALSE)
    setwd(RplK)
    mydata3 <- read.csv(files3[ind3], header = FALSE)
    setwd(RpsT)
    mydata4 <- read.csv(files4[ind4], header = FALSE)
    
    #turn coverage data for each gene into appropriate data vector
    data_vector1 <- rep(mydata1[[1]], mydata1[[2]])
    data_vector2 <- rep(mydata2[[1]], mydata2[[2]])
    data_vector3 <- rep(mydata3[[1]], mydata3[[2]])
    data_vector4 <- rep(mydata4[[1]], mydata4[[2]])
    
    #fit a negative binomial distribution for the coverage data for each gene
    fit.nbinom <- fitdistr(data_vector1, "negative binomial")
    nbinom_quantiles <- qnbinom(ppoints(length(data_vector1)), size=fit.nbinom$estimate[1], mu =fit.nbinom$estimate[2])
    r_nbinom1 <- c(r_nbinom1, cor(nbinom_quantiles, data_vector1))
    fit.nbinom <- fitdistr(data_vector2, "negative binomial")
    nbinom_quantiles <- qnbinom(ppoints(length(data_vector2)), size=fit.nbinom$estimate[1], mu =fit.nbinom$estimate[2])
    r_nbinom2 <- c(r_nbinom2, cor(nbinom_quantiles, data_vector2))
    fit.nbinom <- fitdistr(data_vector3, "negative binomial")
    nbinom_quantiles <- qnbinom(ppoints(length(data_vector3)), size=fit.nbinom$estimate[1], mu =fit.nbinom$estimate[2])
    r_nbinom3 <- c(r_nbinom3, cor(nbinom_quantiles, data_vector3))
    fit.nbinom <- fitdistr(data_vector4, "negative binomial")
    nbinom_quantiles <- qnbinom(ppoints(length(data_vector4)), size=fit.nbinom$estimate[1], mu =fit.nbinom$estimate[2])
    r_nbinom4 <- c(r_nbinom4, cor(nbinom_quantiles, data_vector4))
    
    #combine all the coverage data into one vector to get overall statistics from
    coverage<-c(data_vector1, data_vector2, data_vector3, data_vector4)
    mean<-c(mean, mean(coverage))
    std<-c(std, sd(coverage))
    
    #get the overall negative binomial distribution
    fit.nbinom <- fitdistr(coverage, "negative binomial")
    se_nbinom <- c(se_nbinom, fit.nbinom$sd[2])
    rse_nbinom <- c(rse_nbinom, fit.nbinom$sd[2]/mean(coverage))
    
    #give gene labels for each coverage datapoint, then put together into data frame
    gene<-c(rep("Frr", length(data_vector1)), 
            rep("RplB", length(data_vector2)), 
            rep("RplK", length(data_vector3)), 
            rep("RpsT", length(data_vector4)))
    shakya_data = data.frame(gene, coverage)
    
    #perform anova analysis by the 4 genes groups, store the p-value
    anova <- aov(coverage ~ gene, data = shakya_data)
    anova_p <- c(anova_p, summary(anova)[[1]][["Pr(>F)"]][1])
    
    coverage1 <- c(data_vector1, data_vector2)
    coverage2 <- c(data_vector2, data_vector3)
    coverage3 <- c(data_vector3, data_vector4)
    coverage4 <- c(data_vector4, data_vector1)
    coverage5 <- c(data_vector1, data_vector3)
    coverage6 <- c(data_vector2, data_vector4)
    gene1<-c(rep("Frr", length(data_vector1)), 
            rep("RplB", length(data_vector2)))
    gene2<-c(rep("RplB", length(data_vector2)), 
            rep("RplK", length(data_vector3)))
    gene3<-c(rep("RplK", length(data_vector3)), 
            rep("RpsT", length(data_vector4)))
    gene4<-c(rep("Frr", length(data_vector1)),
            rep("RpsT", length(data_vector4)))
    gene5<-c(rep("Frr", length(data_vector1)),
            rep("RplK", length(data_vector3)))
    gene6<-c(rep("RplB", length(data_vector2)), 
             rep("RpsT", length(data_vector4)))
    
    shakya_data1 = data.frame(gene1, coverage1)
    shakya_data2 = data.frame(gene2, coverage2)
    shakya_data3 = data.frame(gene3, coverage3)
    shakya_data4 = data.frame(gene4, coverage4)
    shakya_data5 = data.frame(gene5, coverage5)
    shakya_data6 = data.frame(gene6, coverage6)
    
    
    anova1 <- aov(coverage1 ~ gene1, data = shakya_data1)
    anova_p1 <- c(anova_p1, summary(anova1)[[1]][["Pr(>F)"]][1])
    anova2 <- aov(coverage2 ~ gene2, data = shakya_data2)
    anova_p2 <- c(anova_p2, summary(anova2)[[1]][["Pr(>F)"]][1])
    anova3 <- aov(coverage3 ~ gene3, data = shakya_data3)
    anova_p3 <- c(anova_p3, summary(anova3)[[1]][["Pr(>F)"]][1])
    anova4 <- aov(coverage4 ~ gene4, data = shakya_data4)
    anova_p4 <- c(anova_p4, summary(anova4)[[1]][["Pr(>F)"]][1])
    anova5 <- aov(coverage5 ~ gene5, data = shakya_data5)
    anova_p5 <- c(anova_p5, summary(anova5)[[1]][["Pr(>F)"]][1])
    anova6 <- aov(coverage6 ~ gene6, data = shakya_data6)
    anova_p6 <- c(anova_p6, summary(anova6)[[1]][["Pr(>F)"]][1])
    
    #fit a linear regression model with gene as the predictor variable, record AIC
    lm.model = lm(coverage ~ gene, data=shakya_data)
    AIC_lm <- c(AIC_lm, AIC(lm.model))
    
    #fit a linear mixed effects model with gene as random variable, record AIC
    lme.model = lmer(coverage ~ (1|gene), data = shakya_data)
    AIC_lme <- c(AIC_lme, AIC(lme.model))
    
    #fit a generalized linear mixed effects model with random gene effects
    #and negative binomial error, record AIC
    glm.model = glmer.nb(coverage ~ (1|gene), data=shakya_data)
    AIC_glm <- c(AIC_glm, AIC(glm.model))
    
    #get the dispersion parameter theta
    theta < c(theta, getME(glm.model, "glmer.nb.theta"))
    
    #get the standard deviation of gene random effects
    rand_std <- c(rand_std, as.data.frame(VarCorr(glm.model))$sdcor)
    
    #get the estimated *random effects*, not estimated individual means, for each gene
    random1<-c(random1, getME(glm.model, "b")[1])
    random2<-c(random2, getME(glm.model, "b")[2])
    random3<-c(random3, getME(glm.model, "b")[3])
    random4<-c(random4, getME(glm.model, "b")[4])
    
    #get the estimated overall mean and standard error of that estimate for the model
    glm_ste_fixed = c(glm_ste_fixed, summary(glm.model)[[10]][[2]])
    glm_est_fixed = c(glm_est_fixed, summary(glm.model)[[10]][[1]])
    
    #boolean of whether GLMM has lower AIC score than the less complex models
    is_glm_lowest <- c(is_glm_lowest, (AIC(glm.model) < AIC(lm.model))
                    & (AIC(glm.model) < AIC(lme.model)))
    
    #print out the index that program is currently working on, to keep track of progress
    cat("working on genome", i, "\n")
  }
  
}

#do conversions of extracted parameters from log space to response scale
#careful: some are straightforward, others need to be interpreted narrowly
glm_mean_converted <- exp(glm_est_fixed)
glm_ste_converted <- exp(glm_est_fixed + glm_ste_fixed) - exp(glm_est_fixed)
rand_std_converted <- exp(glm_est_fixed + rand_std) - exp(glm_est_fixed)
random1_mean_converted <- exp(glm_est_fixed + random1)
random2_mean_converted <-exp(glm_est_fixed + random2)
random3_mean_converted <-exp(glm_est_fixed + random3)
random4_mean_converted <-exp(glm_est_fixed + random4)

#compile results into 5 data frames like a table
results1 = data.frame(names, mean, std, anova_p)
results2 = data.frame(names, r_nbinom1, r_nbinom2, r_nbinom3, r_nbinom4)
results3 = data.frame(names, glm_mean_converted, glm_ste_converted)
results4 = data.frame(names, random1_mean_converted, random2_mean_converted, 
                     random3_mean_converted, random4_mean_converted, rand_std_converted)
results5 = data.frame(names,AIC_glm, AIC_lme, AIC_lm, is_glm_lowest)
results6 = data.frame(names, se_nbinom, rse_nbinom)

#save results into files
setwd("/Users/yancyliao/Dropbox/bioinformatics/Shakya_analyses/")
cat(capture.output(results1), file = 'results1.txt', sep = '\n')
cat(capture.output(results2), file = 'results2.txt', sep = '\n')
cat(capture.output(results3), file = 'results3.txt', sep = '\n')
cat(capture.output(results4), file = 'results4.txt', sep = '\n')
cat(capture.output(results5), file = 'results5.txt', sep = '\n')
cat(capture.output(results6), file = 'results6.txt', sep = '\n')