#########################################################
# (cc) 2012 Angelos Markos, amarkos [at] gmail [.] com  #
######################################################### 

#####################################################
# Chapter 7: CFA with Equality Constraints,         #
#            Multiple Groups, and Mean Structures   #
#####################################################

library(lavaan)

#####################################
# Equality Constraints p.236 - 252  #
#####################################

#Input matrix: the correlation matrix of Figure 7.1
cormat <- '
1.000
0.661  1.000
0.630  0.643  1.000
0.270  0.300  0.268  1.000
0.297  0.265  0.225  0.805  1.000
0.290  0.287  0.248  0.796  0.779  1.000'

Cmat <- getCov(cormat)

neodat.sdev = c(2.610,  2.660,  2.590,  1.940,  2.030,  2.050)
Dmat = diag(neodat.sdev)
covmat = Dmat%*%Cmat%*%Dmat

# assign row and column names to the covariance matrix
rownames(covmat) = c("x1","x2","x3","x4","x5","x6")
colnames(covmat)= rownames(covmat)

############ Congeneric Model #############
#model set up - this is the model presented in Figure 7.1, p.238 
new.model <- 'auditory =~ x1 +  x2 + x3   
visual =~ x4 + x5 + x6'            


#model fit
fit <- cfa(new.model, sample.cov=covmat, sample.nobs=200)

#results same as in Table 7.1, page 241
summary(fit,standardized=TRUE, fit.measures=TRUE, rsq=TRUE, modindices=TRUE)


############# Tau-equivalent Model - Auditory Memory only ############
#model set up - this is the model presented in Figure 7.1, p.238 
new.model <- 'auditory =~ v1*x1 +  v1*x2 + v1*x3
visual =~ x4 + x5 + x6'            

#model fit
fit <- cfa(new.model, sample.cov=covmat, sample.nobs=200)

#results same as in Table 7.3, page 246 (only model fit)
summary(fit,standardized=TRUE, fit.measures=TRUE, rsq=TRUE, modindices=TRUE)


############# Tau-equivalent Model - Auditory & Visual memory ############
#model set up - this is the model presented in Figure 7.1, p.238 
new.model <- 'auditory =~ v1*x1 +  v1*x2 + v1*x3
visual =~ v2*x4 + v2*x5 + v2*x6'

#model fit
fit <- cfa(new.model, sample.cov=covmat, sample.nobs=200)

#results same as in Table 7.3, page 246 (only model fit)
summary(fit,standardized=TRUE, fit.measures=TRUE, rsq=TRUE, modindices=TRUE)


############# Parallel Model - Auditory Memory only ############
#model set up - this is the model presented in Figure 7.1, p.238 
new.model <- 'auditory =~ v1*x1 +  v1*x2 + v1*x3
visual =~ v2*x4 + v2*x5 + v2*x6
#equality constraints of error variances (Auditory)
f1*x1 ~~ f1*x1
f1*x2 ~~ f1*x2
f1*x3 ~~ f1*x3'

#model fit
fit <- cfa(new.model, sample.cov=covmat, sample.nobs=200)

#results same as in Table 7.3, page 246 (only model fit)
summary(fit,standardized=TRUE, fit.measures=TRUE, rsq=TRUE, modindices=TRUE)


############# Parallel Model - Auditory & Visual memory ############
#model set up - this is the model presented in Figure 7.1, p.238 
new.model <- 'auditory =~ v1*x1 +  v1*x2 + v1*x3
visual =~ v2*x4 + v2*x5 + v2*x6
#equality constraints of error variances (Auditory & Visual)
f1*x1 ~~ f1*x1
f1*x2 ~~ f1*x2
f1*x3 ~~ f1*x3
f2*x4 ~~ f2*x4
f2*x5 ~~ f2*x5
f2*x6 ~~ f2*x6'

#model fit
fit <- cfa(new.model, sample.cov=covmat, sample.nobs=200)

#results same as in Table 7.4, page 248
summary(fit,standardized=TRUE, fit.measures=TRUE, rsq=TRUE, modindices=TRUE)

#to see the internal representation of the model (this helps understand the model definition)
inspect(fit,what="list")

########################################
# Longitudinal Invariance p.252 - 268  #
########################################

#Input matrix: the correlation matrix of Figure 7.2
cormat <- '
1.000
0.736  1.000
0.731  0.648  1.000
0.771  0.694  0.700  1.000
0.685  0.512  0.496  0.508  1.000
0.481  0.638  0.431  0.449  0.726  1.000
0.485  0.442  0.635  0.456  0.743  0.672  1.000
0.508  0.469  0.453  0.627  0.759  0.689  0.695  1.000'

Cmat <- getCov(cormat)

neodat.sdev = c(1.940,  2.030,  2.050,  1.990,  2.610,  2.660, 2.590,  2.550)
Dmat = diag(neodat.sdev)
covmat = Dmat%*%Cmat%*%Dmat

# assign row and column names to the covariance matrix
rownames(covmat) = c("A1","B1","C1","D1","A2","B2","C2","D2")
colnames(covmat)= rownames(covmat)

#Means: 1.500  1.320  1.450  1.410  6.600  6.420  6.560  6.310
#I couldn't find a way to bring in the indicator means, 
#so I used the meanstructure = TRUE argument in cfa call (which is not the same thing). 

##################################################################
#Results of the Equal Form Longitudinal Model of Job Satisfaction
##################################################################
#model set up - this is the model presented in Figure 7.2, p.254 
new.model <- 'SATIS1 =~ A1 + B1 + C1 + D1
SATIS2 =~ A2 + B2 + C2 + D2
A1 ~~ A2
B1 ~~ B2
C1 ~~ C2
D1 ~~ D2'

#model fit
fit <- cfa(new.model, sample.cov=covmat, sample.nobs=250, meanstructure = TRUE)

#results same as in Table 7.6, page 260 and fit same as in Table 7.7
summary(fit,standardized=TRUE, fit.measures=TRUE, rsq=TRUE, modindices=TRUE)

##############################################################################
#Results of the Equal Factor Loadings Longitudinal Model of Job Satisfaction #
##############################################################################
#model set up - this is the model presented in Figure 7.2, p.254 
# loadings constrained to be equal using custom labels
new.model <- 'SATIS1 =~ v1*A1 +  v2*B1 + v3*C1 + v4*D1
SATIS2 =~ v1*A2 + v2*B2 + v3*C2 + v4*D2
A1 ~~ A2
B1 ~~ B2
C1 ~~ C2
D1 ~~ D2'

#model fit
fit <- cfa(new.model, sample.cov=covmat, sample.nobs=250,meanstructure = TRUE)

#Fit same as in Table 7.7
summary(fit,standardized=TRUE, fit.measures=TRUE, rsq=TRUE, modindices=TRUE)


##############################################################################
#Equal indicator intercepts Longitudinal Model of Job Satisfaction #
##############################################################################
#model set up - this is the model presented in Figure 7.2, p.254 
# loadings constrained to be equal using custom labels
new.model <- 'SATIS1 =~ v1*A1 +  v2*B1 + v3*C1 + v4*D1
SATIS2 =~ v1*A2 + v2*B2 + v3*C2 + v4*D2
A1 ~~ A2
B1 ~~ B2
C1 ~~ C2
D1 ~~ D2
# intercepts constrained to be equal
# using custom labels
A1 ~ int2 * 1
A2 ~ int2 * 1
B1 ~ int3 * 1
B2 ~ int3 * 1
C1 ~ int4 * 1
C2 ~ int4 * 1
D1 ~ int5 * 1
D2 ~ int5 * 1'

#model fit
fit <- cfa(new.model, sample.cov=covmat, sample.nobs=250,meanstructure=TRUE)

#results are not the same as in Table 7.7 - no way to bring in means in lavaan (?)
summary(fit,standardized=TRUE, fit.measures=TRUE, rsq=TRUE, modindices=TRUE)


##############################################################################
#Equal indicator error variances Longitudinal Model of Job Satisfaction #
##############################################################################
#model set up - this is the model presented in Figure 7.2, p.254
new.model <- 'SATIS1 =~ v1*A1 +  v2*B1 + v3*C1 + v4*D1
SATIS2 =~ v1*A2 + v2*B2 + v3*C2 + v4*D2
A1 ~~ A2
B1 ~~ B2
C1 ~~ C2
D1 ~~ D2
# intercepts constrained to be equal
# using custom labels
A1 ~ int2 * 1
A2 ~ int2 * 1
B1 ~ int3 * 1
B2 ~ int3 * 1
C1 ~ int4 * 1
C2 ~ int4 * 1
D1 ~ int5 * 1
D2 ~ int5 * 1
#error variances constrained to be equal
#using custom labels
f1*A1 ~~ f1*A1
f1*A2 ~~ f1*A2
f2*B1 ~~ f2*B1
f2*B2 ~~ f2*B2
f3*C1 ~~ f3*C1
f3*C2 ~~ f3*C2
f4*D1 ~~ f4*D1
f4*D2 ~~ f4*D2'

#model fit
fit <- cfa(new.model, sample.cov=covmat, sample.nobs=250,meanstructure=TRUE)

#results are not the same as in Table 7.7 - no way to bring in means in lavaan (?)
summary(fit,standardized=TRUE, fit.measures=TRUE, rsq=TRUE, modindices=TRUE)

#to see the internal representation of the model (this helps understand the model definition)
inspect(fit,what="list")

#############################################
#### Multiple Group Analysis pp. 268-299 ####
#############################################

# read Depression example from the book's companion website
depres <- read.table("http://people.bu.edu/tabrown/Ch7/MDDALL.dat")

# assign variable names
names(depres) <- c("Sex","M1","M2","M3","M4","M5","M6","M7","M8","M9")
head(depres)

#model set up - this is the model presented in Figure 7.2, p.238 
new.model <- 'majdep =~ M1 + M2 + M3 + M4 + M5 + M6 + M7 + M8 + M9
M1 ~~ M2'

### Separate analysis of males and females to inspect model fit

#model fit for Female
fit <- cfa(new.model, data=depres[Sex==0,-1])
#results are not the same as in Table xxxx
summary(fit,standardized=TRUE, fit.measures=TRUE, rsq=TRUE, modindices=TRUE)

#model fit for Male
fit2 <- cfa(new.model, data=depres[Sex==1,-1])
#results are not the same as in Table xxxx
summary(fit2,standardized=TRUE, fit.measures=TRUE, rsq=TRUE, modindices=TRUE)

#Tests of Measurement Invariance and Population Heterogeneity

#Measurement Invariance - same model fit as in Table 7.9 
# model 1: configural invariance
fit1 <- cfa(new.model, data=depres, group="Sex")
# model 2: weak invariance
fit2 <- cfa(new.model, data=depres, group="Sex",
            group.equal="loadings")
# model 3: strong invariance
fit3 <- cfa(new.model, data=depres, group="Sex",
            group.equal=c("loadings", "intercepts"))
# model 4: equal loadings + intercepts + residuals
fit4 <- cfa(new.model, data=depres, group="Sex",
            group.equal=c("loadings", "intercepts", "residuals"))
#Population heterogeneity- # same model fit as in Table 7.9 
# model 5: equal loadings + intercepts + residuals + factor variances
fit5 <- cfa(new.model, data=depres, group="Sex",
            group.equal=c("loadings", "intercepts", "residuals","lv.variances"))
# model 6: equal loadings + intercepts + residuals + factor variances + factor means
fit6 <- cfa(new.model, data=depres, group="Sex",
            group.equal=c("loadings", "intercepts", "residuals","lv.variances","means"))

############## Alternatively for Models 1-3 and equal loadings + intercepts + means ##########
############## with chi-square tests of differences ##########
install.packages("semTools")
library(semTools)
measurementInvariance(new.model, data=depres, group="Sex", strict=FALSE)
###########################################

###################################
##### MIMIC Models pp.304-316 #####
###################################

cormat <- '
1.000
0.705   1.000
0.724   0.646   1.000
0.213   0.195   0.190   1.000
0.149   0.142   0.128   0.521   1.000
0.155   0.162   0.135   0.557   0.479   1.000'

Cmat <- getCov(cormat)

neodat.sdev = c(2.26, 2.73, 2.11, 2.32, 2.61, 2.44)
Dmat = diag(neodat.sdev)
covmat = Dmat%*%Cmat%*%Dmat

# assign row and column names to the covariance matrix
rownames(covmat) = c("S1","S2","S3","A1","A2","A3")
colnames(covmat)= rownames(covmat)

# Step 1. Ensure that the two-factor model of Social Phobia
# and Agoraphobia is reasonable and good fitting in the full sample (N = 730).

#model set up - this is the CFA model of the correlation matrix presented in Figure 7.5, p.308 
new.model <- 'socialphobia =~ S1 + S2 + S3   
agoraphobia =~ A1 + A2 + A3
socialphobia ~~ agoraphobia'            

#model fit
fit <- cfa(new.model, sample.cov=covmat, sample.nobs=730)

#results same as p. 308
summary(fit,standardized=TRUE, fit.measures=TRUE, rsq=TRUE, modindices=TRUE)


# Step 2. MIMIC Model of Fig. 7.5, p. 308.

#Introduce sex
cormat <- '
1.000
0.705   1.000
0.724   0.646   1.000
0.213   0.195   0.190   1.000
0.149   0.142   0.128   0.521   1.000
0.155   0.162   0.135   0.557   0.479   1.000
-0.019  -0.024  -0.029  -0.110  -0.074  -0.291   1.000'


Cmat <- getCov(cormat)

neodat.sdev = c(2.26, 2.73, 2.11, 2.32, 2.61, 2.44, 0.50)
Dmat = diag(neodat.sdev)
covmat = Dmat%*%Cmat%*%Dmat

# assign row and column names to the covariance matrix
rownames(covmat) = c("S1","S2","S3","A1","A2","A3","sex")
colnames(covmat)= rownames(covmat)

new.model <- 'socialphobia =~ S1 + S2 + S3   
agoraphobia =~ A1 + A2 + A3
socialphobia ~~ agoraphobia
#regressions
socialphobia ~ sex
agoraphobia ~ sex'            

#model fit
fit <- cfa(new.model, sample.cov=covmat, sample.nobs=730)

#results same as p. 315 - TABLE 7.15. Results of MIMIC Model of Social Phobia and Agoraphobia
summary(fit,standardized=TRUE, fit.measures=TRUE, rsq=TRUE, modindices=TRUE)


### Extra: Evaluation of measurement invariance by fixing all direct effects 
### between the covariate and the indicators to zero and then inspecting modification indices.
### pp.314-316

new.model <- 'socialphobia =~ S1 + S2 + S3   
agoraphobia =~ A1 + A2 + A3
socialphobia ~~ agoraphobia
#regressions
socialphobia ~ sex
agoraphobia ~ sex
#fix direct effects between sex and the indicators to zero
A3 ~ 0*sex
A2 ~ 0*sex
A1 ~ 0*sex'            

#model fit
fit <- cfa(new.model, sample.cov=covmat, sample.nobs=730)

#results same as p. 315 - TABLE 7.15. Results of MIMIC Model of Social Phobia and Agoraphobia
summary(fit,standardized=TRUE, fit.measures=TRUE, rsq=TRUE, modindices=TRUE)