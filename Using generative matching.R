# Chen, Pan, and Xu "Sources of Authoritarian Responsiveness: A Field Experiment in China"
# Replication File -- Table 3 Column (5)

#Data should be downloaded here: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/UMIBSL
#Research paper can be read here: https://onlinelibrary.wiley.com/doi/full/10.1111/ajps.12207

library(foreign)
library (randomForest)
library(rbounds)
library(Matching)
library(rgenoud)
library (rbounds)

# Compliance group
x0<-subset(x,treat==0 & posted == 1)
x1<-subset(x,treat==1 & posted == 1)
x2<-subset(x,treat==2 & posted == 1)
x3<-subset(x,treat==3 & posted == 1)
l0<-dim(x0)[1];l1<- dim(x1)[1];l2<-dim(x2)[1];l3<-dim(x3)[1]

#Change the NA into median
as.numeric(unlist(x0))
as.numeric(unlist(x1))
as.numeric(unlist(x2))
as.numeric(unlist(x3))

NA_change <-function(X){
  for (a in c(1:ncol(X))){
    if (class(X[,a])=='numeric'){
      im=median(X[,a],na.rm=TRUE)
      for (b in c(1:length(X[,a]))){
        if (is.na(X[b,a])==TRUE){
          X[b,a]=im
        }
      }
    }
  }
  return(X)
}

x0 = NA_change(x0)
x1 = NA_change(x1)
x2 = NA_change(x2)
x3 = NA_change(x3)

#Genetic matching between the first treatment and the control group
Data_1 = rbind (x0,x1) #Combine the treatment and control group into a big dataset
Data_1$treat = ifelse(Data_1$tr1 == 1, 1, ifelse(Data_1$treat == 1,1,0)) #If the tr1 == 0, then treatment == 0 
Data_1<-NA_change(Data_1) #Change the NA into a numeric number
glm1  <- glm(treat ~ logpop + logpop00 + popgrow + sexratio + illit + minor+ unemploy + wk_agr +
                wk_ind + wk_ser, income + logincome + loggdp + avggrowth +
                logoutput_agr + logoutput_ind + logoutput_ser + numlfirm + loginvest +
                logsaving + loggov_rev + loggov_exp, data = Data_1, family = binomial)
fit1 = glm1$fitted #Calculate the propensity score
cbind(Data_1,fit1) #Genetic matching with propensity score
Balance_Mat_1= cbind(Data_1$logpop, Data_1$logpop00, Data_1$popgrow, Data_1$sexratio, Data_1$illit, Data_1$minor, Data_1$unemploy, Data_1$wk_agr,
                     Data_1$wk_ind, Data_1$wk_ser, Data_1$income, Data_1$logincome, Data_1$loggdp, Data_1$avggrowth,
                     Data_1$logoutput_agr, Data_1$logoutput_ind + Data_1$logoutput_ser, Data_1$numlfirm, Data_1$loginvest,
                     Data_1$logsaving, Data_1$loggov_rev, Data_1$loggov_exp)
genout_1 <- GenMatch(Tr=Data_1$treat, X=Data_1, BalanceMatrix=Balance_Mat_1, estimand="ATT", M=1, replace = TRUE,
                       pop.size=100, max.generations=15, wait.generations=4)
mout_1 <- Match(Y = Data_1$response_public,Tr = Data_1$treat, X=Data_1, estimand="ATT", Weight.matrix=genout_1, M=1, replace = TRUE)
mout_1$est
mout_1$se
c(mout_1$est - 1.96*mout_1$se, mout_1$est + 1.96*mout_1$se)
mb_1 <- MatchBalance(Data_1$treat ~ Data_1$logpop + Data_1$logpop00 + Data_1$popgrow + Data_1$sexratio + Data_1$illit + Data_1$minor+ Data_1$unemploy + Data_1$wk_agr +
                       Data_1$wk_ind + Data_1$wk_ser + Data_1$income + Data_1$logincome + Data_1$loggdp + Data_1$avggrowth + Data_1$logoutput_agr + Data_1$logoutput_ind + Data_1$logoutput_ser + Data_1$numlfirm + Data_1$loginvest +
                       Data_1$logsaving + Data_1$loggov_rev + Data_1$loggov_exp,
                       match.out=mout_1, nboots=1000)

#Treatment effect: 0.03238095
#Confidence interval: -0.004181026  0.068942931
#Standard error: 0.01865407

#Genetic matching between the second treatment and the control group
Data_2 = rbind (x0,x2) #Combine the treatment and control group into a big dataset
Data_2$treat = ifelse(Data_2$tr2 == 1, 1, ifelse(Data_2$treat == 1,1,0)) #If the tr1 == 0, then treatment == 0 
Data_2<-NA_change(Data_1) #Change the NA into numeric number
glm2  <- glm(treat ~ logpop + logpop00 + popgrow + sexratio + illit + minor+ unemploy + wk_agr +
               wk_ind + wk_ser, income + logincome + loggdp + avggrowth +
               logoutput_agr + logoutput_ind + logoutput_ser + numlfirm + loginvest +
               logsaving + loggov_rev + loggov_exp, data = Data_2, family = binomial)
fit2 = glm2$fitted #Calculate the propensity score
cbind(Data_2,fit2) #Genetic matching with propensity score
Balance_Mat_2= cbind(Data_2$logpop, Data_2$logpop00, Data_2$popgrow, Data_2$sexratio, Data_2$illit, Data_2$minor, Data_2$unemploy, Data_2$wk_agr,
                     Data_2$wk_ind, Data_2$wk_ser, Data_2$income, Data_2$logincome, Data_2$loggdp, Data_2$avggrowth,
                     Data_2$logoutput_agr, Data_2$logoutput_ind + Data_2$logoutput_ser, Data_2$numlfirm, Data_2$loginvest,
                     Data_2$logsaving, Data_2$loggov_rev, Data_2$loggov_exp)
genout_2 <- GenMatch(Tr=Data_2$treat, X=Data_2, BalanceMatrix=Balance_Mat_2, estimand="ATT", M=1, replace = TRUE,
                     pop.size=100, max.generations=20, wait.generations=4)
mout_2 <- Match(Y = Data_2$response_public,Tr = Data_2$treat, X=Data_2, estimand="ATT", Weight.matrix=genout_2, M=1, replace = TRUE)
mout_2$est
mout_2$se
c(mout_2$est - 1.96*mout_2$se, mout_2$est + 1.96*mout_2$se)
#Treatment effect: 0.06285714
#Standard error: 0.02404087
#Confidence interval: 0.01573703 0.10997726

#Genetic matching between the third treatment and control group
Data_3 = rbind (x0,x3) #Combine the treatment and control group into a big dataset
Data_3$treat = ifelse(Data_3$tr3 == 1, 1, ifelse(Data_3$treat == 1,1,0)) #If the tr1 == 0, then treatment == 0 
Data_3<-NA_change(Data_3) #Change the NA into numeric number
glm3  <- glm(treat ~ logpop + logpop00 + popgrow + sexratio + illit + minor+ unemploy + wk_agr +
               wk_ind + wk_ser, income + logincome + loggdp + avggrowth +
               logoutput_agr + logoutput_ind + logoutput_ser + numlfirm + loginvest +
               logsaving + loggov_rev + loggov_exp, data = Data_3, family = binomial)
fit3 = glm3$fitted #Calculate the propensity score
cbind(Data_3,fit3) #Genetic matching with propensity score
Balance_Mat_3= cbind(Data_3$logpop, Data_3$logpop00, Data_3$popgrow, Data_3$sexratio, Data_3$illit, Data_3$minor, Data_3$unemploy, Data_3$wk_agr,
                     Data_3$wk_ind, Data_3$wk_ser, Data_3$income, Data_3$logincome, Data_3$loggdp, Data_3$avggrowth,
                     Data_3$logoutput_agr, Data_3$logoutput_ind + Data_3$logoutput_ser, Data_3$numlfirm, Data_3$loginvest,
                     Data_3$logsaving, Data_3$loggov_rev, Data_3$loggov_exp)
genout_3 <- GenMatch(Tr=Data_3$treat, X=Data_3, BalanceMatrix=Balance_Mat_3, estimand="ATT", M=1, replace = TRUE,
                     pop.size=100, max.generations=20, wait.generations=4)
mout_3 <- Match(Y = Data_3$response_public,Tr = Data_3$treat, X=Data_3, estimand="ATT", Weight.matrix=genout_3, M=1, replace = TRUE)
mout_3$est
mout_3$se
c(mout_3$est - 1.96*mout_3$se, mout_3$est + 1.96*mout_3$se)
#Treatment effect: 0.001893939
#Standard error: 0.01014351
#Confidence interval: -0.01798734  0.02177522

# Calculate the gamma (Rosenbaum's method of sensitivity analyses for 500 observations)
psens(x1$response_public, x0$response_public, Gamma = 1.5, GammaInc = .1) #Gamma = 1.3, p-score = 0.036
psens(x2$response_public, x0$response_public, Gamma = 1.5, GammaInc = .1) #Gamma = 1.0, p-score = 0.03
psens(x3$response_public, x0$response_public, Gamma = 1.5, GammaInc = .1) #Gamma = 1.0, p-score = 0.06

