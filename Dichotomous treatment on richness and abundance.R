Data<-read.csv("D:/University of Southampton/Summer project/Data_Park et al revised1.csv",header = T)
library("geepack")

Data$orchard<-Data$ï..orchard
Data<-Data[,-1]
detach(Data)
attach(Data)
summary(Data)

Data$eiqB11d<-ifelse(eiqB11<=median(eiqB11), 0, 1)
hist(Data$eiqB11d)
Data$regionnew<-ifelse(region=="LO" | region=="GV", 0, 1)
tapply(long,as.factor(region),mean)
tapply(lat,as.factor(region),mean)
# GV is closer to LO

tablepositivityregion<-table(Data$regionnew,Data$eiqB11d)
colnames(tablepositivityregion)<-c("Total PUI <= 224.65","Total PUI > 224.5")
rownames(tablepositivityregion)<-c("GV+LO","S")
prop.table(tablepositivityregion,1)

#      Total PUI <= 224.65 Total PUI > 224.5
#GV+LO               0.361             0.639
#S                   0.731             0.269
#Positivity True

summary(temp[Data$eiqB11d==0])
#Min.  1st Qu.  Median    Mean 3rd Qu.    Max. 
#14.44   17.73   22.11   21.09   23.89   29.17 

summary(temp[Data$eiqB11d==1])
#Min.  1st Qu.  Median    Mean 3rd Qu.    Max. 
#14.18   18.98   20.68   21.00   23.10   29.44 

summary(X2000nat[Data$eiqB11d==0])
#Min.  1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1915  0.2304  0.3410  0.3572  0.4692  0.6610 

summary(X2000nat[Data$eiqB11d==1])
#Min.  1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2184  0.3213  0.4532  0.3974  0.4885  0.5594 

dtreatmentcharacter<-data.frame(matrix(c(0.36,0.63,21.09,21.00,0.36,0.40),ncol = 2,byrow = T))
colnames(dtreatmentcharacter)<-c("Total PUI <= 224.65","Total PUI > 224.5")
rownames(dtreatmentcharacter)<-c("GV+LO(%)","Temp(°C)","% natural area(2km radius)")
dtreatmentcharacter

#                            Total PUI <= 224.65 Total PUI > 224.5
#GV+LO (%)                                  0.36              0.63
#Temp (°C)                                 21.09             21.00
#% natural area (2km radius)                0.36              0.40

#Imbalance covariates: region and % natural area. Temperature is balanced after dichotomising the treatment with the median as the threshold
#Effect modification: Region, Confounding: % natural area

#Assuming exchangeability true
#------------------------------------------------------Stable IP weighting with effect modification under assumptions---------------------------------------------------
# estimation of denominator of ip weights continuous treatment
dIpwDicho <- glm(
  eiqB11d ~ X2000nat+I(X2000nat^2)+regionnew,
  data = Data,
  family = binomial()
)
summary(dIpwDicho)


pdDicho <- predict(dIpwDicho, type = "response")


# estimation of numerator of ip weights (no effect modification)
nIpwDicho <- glm(eiqB11d ~ 1, family = binomial(), data = Data)
summary(nIpwDicho)

pnDicho <- predict(nIpwDicho, type = "response")

Data$sw <-
  ifelse(Data$eiqB11d == 0, ((1 - pnDicho) / (1 - pdDicho)),
         (pnDicho / pdDicho))

summary(Data$sw)
prop.table(xtabs(Data$sw ~ Data$regionnew + Data$eiqB11d),1)
#Balanced region

Data.sw <- geeglm(
  apisAb ~ eiqB11d,
  data = Data,
  weights = sw,
  id = as.factor(orchard),
  corstr = "independence"
)
summary(Data.sw)

# estimation of numerator of ip weights
nIpwDichoeffectmodification <-
  glm(eiqB11d ~ as.factor(regionnew), family = binomial(), data = Data)
summary(nIpwDichoeffectmodification)

pnDichoeffectmodification <- predict(nIpwDichoeffectmodification, type = "response")

Data$sw.a <-
  ifelse(Data$eiqB11d == 0, ((1 - pnDichoeffectmodification) / (1 - pdDicho)),
         (pnDichoeffectmodification / pdDicho))

summary(Data$sw.a)
Data.emm <- geeglm(
  apisAb ~ eiqB11d + as.factor(regionnew)
  + eiqB11d:as.factor(regionnew),
  data = Data,
  weights = sw.a,
  id = as.factor(orchard),
  corstr = "independence"
)
summary(Data.emm)

beta <- coef(Data.emm)
SE <- coef(summary(Data.emm))[, 2]
lcl <- beta - qnorm(0.975) * SE
ucl <- beta + qnorm(0.975) * SE
cbind(beta, lcl, ucl)

#No causal effect of dichotomous treatment total PUI on the outcome honeybee abundance in no effect modification model
#There is the causal effect of dichotomous treatment total PUI on the outcome honeybee abundance in effect modification model
#Under assumption, we estimated that Total PUI >= 244.5 in the region(LO+GV) increases honeybee abundance by 6.36 unit. If it is in region S, the effect is -6.06 unit.

#  Coefficients:
#                                Estimate Std.err Wald Pr(>|W|)    
#  (Intercept)                       2.97    0.37 64.2  1.1e-15 ***
#  eiqB11d                           6.36    1.44 19.6  9.5e-06 ***
#  as.factor(regionnew)1             6.50    1.71 14.5  0.00014 ***
#  eiqB11d:as.factor(regionnew)1   -12.42    2.21 31.7  1.8e-08 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Correlation structure = independence 
#Estimated Scale Parameters:
  
#  Estimate Std.err
#(Intercept)     23.4     6.1
#Number of clusters:   19  Maximum cluster size: 4 


#------------------------------------------------------Stratifying the treatment by blooming period (FUI.post blooming)---------------------------------------------
Data$eiqB11.funposd<-ifelse(eiqB11F.pos<=median(eiqB11F.pos), 0, 1)
hist(Data$eiqB11.funposd)

Data$regionnew<-ifelse(region=="LO" | region=="GV", 0, 1)
tapply(long,as.factor(region),mean)
tapply(lat,as.factor(region),mean)
# GV is closer to LO
tablepositivityregion.funpos<-table(Data$regionnew,Data$eiqB11.funposd)
colnames(tablepositivityregion.funpos)<-c("Total PUI <= 25.5","Total PUI > 25.5")
rownames(tablepositivityregion.funpos)<-c("GV+LO","S")
prop.table(tablepositivityregion.funpos,1)

summary(temp[Data$eiqB11.funposd==0])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#14.2    17.7    20.0    20.2    22.8    27.8

summary(temp[Data$eiqB11.funposd==1])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#14.4    19.5    21.1    21.9    24.0    29.4 

summary(X2000nat[Data$eiqB11.funposd==0])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.212   0.230   0.321   0.322   0.343   0.661 

summary(X2000nat[Data$eiqB11.funposd==1])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.192   0.376   0.474   0.431   0.489   0.559 

dtreatmentcharacter.funpos<-data.frame(matrix(c(0.39,0.61,20.20,21.90,0.32,0.43),ncol = 2,byrow = T))
colnames(dtreatmentcharacter.funpos)<-c("Total PUI <= 25.5","Total PUI > 25.5")
rownames(dtreatmentcharacter.funpos)<-c("GV+LO(%)","Temp(°C)","% natural area(2km radius)")
dtreatmentcharacter.funpos

#                         Total PUI <= 25.5 Total PUI > 25.5
#GV+LO(%)                               0.39            0.61
#Temp(°C)                              20.20           21.90
#% natural area(2km radius)             0.32            0.43

#Imbalance covariates: temp, region and % natural area.
#Effect modification: Region and temp, Confounding: % natural area

#------------------------------------------------------Stable IP weighting with effect modification---------------------------------------------------
# estimation of denominator of ip weights continuous treatment(FUI.pos)
dIpwDicho.funpos <- glm(
  eiqB11.funposd ~ temp+I(temp^2)+X2000nat+I(X2000nat^2)+regionnew,
  data = Data,
  family = binomial()
)
summary(dIpwDicho.funpos)


pdDicho.funpos <- predict(dIpwDicho.funpos, type = "response")


# estimation of numerator of ip weights
nIpwDicho.funpos <- glm(eiqB11.funposd ~ 1, family = binomial(), data = Data)
summary(nIpwDicho.funpos)

pnDicho.funpos <- predict(nIpwDicho.funpos, type = "response")

Data$sw.funpos <-
  ifelse(Data$eiqB11.funposd == 0, ((1 - pnDicho.funpos) / (1 - pdDicho.funpos)),
         (pnDicho.funpos / pdDicho.funpos))

summary(Data$sw.funpos)
prop.table(xtabs(Data$sw.funpos ~ Data$regionnew + Data$eiqB11.funposd),1)
#Still imbalanced adjusted region
#                 Data$eiqB11.funposd
#Data$regionnew        0     1
#                0 0.341 0.659
#                1 0.423 0.577


Data.sw.funpos <- geeglm(
  socialRichF ~ eiqB11.funposd,
  data = Data,
  weights = sw.funpos,
  id = as.factor(orchard),
  corstr = "independence"
)
summary(Data.sw.funpos)

nIpwDichoeffectmodification.funpos <-
  glm(eiqB11.funposd ~ as.factor(regionnew)+temp+I(temp^2), family = binomial(), data = Data)
summary(nIpwDichoeffectmodification.funpos)

pnDichoeffectmodification.funpos <- predict(nIpwDichoeffectmodification.funpos, type = "response")

Data$sw.a.funpos <-
  ifelse(Data$eiqB11.funposd == 0, ((1 - pnDichoeffectmodification.funpos) / (1 - pdDicho.funpos)),
         (pnDichoeffectmodification.funpos / pdDicho.funpos))

summary(Data$sw.a.funpos)
Data.emm.funpos <- geeglm(
  socialRichF ~ eiqB11.funposd + as.factor(regionnew)+temp+I(temp^2)+temp:eiqB11.funposd
  + eiqB11.funposd:as.factor(regionnew),
  data = Data,
  weights = sw.a.funpos,
  id = as.factor(orchard),
  corstr = "independence"
)
summary(Data.emm.funpos)

#No causal effect of dichotomous treatment Post blooming FUI on the outcome social bee richness in no effect modification model
#There is the causal effect of dichotomous treatment Post blooming FUI on the outcome social bee richness in effect modification model, after adding covariates temp and region
#Under assumption, we estimated that Post blooming FUI >= 25.5 in the region(LO+GV) increases social bee richness by 2.64 unit at the same temperature. 1 degree increased in temperature decreases the effect by 0.112 units.  
#If it is in region S, the effect is 1.2 unit at the same temperature. 1 degree increased in temperature decreases the effect by 0.112 units.


#------------------------------------------------------Stratifying the treatment by blooming period (IUI.post blooming)---------------------------------------------
Data$eiqB11.iuiposd<-ifelse(eiqB11I.pos<=median(eiqB11I.pos), 0, 1)
hist(Data$eiqB11.iuiposd)

Data$regionnew<-ifelse(region=="LO" | region=="GV", 0, 1)
tapply(long,as.factor(region),mean)
tapply(lat,as.factor(region),mean)
# GV is closer to LO
tablepositivityregion.iuipos<-table(Data$regionnew,Data$eiqB11.iuiposd)
colnames(tablepositivityregion.iuipos)<-c("Total PUI <= 20.8","Total PUI > 20.8")
rownames(tablepositivityregion.iuipos)<-c("GV+LO","S")
prop.table(tablepositivityregion.iuipos,1)

summary(temp[Data$eiqB11.iuiposd==0])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#14.4    17.8    22.2    21.4    23.9    29.2 

summary(temp[Data$eiqB11.iuiposd==1])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#14.2    18.9    20.0    20.6    22.7    29.4 

summary(X2000nat[Data$eiqB11.iuiposd==0])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.192   0.230   0.376   0.373   0.474   0.661

summary(X2000nat[Data$eiqB11.iuiposd==1])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.218   0.321   0.343   0.381   0.486   0.554 

dtreatmentcharacter.iuipos<-data.frame(matrix(c(0.42,0.58,21.40,20.60,0.37,0.38),ncol = 2,byrow = T))
colnames(dtreatmentcharacter.iuipos)<-c("Total PUI <= 20.8","Total PUI > 20.8")
rownames(dtreatmentcharacter.iuipos)<-c("GV+LO(%)","Temp(°C)","% natural area(2km radius)")
dtreatmentcharacter.iuipos

#                         Total PUI <= 20.8 Total PUI > 20.8
#GV+LO(%)                               0.42            0.58
#Temp(°C)                              21.40           20.60
#% natural area(2km radius)             0.37            0.38

#Imbalance covariates: region and temp and % natural area(2km radius)
#Effect modification: Region and temp, Confounding: % natural area

#------------------------------------------------------Stable IP weighting with effect modification---------------------------------------------------
# estimation of denominator of ip weights continuous treatment(IUI.pos)
dIpwDicho.iuipos <- glm(
  eiqB11.iuiposd ~ temp+I(temp^2)+regionnew+X2000nat+I(X2000nat^2),
  data = Data,
  family = binomial()
)
summary(dIpwDicho.iuipos)


pdDicho.iuipos <- predict(dIpwDicho.iuipos, type = "response")


# estimation of numerator of ip weights
nIpwDicho.iuipos <- glm(eiqB11.iuiposd ~ 1, family = binomial(), data = Data)
summary(nIpwDicho.iuipos)

pnDicho.iuipos <- predict(nIpwDicho.iuipos, type = "response")

Data$sw.iuipos <-
  ifelse(Data$eiqB11.iuiposd == 0, ((1 - pnDicho.iuipos) / (1 - pdDicho.iuipos)),
         (pnDicho.iuipos / pdDicho.iuipos))

summary(Data$sw.iuipos)
prop.table(xtabs(Data$sw.iuipos ~ Data$regionnew + Data$eiqB11.iuiposd),1)

Data.sw.iuipos <- geeglm(
  socialAbF ~ eiqB11.iuiposd,
  data = Data,
  weights = sw.iuipos,
  id = as.factor(orchard),
  corstr = "independence"
)
summary(Data.sw.iuipos)

nIpwDichoeffectmodification.iuipos <-
  glm(eiqB11.iuiposd ~ as.factor(regionnew)+temp+I(temp^2), family = binomial(), data = Data)
summary(nIpwDichoeffectmodification.iuipos)

pnDichoeffectmodification.iuipos <- predict(nIpwDichoeffectmodification.iuipos, type = "response")

Data$sw.a.iuipos <-
  ifelse(Data$eiqB11.iuiposd == 0, ((1 - pnDichoeffectmodification.iuipos) / (1 - pdDicho.iuipos)),
         (pnDichoeffectmodification.iuipos / pdDicho.iuipos))

summary(Data$sw.a.iuipos)
Data.emm.iuipos <- geeglm(
  socialAbF ~ eiqB11.iuiposd + as.factor(regionnew)+temp+I(temp^2)+temp:eiqB11.iuiposd
  + eiqB11.iuiposd:as.factor(regionnew),
  data = Data,
  weights = sw.a.iuipos,
  id = as.factor(orchard),
  corstr = "independence"
)
summary(Data.emm.iuipos)

#No causal effect of dichotomous treatment Post blooming IUI on the outcome social bee abundance in no effect modification model
#There is the causal effect of dichotomous treatment Post blooming IUI on the outcome social bee abundance in effect modification model, after adding covariates temp and region
#Under assumption, we estimated that Post blooming IUI >= 20.8 decreases social bee abundance by 3.64 unit at the same temperature. 1 degree increased in temperature increases the effect by 0.173 units.  

