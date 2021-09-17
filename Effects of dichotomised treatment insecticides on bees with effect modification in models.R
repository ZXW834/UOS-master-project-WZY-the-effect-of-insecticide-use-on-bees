Data<-read.csv("D:/University of Southampton/Summer project/Data_Park et al revised1.csv",header = T)
library("geepack")

Data$orchard<-Data$ï..orchard
Data<-Data[,-1]
Outcomes<-cbind(Data[,5:11],Data[,5:11],Data[,5:11])
detach(Data)
attach(Data)
summary(Data)

#---------------------------Total IUI------------------------------------
Data$eiqB11iuid<-ifelse(eiqB11.ins<=median(eiqB11.ins), 0, 1)

barplot(table(Data$eiqB11iuid),ylim = c(0,40),xlab = "Dichotomised IUI",ylab="Frequency",main = "Histogram of Dichotomised IUI")
Data$regionnew<-ifelse(region=="LO" | region=="GV", 0, 1)
tapply(long,as.factor(region),mean)
tapply(lat,as.factor(region),mean)
# GV is closer to LO

tablepositivityregion<-table(Data$regionnew,Data$eiqB11iuid)
colnames(tablepositivityregion)<-c("IUI <= 28.7","IUI > 28.7")
rownames(tablepositivityregion)<-c("GV+LO","S")
prop.table(tablepositivityregion,1)

#      IUI <= 28.7 IUI > 28.7
#GV+LO       0.528      0.472
#S           0.538      0.462

summary(temp[Data$eiqB11iuid==0])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#14.2    18.3    20.9    20.8    23.2    27.8 

summary(temp[Data$eiqB11iuid==1])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#14.4    18.9    21.0    21.3    23.9    29.4 

summary(X2000nat[Data$eiqB11iuid==0])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.212   0.230   0.276   0.346   0.474   0.661 

summary(X2000nat[Data$eiqB11iuid==1])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.192   0.343   0.453   0.411   0.486   0.554 

dtreatmentcharacter<-data.frame(matrix(c(0.53,0.47,20.80,21.30,0.35,0.41),ncol = 2,byrow = T))
colnames(dtreatmentcharacter)<-c("IUI <= 28.7","IUI > 28.7")
rownames(dtreatmentcharacter)<-c("GV+LO(%)","Temp(°C)","% natural area(2km radius)")
dtreatmentcharacter

#                           IUI <= 28.7 IUI > 28.7
#GV+LO(%)                          0.53       0.47
#Temp(°C)                         20.80      21.30
#% natural area(2km radius)        0.35       0.41

#Imbalance covariates: region and temp and % natural area(2km radius)
#Effect modification: Region and temp, Confounding: % natural area

#Assuming exchangeability true
#------------------------------------------------------Stable IP weighting with effect modification under assumptions---------------------------------------------------
# estimation of denominator of ip weights continuous treatment
dIpwDicho <- glm(
  eiqB11iuid ~ X2000nat+I(X2000nat^2)+regionnew+temp+I(temp^2),
  data = Data,
  family = binomial()
)
summary(dIpwDicho)


pdDicho <- predict(dIpwDicho, type = "response")


# estimation of numerator of ip weights (no effect modification)
nIpwDicho <- glm(eiqB11iuid ~ 1, family = binomial(), data = Data)
summary(nIpwDicho)

pnDicho <- predict(nIpwDicho, type = "response")

Data$w <-
  ifelse(Data$eiqB11iuid == 0, (1 / (1 - pdDicho)),
         (1 / pdDicho))

summary(Data$w)

Data$sw <-
  ifelse(Data$eiqB11iuid == 0, ((1 - pnDicho) / (1 - pdDicho)),
         (pnDicho / pdDicho))

summary(Data$sw)

# estimation of numerator of ip weights sw_v
nIpwDichoeffectmodification <-
  glm(eiqB11iuid ~ as.factor(regionnew)+temp+I(temp^2), family = binomial(), data = Data)
summary(nIpwDichoeffectmodification)

pnDichoeffectmodification <- predict(nIpwDichoeffectmodification, type = "response")

Data$sw.a <-
  ifelse(Data$eiqB11iuid == 0, ((1 - pnDichoeffectmodification) / (1 - pdDicho)),
         (pnDichoeffectmodification / pdDicho))

summary(Data$sw.a)
weights.iuid<-cbind(w.iuid=Data$w,sw.iuid=Data$sw,swa.iuid=Data$sw.a)

Data.iuid<-list()
j.iuid<-1
for (i in 1:21){
  Data.iuid[[i]] <- geeglm(
  Outcomes[,i] ~ eiqB11iuid + as.factor(regionnew)+temp+I(temp^2)+eiqB11iuid:temp
  + eiqB11iuid:as.factor(regionnew),
  data = Data,
  weights = weights.iuid[,j.iuid],
  id = as.factor(orchard),
  corstr = "independence"
  )
  if (i %% 7 ==0){
    j.iuid<-j.iuid+1
  }
}
beta.iuid<-list()
SE.iuid<-list()
lcl.iuid<-list()
ucl.iuid<-list()
result.iuid<-list()
for (z in 1:21){
beta.iuid[[z]] <- coef(Data.iuid[[z]])
SE.iuid[[z]] <- coef(summary(Data.iuid[[z]]))[, 2]
lcl.iuid[[z]] <- beta.iuid[[z]] - qnorm(0.975) * SE.iuid[[z]]
ucl.iuid[[z]] <- beta.iuid[[z]] + qnorm(0.975) * SE.iuid[[z]]
result.iuid[[z]]<-cbind(beta.iuid[[z]], lcl.iuid[[z]], ucl.iuid[[z]])
colnames(result.iuid[[z]])<-c("beta", "lcl", "ucl")
}

#---------------------------Pre-bloom IUI------------------------------------
Data$eiqB11.iuipred<-ifelse(eiqB11I.pre<=median(eiqB11I.pre), 0, 1)
hist(Data$eiqB11.iuipred)

barplot(table(Data$eiqB11.iuipred),ylim = c(0,50),xlab = "Dichotomised pre-bloom IUI",ylab="Frequency",main = "Histogram of Dichotomised pre-bloom IUI")
Data$regionnew<-ifelse(region=="LO" | region=="GV", 0, 1)
tapply(long,as.factor(region),mean)
tapply(lat,as.factor(region),mean)
# GV is closer to LO
tablepositivityregion.iuipre<-table(Data$regionnew,Data$eiqB11.iuipred)
colnames(tablepositivityregion.iuipre)<-c("Total PUI <= 20.8","Total PUI > 20.8")
rownames(tablepositivityregion.iuipre)<-c("GV+LO","S")
prop.table(tablepositivityregion.iuipre,1)

summary(temp[Data$eiqB11.iuipred==0])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#14.4    18.7    22.0    21.2    23.9    29.4 

summary(temp[Data$eiqB11.iuipred==1])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#14.2    19.3    20.0    20.6    22.0    29.2 

summary(X2000nat[Data$eiqB11.iuipred==0])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.212   0.321   0.468   0.408   0.489   0.661 

summary(X2000nat[Data$eiqB11.iuipred==1])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.192   0.218   0.231   0.292   0.339   0.453

dtreatmentcharacter.iuipre<-data.frame(matrix(c(0.61,0.39,21.20,20.60,0.41,0.29),ncol = 2,byrow = T))
colnames(dtreatmentcharacter.iuipre)<-c("preBloom IUI <= 9.36","preBloom IUI > 9.36")
rownames(dtreatmentcharacter.iuipre)<-c("GV+LO(%)","Temp(°C)","% natural area(2km radius)")
dtreatmentcharacter.iuipre

#                           preBloom IUI <= 9.36 preBloom IUI > 9.36
#GV+LO(%)                                   0.61                0.39
#Temp(°C)                                  21.20               20.60
#% natural area(2km radius)                 0.41                0.29

#Imbalance covariates: region and temp and % natural area(2km radius)
#Effect modification: Region and temp, Confounding: % natural area

#------------------------------------------------------Stable IP weighting with effect modification---------------------------------------------------
# estimation of denominator of ip weights continuous treatment(IUI.pre)
dIpwDicho.iuipre <- glm(
  eiqB11.iuipred ~ temp+I(temp^2)+regionnew+X2000nat+I(X2000nat^2),
  data = Data,
  family = binomial()
)
summary(dIpwDicho.iuipre)


pdDicho.iuipre <- predict(dIpwDicho.iuipre, type = "response")


# estimation of numerator of ip weights
nIpwDicho.iuipre <- glm(eiqB11.iuipred ~ 1, family = binomial(), data = Data)
summary(nIpwDicho.iuipre)

pnDicho.iuipre <- predict(nIpwDicho.iuipre, type = "response")

Data$w.iuipre <-
  ifelse(Data$eiqB11.iuipred == 0, (1  / (1 - pdDicho.iuipre)),
         (1 / pdDicho.iuipre))

summary(Data$w.iuipre)

Data$sw.iuipre <-
  ifelse(Data$eiqB11.iuipred == 0, ((1 - pnDicho.iuipre) / (1 - pdDicho.iuipre)),
         (pnDicho.iuipre / pdDicho.iuipre))

summary(Data$sw.iuipre)

nIpwDichoeffectmodification.iuipre <-
  glm(eiqB11.iuipred ~ as.factor(regionnew)+temp+I(temp^2), family = binomial(), data = Data)
summary(nIpwDichoeffectmodification.iuipre)

pnDichoeffectmodification.iuipre <- predict(nIpwDichoeffectmodification.iuipre, type = "response")

Data$sw.a.iuipre <-
  ifelse(Data$eiqB11.iuipred == 0, ((1 - pnDichoeffectmodification.iuipre) / (1 - pdDicho.iuipre)),
         (pnDichoeffectmodification.iuipre / pdDicho.iuipre))

summary(Data$sw.a.iuipre)

weights.iui.pre<-cbind(w.iui.pre=Data$w.iuipre,sw.iui.pre=Data$sw.iuipre,swa.iui.pre=Data$sw.a.iuipre)

Data.iui.pre<-list()
j.iui.pre<-1
for (i in 1:21){
  Data.iui.pre[[i]] <- geeglm(
    Outcomes[,i] ~ eiqB11.iuipred + as.factor(regionnew)+temp+I(temp^2)+temp:eiqB11.iuipred
    + eiqB11.iuipred:as.factor(regionnew),
    data = Data,
    weights = weights.iui.pre[,j.iui.pre],
    id = as.factor(orchard),
    corstr = "independence"
  )
  if (i %% 7 ==0){
    j.iui.pre<-j.iui.pre+1
  }
}
beta.iui.pre<-list()
SE.iui.pre<-list()
lcl.iui.pre<-list()
ucl.iui.pre<-list()
result.iui.pre<-list()
for (z in 1:21){
  beta.iui.pre[[z]] <- coef(Data.iui.pre[[z]])
  SE.iui.pre[[z]] <- coef(summary(Data.iui.pre[[z]]))[, 2]
  lcl.iui.pre[[z]] <- beta.iui.pre[[z]] - qnorm(0.975) * SE.iui.pre[[z]]
  ucl.iui.pre[[z]] <- beta.iui.pre[[z]] + qnorm(0.975) * SE.iui.pre[[z]]
  result.iui.pre[[z]]<-cbind(beta.iui.pre[[z]], lcl.iui.pre[[z]], ucl.iui.pre[[z]])
  colnames(result.iui.pre[[z]])<-c("beta", "lcl", "ucl")
}

#---------------------------Bloom IUI------------------------------------
Data$eiqB11.iuiblmd<-ifelse(eiqB11I.blm<=median(eiqB11I.blm), 0, 1)
hist(Data$eiqB11.iuiblmd)

barplot(table(Data$eiqB11.iuiblmd),ylim = c(0,40),xlab = "Dichotomised bloom IUI",ylab="Frequency",main = "Histogram of Dichotomised bloom IUI")
Data$regionnew<-ifelse(region=="LO" | region=="GV", 0, 1)
tapply(long,as.factor(region),mean)
tapply(lat,as.factor(region),mean)
# GV is closer to LO
tablepositivityregion.iuiblm<-table(Data$regionnew,Data$eiqB11.iuiblmd)
colnames(tablepositivityregion.iuiblm)<-c("Total PUI <= 20.8","Total PUI > 20.8")
rownames(tablepositivityregion.iuiblm)<-c("GV+LO","S")
prop.table(tablepositivityregion.iuiblm,1)

summary(temp[Data$eiqB11.iuiblmd==0])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#14.2    19.4    20.9    21.2    23.7    29.4 

summary(temp[Data$eiqB11.iuiblmd==1])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#14.4    17.7    21.5    20.8    23.3    29.2

summary(X2000nat[Data$eiqB11.iuiblmd==0])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.212   0.231   0.474   0.393   0.489   0.559

summary(X2000nat[Data$eiqB11.iuiblmd==1])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.192   0.264   0.332   0.357   0.453   0.661 

dtreatmentcharacter.iuiblm<-data.frame(matrix(c(0.69,0.31,21.20,20.80,0.39,0.36),ncol = 2,byrow = T))
colnames(dtreatmentcharacter.iuiblm)<-c("Bloom IUI <= 9.36","Bloom IUI > 9.36")
rownames(dtreatmentcharacter.iuiblm)<-c("GV+LO(%)","Temp(°C)","% natural area(2km radius)")
dtreatmentcharacter.iuiblm

#                           Bloom IUI <= 9.36 Bloom IUI > 9.36
#GV+LO(%)                                0.69             0.31
#Temp(°C)                               21.20            20.80
#% natural area(2km radius)              0.39             0.36

#Imbalance covariates: region and temp and % natural area(2km radius)
#Effect modification: Region and temp, Confounding: % natural area

#------------------------------------------------------Stable IP weighting with effect modification---------------------------------------------------
# estimation of denominator of ip weights continuous treatment(IUI.blm)
dIpwDicho.iuiblm <- glm(
  eiqB11.iuiblmd ~ temp+I(temp^2)+regionnew+X2000nat+I(X2000nat^2),
  data = Data,
  family = binomial()
)
summary(dIpwDicho.iuiblm)


pdDicho.iuiblm <- predict(dIpwDicho.iuiblm, type = "response")


# estimation of numerator of ip weights
nIpwDicho.iuiblm <- glm(eiqB11.iuiblmd ~ 1, family = binomial(), data = Data)
summary(nIpwDicho.iuiblm)

pnDicho.iuiblm <- predict(nIpwDicho.iuiblm, type = "response")

Data$w.iuiblm <-
  ifelse(Data$eiqB11.iuiblmd == 0, (1/ (1 - pdDicho.iuiblm)),
         (1 / pdDicho.iuiblm))

summary(Data$w.iuiblm)

Data$sw.iuiblm <-
  ifelse(Data$eiqB11.iuiblmd == 0, ((1 - pnDicho.iuiblm) / (1 - pdDicho.iuiblm)),
         (pnDicho.iuiblm / pdDicho.iuiblm))

summary(Data$sw.iuiblm)
nIpwDichoeffectmodification.iuiblm <-
  glm(eiqB11.iuiblmd ~ as.factor(regionnew)+temp+I(temp^2), family = binomial(), data = Data)
summary(nIpwDichoeffectmodification.iuiblm)

pnDichoeffectmodification.iuiblm <- predict(nIpwDichoeffectmodification.iuiblm, type = "response")

Data$sw.a.iuiblm <-
  ifelse(Data$eiqB11.iuiblmd == 0, ((1 - pnDichoeffectmodification.iuiblm) / (1 - pdDicho.iuiblm)),
         (pnDichoeffectmodification.iuiblm / pdDicho.iuiblm))

summary(Data$sw.a.iuiblm)

weights.iui.blm<-cbind(w.iui.blm=Data$w.iuiblm,sw.iui.blm=Data$sw.iuiblm,swa.iui.blm=Data$sw.a.iuiblm)

Data.iui.blm<-list()
j.iui.blm<-1
for (i in 1:21){
  Data.iui.blm[[i]] <- geeglm(
    Outcomes[,i] ~ eiqB11.iuiblmd + as.factor(regionnew)+temp+I(temp^2)+temp:eiqB11.iuiblmd
    + eiqB11.iuiblmd:as.factor(regionnew),
    data = Data,
    weights = weights.iui.blm[,j.iui.blm],
    id = as.factor(orchard),
    corstr = "independence"
  )
  if (i %% 7 ==0){
    j.iui.blm<-j.iui.blm+1
  }
}
beta.iui.blm<-list()
SE.iui.blm<-list()
lcl.iui.blm<-list()
ucl.iui.blm<-list()
result.iui.blm<-list()
for (z in 1:21){
  beta.iui.blm[[z]] <- coef(Data.iui.blm[[z]])
  SE.iui.blm[[z]] <- coef(summary(Data.iui.blm[[z]]))[, 2]
  lcl.iui.blm[[z]] <- beta.iui.blm[[z]] - qnorm(0.975) * SE.iui.blm[[z]]
  ucl.iui.blm[[z]] <- beta.iui.blm[[z]] + qnorm(0.975) * SE.iui.blm[[z]]
  result.iui.blm[[z]]<-cbind(beta.iui.blm[[z]], lcl.iui.blm[[z]], ucl.iui.blm[[z]])
  colnames(result.iui.blm[[z]])<-c("beta", "lcl", "ucl")
}

#---------------------------Post-bloom IUI------------------------------------
Data$eiqB11.iuiposd<-ifelse(eiqB11I.pos<=median(eiqB11I.pos), 0, 1)
hist(Data$eiqB11.iuiposd)

barplot(table(Data$eiqB11.iuiposd),ylim = c(0,40),xlab = "Dichotomised post-bloom IUI",ylab="Frequency",main = "Histogram of Dichotomised post-bloom IUI")
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

Data$w.iuipos <-
  ifelse(Data$eiqB11.iuiposd == 0, (1 / (1 - pdDicho.iuipos)),
         (1 / pdDicho.iuipos))
summary(Data$w.iuipos)

# estimation of numerator of ip weights
nIpwDicho.iuipos <- glm(eiqB11.iuiposd ~ 1, family = binomial(), data = Data)
summary(nIpwDicho.iuipos)

pnDicho.iuipos <- predict(nIpwDicho.iuipos, type = "response")

Data$sw.iuipos <-
  ifelse(Data$eiqB11.iuiposd == 0, ((1 - pnDicho.iuipos) / (1 - pdDicho.iuipos)),
         (pnDicho.iuipos / pdDicho.iuipos))
summary(Data$sw.iuipos)
nIpwDichoeffectmodification.iuipos <-
  glm(eiqB11.iuiposd ~ as.factor(regionnew)+temp+I(temp^2), family = binomial(), data = Data)
summary(nIpwDichoeffectmodification.iuipos)

pnDichoeffectmodification.iuipos <- predict(nIpwDichoeffectmodification.iuipos, type = "response")

Data$sw.a.iuipos <-
  ifelse(Data$eiqB11.iuiposd == 0, ((1 - pnDichoeffectmodification.iuipos) / (1 - pdDicho.iuipos)),
         (pnDichoeffectmodification.iuipos / pdDicho.iuipos))
summary(Data$sw.a.iuipos)

weights.iui.pos<-cbind(w.iui.pos=Data$w.iuipos,sw.iui.pos=Data$sw.iuipos,swa.iui.pos=Data$sw.a.iuipos)

Data.iui.pos<-list()
j.iui.pos<-1
for (i in 1:21){
  Data.iui.pos[[i]] <- geeglm(
    Outcomes[,i] ~ eiqB11.iuiposd + as.factor(regionnew)+temp+I(temp^2)+temp:eiqB11.iuiposd
    + eiqB11.iuiposd:as.factor(regionnew),
    data = Data,
    weights = weights.iui.pos[,j.iui.pos],
    id = as.factor(orchard),
    corstr = "independence"
  )
  if (i %% 7 ==0){
    j.iui.pos<-j.iui.pos+1
  }
}
beta.iui.pos<-list()
SE.iui.pos<-list()
lcl.iui.pos<-list()
ucl.iui.pos<-list()
result.iui.pos<-list()
for (z in 1:21){
  beta.iui.pos[[z]] <- coef(Data.iui.pos[[z]])
  SE.iui.pos[[z]] <- coef(summary(Data.iui.pos[[z]]))[, 2]
  lcl.iui.pos[[z]] <- beta.iui.pos[[z]] - qnorm(0.975) * SE.iui.pos[[z]]
  ucl.iui.pos[[z]] <- beta.iui.pos[[z]] + qnorm(0.975) * SE.iui.pos[[z]]
  result.iui.pos[[z]]<-cbind(beta.iui.pos[[z]], lcl.iui.pos[[z]], ucl.iui.pos[[z]])
  colnames(result.iui.pos[[z]])<-c("beta", "lcl", "ucl")
}
