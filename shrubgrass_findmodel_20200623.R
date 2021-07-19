# N:\code\tkohler\R\shrubgrass_findmodel_20200623.R

# Written for the Southern Cal invasive grasses or SoCal_ShrubGrass work
# Experiment with various models (regsubsets, lm, random forest, & betareg), generate plots
# & explore relationships. Find best 2 or 3 variables for each model

# Note: this is a large script and should not be run from start to end. It was developed as needed for exploration of different models with 
# various sets of data, plotting, and predicting outputs from various combinations and models. It is meant to be run in sections, depending
# on which stage of the process is needed.
# Sections are noted with comments, and can be collapsed as needed in Rstudio. 
#
library(HH)
library(leaps) # regsubsets
library(ggplot2)
library(gghighlight)
library(GGally) #ggpairs
library(randomForest)
library(randomForestExplainer)
library(psych) # for corr.test
library(betareg)
library(car)
library(broom) #augment
library(ggRandomForests)
library(data.table)
library(stringr) # str_replace_all
library(Metrics)
library(lmtest)
library(raster)
library(cluster)
library(parallel)
library(rfUtilities)
library(mgsub) # string substitution
library(stringr) # reverse string selection


Sys.setenv(TMP = "D:\\p_data\\rworking")
rasterOptions(tmpdir="D:\\p_data\\rworking")

# ------------------- Set up tables with which to work -----------------------------
#setwd("N:\\project\\monitoring_vol\\SoCalShrubs\\statistics")
setwd("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions")

#incsv <- read.table("u1011_2009_2018_mar2novby3_0103.csv", header = TRUE, sep=",") # GEE data
incsv <- read.table("N:/project/monitoring_vol/SoCalShrubs/points/u1011_2009_2018_mar2novby3_20201104.csv", header = TRUE, sep=",")



colnames(incsv)

# Including the ts value for that year ie ts2009
y09 <- incsv[c(2:77,453,459,463,1)] 
y11 <- incsv[c(2,78:152,453,458,464,1)]
y13 <- incsv[c(2,153:227,453,457,465,1)]
y15 <- incsv[c(2,228:302,453,456,466,1)]
y17 <- incsv[c(2,303:377,453,455,467,1)]
y18 <- incsv[c(2,378:454,468,1)]

colnames(y17)

# # Take subset for recs with no fires for up to 5 years before date
# y09_nofires <- y09[which(y09$fireyr != 2005 & y09$fireyr != 2006 & y09$fireyr != 2007 & y09$fireyr != 2008 & y09$fireyr != 2009),]
# y11_nofires <- y11[which(y11$fireyr != 2007 & y11$fireyr != 2008 & y11$fireyr != 2009 & y11$fireyr != 2010 & y11$fireyr != 2011),]
# y13_nofires <- y13[which(y13$fireyr != 2009 & y13$fireyr != 2010 & y13$fireyr != 2011 & y13$fireyr != 2012 & y13$fireyr != 2013),]
# y15_nofires <- y15[which(y15$fireyr != 2011 & y15$fireyr != 2012 & y15$fireyr != 2013 & y15$fireyr != 2014 & y15$fireyr != 2015),]
# y17_nofires <- y17[which(y17$fireyr != 2013 & y17$fireyr != 2014 & y17$fireyr != 2015 & y17$fireyr != 2016 & y17$fireyr != 2017),]
# y18_nofires <- y18[which(y18$fireyr != 2014 & y18$fireyr != 2015 & y18$fireyr != 2016 & y18$fireyr != 2017 & y18$fireyr != 2018),]
# 
# colnames(y09_nofires)
# 
# # Create dataset for each year wirh no fires at all in history, so fireyr == 0
# y09_nofiresEver <- y09[which(y09$fireyr == 0),]
# 
# 
# # Faster way to create string for lm 
# s09str <- toString(colnames(y09_nofires[1:75]))
# s11str <- toString(colnames(y11_nofires[1:75]))
# s13str <- toString(colnames(y13_nofires[1:75]))
# s15str <- toString(colnames(y15_nofires[1:75]))
# s17str <- toString(colnames(y17_nofires[1:75]))
# s18str <- toString(colnames(y18_nofires[1:75]))
# s11str
# 
# str2useinmodel <- str_replace_all(s09str,", "," + ")
# str2useinmodel
# 
# y09_nofires <- y09_nofires[complete.cases(y09_nofires),]
# y11_nofires <- y11_nofires[complete.cases(y11_nofires),]
# y13_nofires <- y13_nofires[complete.cases(y13_nofires),]
# y15_nofires <- y15_nofires[complete.cases(y15_nofires),]
# y17_nofires <- y17_nofires[complete.cases(y17_nofires),]
# y18_nofires <- y18_nofires[complete.cases(y18_nofires),]

# Get complete cases to remove where herb cover is not known (or is NA)
y09 <- y09[complete.cases(y09),]
y11 <- y11[complete.cases(y11),]
y13 <- y13[complete.cases(y13),]
y15 <- y15[complete.cases(y15),]
y17 <- y17[complete.cases(y17),]
y18 <- y18[complete.cases(y18),]


# ------------------- 2018 data with edart-derived products ----------------------
setwd("N:\\project\\monitoring_vol\\SoCalShrubs\\statistics")

incsv <- read.table("N:\\project\\monitoring_vol\\SoCalShrubs\\shrubgrass_utm1011_fire_formodel_Michele.csv", header = TRUE, sep=",")
# Remove fires
nofires <- incsv[ which(incsv$fire1_0 ==0),] # currently set for 2014 & onwards

# Random Forest
rf <- randomForest(Herb_PCT ~ jja_RGA_mean + mam_NDVI_max + mam_RGA_sd,
                   importance=TRUE,na.action=na.exclude,type="regression",data=nofires)
rf
round(importance(rf),2)
bptest(rf, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test

rf_p <- predict(rf, nofires) #, norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
rf_pdf <- as.data.frame(rf_p); colnames(rf_pdf) <- "rfresponse"
aarf <- cbind(rf_pdf,nofires); paste("rmse = ", round(rmse(aarf$Herb_PCT,aarf$rfresponse),4))

# betareg for 2018 edart-data
nofires$hpsc <- nofires$Herb_PCT * .01
nofires$hpsc[nofires$hpsc == 1] <- .99
nofires$hpsc[nofires$hpsc == 0] <- .001
beta_model <- betareg(hpsc ~ jja_RGA_mean + mam_NDVI_max + mam_RGA_sd, na.action=na.exclude, link = "probit",data=nofires)
summary(beta_model)
bptest(beta_model, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test 
betareg_predict <- cbind(
  predict(beta_model, type = "response")
)
colnames(betareg_predict) <- c("brresponse")
aa <- cbind(betareg_predict,nofires)
paste("rmse = ", round(rmse(aa$hpsc,aa$brresponse),4))

# ------------------- Linear model section --------------------------------
lmodel <- lm(heb_2015 ~ ARGmax1503+ARGmea1503+ARGmed1503+ARGmin1503+ARGstd1503+B3max1503+B3mea1503+B3med1503+B3min1503+B3std1503+NBRmax1503+NBRmea1503+NBRmed1503+NBRmin1503+NBRstd1503+NDImax1503+NDImea1503+NDImed1503+NDImin1503+NDIstd1503+NDVmax1503+NDVmea1503+NDVmed1503+NDVmin1503+NDVstd1503+ARGmax1506+ARGmea1506+ARGmed1506+ARGmin1506+ARGstd1506+B3max1506+B3mea1506+B3med1506+B3min1506+B3std1506+NBRmax1506+NBRmea1506+NBRmed1506)

lmodel <- lm(heb_2015 ~ ARGmea1506 + NDVmax1503 + ARGstd1503)

summary(lmodel)
ggplot(y15_nofires, aes(x=heb_2015, y=ARGmea1506 + NDVmax1503 + ARGstd1503))+
  geom_smooth(method = lm, se = FALSE) +
  geom_point() + geom_text(aes(label=pltnum),check_overlap = TRUE)

lmodel <- lm(heb_2015 ~ NDVstd1503 + B3max1506 + NDIstd1509)
dev.off()

#lmodel <- lm(heb_2015 ~ NDVstd1503 + ARGstd1503 + B3max1506)
# This lm based on regsubsets below
#lmodel <- lm(Herb_PCT ~  son_RGA_max + mam_RGA_sd)# + mam_NDVI_max + mam_RGA_sd)  
summary(lmodel)
vif(lmodel)
par(mfrow=c(2,2))
plot(Herb_PCT ~ jja_RGA_max )
model.diag.metrics <- augment(lmodel)
head(model.diag.metrics)


# ------------------- Regsubsets to get top 5 single,double,& triple variable combinations ---------
str2useinmodel

Subsets.1 = regsubsets(heb_2015 ~ ARGmax1503 + B3std1503 + ARGmax1506 + NDVstd1503 + ARGmed1506  + B3max1506    + B3mea1506   + B3med1506  + B3min1506 + ARGmea1503 + NDVmax1503 +  NBRmax1503 + NBRmed1503 + NDImea1503 +  ARGmed1509    + ARGmax1509  + ARGmea1509  + NDVmin1509 + NDImin1509 + NBRstd1509  + NDIstd1509 + NBRmea1509 +  B3med1509 + NBRmed1509 + NDImed1509, 
                       data=y15_nofires,nbest=5)


summary( lm.regsubsets(Subsets.1,4))
# create an empty matrix
explore = matrix(data=NA,nrow=20,ncol=1)
# Top 1
model.1 = summary( lm.regsubsets(Subsets.1, 1))
tmp = round(model.1$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2], sep=" ")
explore[1] = (paste("model.1", model.1$call[2], "    adj r^2 = ",round(model.1$adj.r.squared,digits=3), " P values = ",tmp.1,sep=" "))

model.2 = summary( lm.regsubsets(Subsets.1, 2))
tmp = round(model.2$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2], sep=" ")
explore[2] = (paste("model.2", model.2$call[2], "    adj r^2 = ",round(model.2$adj.r.squared,digits=3), " P values = ",tmp.1,sep=" "))

model.3 = summary( lm.regsubsets(Subsets.1, 3))
tmp = round(model.3$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2], sep=" ")
explore[3] = (paste("model.3", model.3$call[2], "    adj r^2 = ",round(model.3$adj.r.squared,digits=3), " P values = ",tmp.1,sep=" "))

model.4 = summary( lm.regsubsets(Subsets.1, 4))
tmp = round(model.4$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2], sep=" ")
explore[4] = (paste("model.4", model.4$call[2], "    adj r^2 = ",round(model.4$adj.r.squared,digits=3), " P values = ",tmp.1,sep=" "))

model.5 = summary( lm.regsubsets(Subsets.1, 5))
tmp = round(model.5$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2], sep=" ")
explore[5] = (paste("model.5", model.5$call[2], "    adj r^2 = ",round(model.5$adj.r.squared,digits=3), " P values = ",tmp.1,sep=" "))

# Top 2
model.6 = summary( lm.regsubsets(Subsets.1, 6))
tmp = round(model.6$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2], sep=" ")
explore[6] = (paste("model.6", model.6$call[2], "    adj r^2 = ",round(model.6$adj.r.squared,digits=3), " P values = ",tmp.1,sep=" "))

model.7 = summary( lm.regsubsets(Subsets.1, 7))
tmp = round(model.7$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2], sep=" ")
explore[7] = (paste("model.7", model.7$call[2], "    adj r^2 = ",round(model.7$adj.r.squared,digits=3), " P values = ",tmp.1,sep=" "))

model.8 = summary( lm.regsubsets(Subsets.1, 8))
tmp = round(model.8$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2], sep=" ")
explore[8] = (paste("model.8", model.8$call[2], "    adj r^2 = ",round(model.8$adj.r.squared,digits=3), " P values = ",tmp.1,sep=" "))

model.9 = summary( lm.regsubsets(Subsets.1, 9))
tmp = round(model.9$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2], sep=" ")
explore[9] = (paste("model.9", model.9$call[2], "    adj r^2 = ",round(model.9$adj.r.squared,digits=3), " P values = ",tmp.1,sep=" "))

model.10 = summary( lm.regsubsets(Subsets.1, 10))
tmp = round(model.9$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2], sep=" ")
explore[10] = (paste("model.10", model.10$call[2], "    adj r^2 = ",round(model.10$adj.r.squared,digits=3), " P values = ",tmp.1,sep=" "))

# Top 3
model.11 = summary( lm.regsubsets(Subsets.1, 11))
tmp = round(model.11$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2],tmp[3], sep=" ")
explore[11] = (paste("model.11", model.11$call[2], "    adj r^2 = ",round(model.11$adj.r.squared,digits=3), " P values = ",tmp.1,sep=" "))

model.12 = summary( lm.regsubsets(Subsets.1, 12))
tmp = round(model.12$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2],tmp[3], sep=" ")
explore[12] = (paste("model.12", model.12$call[2], "    adj r^2 = ",round(model.12$adj.r.squared,digits=3), " P values = ",tmp.1,sep=" "))

model.13 = summary( lm.regsubsets(Subsets.1, 13))
tmp = round(model.13$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2],tmp[3], sep=" ")
explore[13] = (paste("model.13", model.13$call[2], "    adj r^2 = ",round(model.13$adj.r.squared,digits=3), " P values = ",tmp.1,sep=" "))
model.14 = summary( lm.regsubsets(Subsets.1, 14))
tmp = round(model.14$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2],tmp[3], sep=" ")
explore[14] = (paste("model.14", model.14$call[2], "    adj r^2 = ",round(model.14$adj.r.squared,digits=3), " P values = ",tmp.1,sep=" "))
model.15 = summary( lm.regsubsets(Subsets.1, 15))
tmp = round(model.15$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2],tmp[3], sep=" ")
explore[15] = (paste("model.15", model.15$call[2], "   adj r^2 = ",round(model.15$adj.r.squared,digits=3), " P values = ",tmp.1,sep=" "))

# Top 4
model.16 = summary( lm.regsubsets(Subsets.1, 16))
tmp = round(model.16$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2],tmp[3],tmp[4], sep=" ")
explore[16] = (paste("model.16", model.16$call[2], "    adj r^2 = ",round(model.16$adj.r.squared,digits=3), " P values = ",tmp.1, sep=" "))

model.17 = summary( lm.regsubsets(Subsets.1, 17))
tmp = round(model.17$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2],tmp[3],tmp[4], sep=" ")
explore[17] = (paste("model.17", model.17$call[2], "    adj r^2 = ",round(model.17$adj.r.squared,digits=3), " P values = ",tmp.1,sep=" "))

model.18 = summary( lm.regsubsets(Subsets.1, 18))
tmp = round(model.18$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2],tmp[3],tmp[4], sep=" ")
explore[18] = (paste("model.18", model.18$call[2], "    adj r^2 = ",round(model.18$adj.r.squared,digits=3), " P values = ",tmp.1,sep=" "))

# model.19 = summary( lm.regsubsets(Subsets.1, 19))
# tmp = round(model.19$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2],tmp[3],tmp[4], sep=" ")
# explore[19] = (paste("model.19", model.19$call[2], "    adj r^2 = ",round(model.19$adj.r.squared,digits=3), " P values = ",tmp.1,sep=" "))
# 
# model.20 = summary( lm.regsubsets(Subsets.1, 20))
# tmp = round(model.20$coefficients[,4], digits=2); tmp.1 = paste(tmp[1], tmp[2],tmp[3],tmp[4], sep=" ")
# explore[20] = (paste("model.20", model.20$call[2], "   adj r^2 = ",round(model.20$adj.r.squared,digits=3), " P values = ",tmp.1,sep=" "))

print(explore)

lmodel <- lm(heb_2015 ~ NDVstd1503 + NBRmea1509 + B3med1509 + NDImed1509)
summary(lmodel)
vif(lmodel)
ggplot(y15_nofires, aes(x=heb_2015, y=NDVstd1503 + NDIstd1509 + B3med1509 + NDImed1509))+
  geom_smooth(method = lm, se = FALSE) +
  geom_point() + geom_text(aes(label=pltnum))



# ------------------- Plotting ggplot  ----------------------------

nofires2 <- nofires[complete.cases(nofires$jja_RGA_mean,nofires$mam_NDVI_max,nofires$mam_RGA_sd),]
nrow(nofires2)
detach(nofires2)
attach(y09_nofires)

c <- lm(heb_2015 ~ ARGmax1503 + ARGmax1506)
summary(c)
plot(c, which = 2)
model.diag.metrics <- augment(c)
head(model.diag.metrics)
coeff=coefficients(c)
y15_nofires$fits <- fitted.values(c)
eq = "ARGmax1503 + ARGmax1506"
plot.new()
pdf("NBRstd0903ARGmed0906B3max0909_y09_nofires.pdf") 
ggplot(y15_nofires, aes(x=heb_2015, y=fits)) + # + mam_NDVI_max + mam_RGA_sd,colour=Herb_PCT)) +
  geom_smooth(method = lm, se = FALSE) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  ggtitle(eq)
#geom_text(aes(label=NewID)) 
dev.off() 

# Plot graph in color for newsletter with 
rf <- randomForest(heb_2009 ~  ARGmax0903 + NBRstd0903 + ARGmed0906 + B3std0909, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)

rf_p <- predict(rf, y09_nofires) #, norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
rf_pdf <- as.data.frame(rf_p); colnames(rf_pdf) <- "rfresponse"
aarf <- cbind(rf_pdf,y09_nofires); paste("rmse = ", round(rmse(aarf$heb_2009,aarf$rfresponse),4))
#  plot betareg output w ggplot
setwd("N:\\project\\monitoring_vol\\SoCalShrubs\\figures")
pdf("heb_2009_best4_rf.pdf") 
ggplot(aarf, aes(x=heb_2009, y=rfresponse)) +
  geom_point() + #geom_text(aes(label=pltnum)) +
  scale_fill_grey() +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  #geom_abline(intercept=0,slope=1,color="red",size=1) +
  geom_segment(aes(x=0,y=0,xend=100,yend=100),color="red",size=1) +
  #geom_smooth(method = "lm", se = TRUE,colour = "red") +
  xlab("Observed Herbaceous Cover (%)") + 
  theme(axis.title.x = element_text(face="bold",color="blue",size=16)) +
  ylab("Predicted Herbaceous Cover (%)") +
  theme(axis.title.y = element_text(face="bold",color="blue",size=16)) 
#ggtitle("Predicted vs Observed Herbaceous Cover for 2018") +
#theme(plot.title = element_text(vjust = -22,hjust = .2,face="bold",color="blue"))
#gghighlight(pltnum %in% c(11188,11206),label_key = pltnum)
dev.off()


c2 <- lm(heb_2013 ~ ARGmax1303 + NDVmax1303 + B3med1309)
summary(c2)
plot(c2, which = 2)

ggplot(y11_nofires, aes(x=heb_2011, y=ARGmea1106+NDVmax1103+ARGstd1103))+
  geom_smooth(method = lm, se = FALSE) +
  geom_point()

predicted_df <- data.frame(herb_pct_pred = predict(c, nofires), hp=df$hp)

ggplot(model.diag.metrics, aes(heb_2009, NBRstd0903 + ARGmed0906 + B3max0909)) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  geom_segment(aes(xend = model.diag.metrics$.sigma, yend = model.diag.metrics$.fitted), color = "red", size = 0.3)
model.diag.metrics$.resid

ggplot(incsv, aes(x=heb_2009, y=ARGmax0903 ,colour=heb_2009)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point() +
  #gghighlight(NewID == 184,label_key = NewID) +
  gghighlight(NewID %in% c(450, 413, 447, 184, 116, 152),label_key = NewID)
#geom_abline(intercept = 0, slope = 0) 
#geom_smooth(method = "lm", se = FALSE)

model.diag.metrics <- augment(lmodel)
head(model.diag.metrics)

ggplot(model.diag.metrics, aes(youtube, sales)) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  geom_segment(aes(xend = youtube, yend = .fitted), color = "red", size = 0.3)

# -------  plot a - jja_RGA_mean
a <- lm(heb_2013 ~ ARG_mx2013MAM + NDVI_std2013MAM +      ARG_mx2013JJA + ARG_mn2013JJA + ARG_mn2013SON + B3_mx2013SON)
summary(a)
coeff=coefficients(a)
plot.new()
pdf("heb_2013_best4_rf.pdf") 
ggplot(y11_nofires, aes(x=heb_2011, y=ARGmed1106+NDVmax1103+ARGstd1103,colour=heb_2011)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)
dev.off() 

heb_2011 ~ ARGmed1106+NDVmax1103+ARGstd1103

# ------------------- Pairs panels -----------------

colnames(y09_nofires)

#pdf("pairs_13.pdf")
#plot.new()
pairs(y09_nofires[c(77,27,21,5)]) # 2018 vars in 2009
pairs(y09_nofires[c(77,1,15,28,60)]) # random forest 4 best ARGmax0903 + NBRstd0903 + ARGmed0906 + B3std0909
pairs(y09_nofires[c(77,1,15,28)]) # random forest 3 best ARGmax0903 + NBRstd0903 + ARGmed0906
pairs(y09_nofires[c(79,15,28,60)]) # betareg 3 best NBRstd0903 + ARGmed0906 + B3std0909

colnames(y11_nofires)
pairs(y11_nofires[c(77,1:10)])
pairs(y11_nofires[c(77,27,21,5)]) # 2018 vars in 2011 heb_2011 ~ ARGmea1106 + NDVmax1103 + ARGstd1103
pairs(y11_nofires[c(77,20,15,26,36)]) # rf 4 best heb_2011 ~ NDIstd1103 + NBRstd1103 + ARGmax1106 + B3max1106
pairs(y11_nofires[c(77,25,27,35)]) # rf 3 best heb_2011 ~ NDVstd1103 + ARGmea1106 + B3std1106
pairs(y11_nofires[c(79,12,19,33,18)]) # betareg 4 best NBRmea1103 + NDImin1103 + B3med1106 + NDImed1103
pairs(y11_nofires[c(79,25,27,35)]) # betareg 3 best hpsc ~ NDVstd1103 + ARGmea1106 + B3std1106

colnames(y13_nofires)
pairs(y13_nofires[c(77,27,21,5)]) # 2018 vars in 2013 heb_2011 ~ ARGmea1106 + NDVmax1103 + ARGstd1103
pairs(y13_nofires[c(77,52,25,56,27)]) # rf 4 best heb_2013 ~ ARGmea1309 + NDVstd1303 + B3max1309 +  ARGmea1306
pairs(y13_nofires[c(77,27,25,56)]) # rf 3 best heb_2013 ~ ARGmea1306 + NDVstd1303 + B3max1309

pairs(y13_nofires[c(79,25,56,27)]) # BR 3 best hpsc ~ NDVstd1303 + B3max1309 + ARGmea1306

colnames(y15_nofires)
pairs(y15_nofires[c(77,27,21,5)]) # 2018 vars in 2015 heb_2015~ ARGmea1506 + NDVmax1503 + ARGstd1503
pairs(y15_nofires[c(77,1,28,26,15)]) # rf 4 best heb_2015 ~ ARGmax1503 + ARGmed1506 + ARGmax1506 +NBRstd1503
pairs(y15_nofires[c(77,1,26,74)]) # rf 3 best heb_2015 ~ ARGmax1503 + ARGmax1506 + NDVmin1509
pairs(y15_nofires[c(79,13,19,21,57)]) # BR 4 best hpsc ~ NBRmed1503 + NDImin1503 + NDVmax1503 + B3mea1509
pairs(y15_nofires[c(79,19,21,57)]) # BR 3 best hpsc ~ NDImin1503 + NDVmax1503 + B3mea1509

colnames(y17_nofires)
pairs(y17_nofires[c(77,27,21,5)]) # 2018 vars in 2015 heb_2015~ ARGmea1706 + NDVmax1703 + ARGstd1703
pairs(y17_nofires[c(77,26,28,25,53)]) # rf 4 best ARGmax1706 + ARGmed1706 + NBRstd1703 + ARGmed1709
pairs(y17_nofires[c(77,26,25,53)]) # rf 3 best heb_2017 ~ ARGmax1706 + NBRstd1703 + ARGmed1709

colnames(y18_nofires)
pairs(y18_nofires[c(77,27,21,5)]) # 2018 vars in 2018 Herb_2018~ ARGmea1806 + NDVmax1803 + ARGstd1803
pairs(y18_nofires[c(77,26,51,15,10)]) # rf 4 best Herb_2018 ~ ARGmax1806 + ARGmax1809 +NBRstd1803 + B3std1803
pairs(y18_nofires[c(77,26,15,10)]) # rf 3 best Herb_2018 ~ ARGmax1806 + NBRstd1803 +      B3std1803
pairs(y18_nofires[c(79,26,25,15,10)]) # BR 4 best hpsc ~ ARGmax1806 + NDVstd1803 + NBRstd1803 + B3std1803


#dev.off()
#graphics.off()


# 2018 vars in pairs
# ARGstd1303, NDVmax1303, ARGmea1306, heb_2013
# 5, 21, 27, 77
# -----------------------------------------------------------------------------------------------------------------------#
# ------------------- Random forest -------------------------------------------------

rf <- randomForest(heb_2009 ~ ARGmax0903 + NBRstd0903 + ARGmed0906 +  B3std0909, 
                    data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
rf
bptest(rf, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test
#plot(gg_vimp(rf)); round(importance(rf),2)

y09_nofires$predict <- predict(rf)
y09_nofires$resids <- y09_nofires$predict - y09_nofires$heb_2009
# generate heatmap based on mean of resids
# add column re whether rec is inside historic fire & separate heatmaps based on whether rec is inside historic fire

rf_p <- predict(rf, y09_nofires) #, norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
rf_pdf <- as.data.frame(rf_p); colnames(rf_pdf) <- "rfresponse"
aarf <- cbind(rf_pdf,y09_nofires); paste("rmse = ", round(rmse(aarf$heb_2009,aarf$rfresponse),4))

#  plot betareg output w ggplot
ggplot(aarf, aes(x=heb_2009, y=rfresponse)) +
  geom_point() + #geom_text(aes(label=pltnum)) +
  scale_fill_grey() +
  #geom_abline(intercept=0,slope=1,color="red",size=1) +
  geom_segment(aes(x=0,y=0,xend=100,yend=100),color="red",size=1) +
  #geom_smooth(method = "lm", se = TRUE,colour = "red") +
  xlab("Observed Herbaceous Cover") + 
  theme(axis.title.x = element_text(face="bold",color="blue")) +
  ylab("Predicted Herbaceous Cover") +
  theme(axis.title.y = element_text(face="bold",color="blue")) 
  #ggtitle("Predicted vs Observed Herbaceous Cover for 2018") +
  #theme(plot.title = element_text(vjust = -22,hjust = .2,face="bold",color="blue"))
  #gghighlight(pltnum %in% c(11188,11206),label_key = pltnum)

partialPlot(rf, y09_nofires, NDIstd0903, "versicolor")

head(aarf)
detach(y18_nofires)
detach()
attach(y09_nofires)
ss <- lm(heb_2009 ~ ARGmed0906 )
summary(ss); vif(ss)
ggplot(y09_nofires, aes(x=heb_2009, y=ARGmed0906  ,colour=heb_2009)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point() + geom_text(aes(label=pltnum)) +
  gghighlight(pltnum == 10060,label_key = pltnum) 


# ------------------- RF Cross Validation ---------------------------------
# ------ Generate a (cross validation) table of statistics for each rf run ----------
setwd("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables")

# 2009 origmodelyr: heb_2009 ~ ARGmax0903 + NBRstd0903 + ARGmed0906 +  B3std0909", "heb_2009 ~ ARGmax0903 + NBRstd0903 + ARGmed0906"
# 2011 origmodelyr: "heb_2009 ~ NDVstd0903 + ARGmea0906 + B3std0906"
# 2013 origmodelyr: "heb_2009 ~ ARGmea0906 + NDVstd0903 + B3max0909" 
# 2015 origmodelyr: "heb_2009 ~ ARGmax0903 + ARGmax0906 + NDVmin0909"
# 2017 origmodelyr: "heb_2009 ~ ARGmax0906 + NBRstd0903 + ARGmed0909"
# 2018 origmodelyr: "heb_2009 ~ ARGmax0906 + NBRstd0903 +  B3std0903"

intable <- y09_nofires
listofstrings <- c("heb_2009 ~ ARGmax0903 + NBRstd0903 + ARGmed0906 +  B3std0909",
                   "heb_2009 ~ ARGmax0903 + NBRstd0903 + ARGmed0906",
                   "heb_2009 ~ NDVstd0903 + ARGmea0906 + B3std0906",
                   "heb_2009 ~ ARGmea0906 + NDVstd0903 + B3max0909",
                   "heb_2009 ~ ARGmax0903 + ARGmax0906 + NDVmin0909",
                   "heb_2009 ~ ARGmax0906 + NBRstd0903 + ARGmed0909",
                   "heb_2009 ~ ARGmax0906 + NBRstd0903 +  B3std0903")

intable <- y11_nofires
listofstrings <- c("heb_2011 ~ ARGmax1103 + NBRstd1103 + ARGmed1106 +  B3std1109",
                   "heb_2011 ~ ARGmax1103 + NBRstd1103 + ARGmed1106",
                   "heb_2011 ~ NDVstd1103 + ARGmea1106 + B3std1106",
                   "heb_2011 ~ ARGmea1106 + NDVstd1103 + B3max1109",
                   "heb_2011 ~ ARGmax1103 + ARGmax1106 + NDVmin1109",
                   "heb_2011 ~ ARGmax1106 + NBRstd1103 + ARGmed1109",
                   "heb_2011 ~ ARGmax1106 + NBRstd1103 +  B3std1103")

intable <- y13_nofires
listofstrings <- c("heb_2013 ~ ARGmax1303 + NBRstd1303 + ARGmed1306 +  B3std1309",
                   "heb_2013 ~ ARGmax1303 + NBRstd1303 + ARGmed1306",
                   "heb_2013 ~ NDVstd1303 + ARGmea1306 + B3std1306",
                   "heb_2013 ~ ARGmea1306 + NDVstd1303 + B3max1309",
                   "heb_2013 ~ ARGmax1303 + ARGmax1306 + NDVmin1309",
                   "heb_2013 ~ ARGmax1306 + NBRstd1303 + ARGmed1309",
                   "heb_2013 ~ ARGmax1306 + NBRstd1303 +  B3std1303")

intable <- y15_nofires
listofstrings <- c("heb_2015 ~ ARGmax1503 + NBRstd1503 + ARGmed1506 +  B3std1509",
                   "heb_2015 ~ ARGmax1503 + NBRstd1503 + ARGmed1506",
                   "heb_2015 ~ NDVstd1503 + ARGmea1506 + B3std1506",
                   "heb_2015 ~ ARGmea1506 + NDVstd1503 + B3max1509",
                   "heb_2015 ~ ARGmax1503 + ARGmax1506 + NDVmin1509",
                   "heb_2015 ~ ARGmax1506 + NBRstd1503 + ARGmed1509",
                   "heb_2015 ~ ARGmax1506 + NBRstd1503 +  B3std1503")

intable <- y17_nofires
listofstrings <- c("heb_2017 ~ ARGmax1703 + NBRstd1703 + ARGmed1706 +  B3std1709",
                   "heb_2017 ~ ARGmax1703 + NBRstd1703 + ARGmed1706",
                   "heb_2017 ~ NDVstd1703 + ARGmea1706 + B3std1706",
                   "heb_2017 ~ ARGmea1706 + NDVstd1703 + B3max1709",
                   "heb_2017 ~ ARGmax1703 + ARGmax1706 + NDVmin1709",
                   "heb_2017 ~ ARGmax1706 + NBRstd1703 + ARGmed1709",
                   "heb_2017 ~ ARGmax1706 + NBRstd1703 +  B3std1703") 

intable <- y18_nofires
listofstrings <- c("Herb_2018 ~ ARGmax1803 + NBRstd1803 + ARGmed1806 +  B3std1809",
                   "Herb_2018 ~ ARGmax1803 + NBRstd1803 + ARGmed1806",
                   "Herb_2018 ~ NDVstd1803 + ARGmea1806 + B3std1806",
                   "Herb_2018 ~ ARGmea1806 + NDVstd1803 + B3max1809",
                   "Herb_2018 ~ ARGmax1803 + ARGmax1806 + NDVmin1809",
                   "Herb_2018 ~ ARGmax1806 + NBRstd1803 + ARGmed1809",
                   "Herb_2018 ~ ARGmax1806 + NBRstd1803 +  B3std1803")

for (instring in listofstrings){
  setwd("N:/project/monitoring_vol/SoCalShrubs/EE_Extractions/cvtables")
  rf <- randomForest(as.formula(instring),data = intable, importance = TRUE, type = "regression", na.action = na.exclude)
  outstr <- mgsub::mgsub(instring, c("\\+","\\,", "\\~","\\ ","\\(","\\)"), c("_","","_","","","")); outstr
  cvtab <- data.frame()
  #set10 <- c(1,2,3,4,5)
  for (i in 1:5){
    cvout <- rf.crossValidation(rf, intable, p=i/10, n=100, ntree=501)
    cvout
    cvtab[i,1] <- cvout$fit.var.exp
    cvtab[i,2] <- cvout$fit.mse
    cvtab[i,3] <- median(cvout$y.rmse)
    cvtab[i,4] <- var(cvout$y.rmse)
    cvtab[i,5] <- median(cvout$y.mbe)
    cvtab[i,6] <- median(cvout$y.mae)
    cvtab[i,7] <- median(cvout$model.mse)
    colnames(cvtab) <- c("fit.var.exp","fit.mse","med.y.rmse","var.y.rmse","med.y.mbe","med.y.mae","med.model.mse")
    
    tname <- paste("rfcv_", outstr,"_iter100.csv",sep='')
    #write.csv(cvtab,tname,col.names=TRUE)
  }
}

#------- Append model results into one file ------------------------------
library(purrr)
library(readr)
library(mgsub)
# substrRight <- function(x, n){
#   substr(x, nchar(x)-n+1, nchar(x))
# }

list_of_files <- list.files(path = "N:/project/monitoring_vol/SoCalShrubs/EE_Extractions/cvtables5_100_3var",
                            pattern = "rfcv",
                            full.names = FALSE)
all_models_tempint <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read_csv(.x, col_types = cols(), col_names = TRUE), .id = "file_name")
all_models_tempint$year <- gsub("_", "",substring(as.character(all_models_tempint$file_name), 10, 14))
all_models_tempint$var1 <- gsub("_", "",substring(as.character(all_models_tempint$file_name), 15, 25))
all_models_tempint$var2 <- gsub("_", "",substring(as.character(all_models_tempint$file_name), 26, 36))
all_models_tempint$var3 <- gsub("_", "",substring(as.character(all_models_tempint$file_name), 37, 46))
all_models_tempint$iter <- gsub("\\.", "",str_sub(as.character(all_models_tempint$file_name), -7,-4)) # all 100 at this point

# 20200529 - Create BaseModel field with no data year info
all_models_tempint$BaseModel <- mgsub::mgsub(all_models_tempint$file_name,c("rfcv_Herb_2018_","rfcv_heb_2017_","rfcv_heb_2015_","rfcv_heb_2013_","rfcv_heb_2011_","rfcv_heb_2009_","_iter100.csv"),c("","","","","","",""))
all_models_tempint$BaseModel <- mgsub::mgsub(all_models_tempint$BaseModel, c("0903","1103","1303","1503","1703","1803"),c("YR03","YR03","YR03","YR03","YR03","YR03"))
all_models_tempint$BaseModel <- mgsub::mgsub(all_models_tempint$BaseModel, c("0906","1106","1306","1506","1706","1806"),c("YR06","YR06","YR06","YR06","YR06","YR06"))
all_models_tempint$BaseModel <- mgsub::mgsub(all_models_tempint$BaseModel, c("0909","1109","1309","1509","1709","1809"),c("YR09","YR09","YR09","YR09","YR09","YR09"))
unique(all_models_tempint$BaseModel)

colnames(all_models_tempint) <- c("file_name","p_withheld", "fitvarexp","fitmse","medyrmse","varyrmse","medymbe","medymae","medmodelmse","DataYear","var1","var2","var3","iter","BaseModel")
all_models_tempint$p_withheld <- all_models_tempint$p_withheld * 10
write.csv(all_models_tempint, file = "allrfcv_20200528.csv")
# A lookup table is used to populate the OrigModelYr column. 2 additional cols are added

# ------ Set up for graphing the combination table --------------------
require(reshape2)
library(reshape)
require(ggplot2)
require(scales)
require(plot3D)
require(akima) #Interpolation of Irregularly and Regularly Spaced Data
require(rgl) # 3D Visualization Using OpenGL
library(manipulate)
library(data.table)

setwd("N:/project/monitoring_vol/SoCalShrubs/EE_Extractions/cvtables5_100_3var")
all_models_2more <- read.table("allrfcv_20200528.csv", header=T,sep=',',stringsAsFactors = FALSE)
head(all_models_2more); tail(all_models_2more)
colnames(all_models_2more) <- c("id","file_name","IndVars","OrigModelYr","p_withheld","fitvarexp","fitmse","medyrmse","varyrmse","medymbe","medymae","medmodelmse","DataYear","var1","var2","var3","iter")
# if this wasn't done in all_models_tempint then do it here
all_models_2more$BaseModel <- mgsub::mgsub(all_models_2more$IndVars,c("Herb_2018_","heb_2017_","heb_2015_","heb_2013_","heb_2011_","heb_2009_"),c("","","","","",""))
all_models_2more$BaseModel <- mgsub::mgsub(all_models_2more$BaseModel, c("0903","1103","1303","1503","1703","1803"),c("YR03","YR03","YR03","YR03","YR03","YR03"))
all_models_2more$BaseModel <- mgsub::mgsub(all_models_2more$BaseModel, c("0906","1106","1306","1506","1706","1806"),c("YR06","YR06","YR06","YR06","YR06","YR06"))
all_models_2more$BaseModel <- mgsub::mgsub(all_models_2more$BaseModel, c("0909","1109","1309","1509","1709","1809"),c("YR09","YR09","YR09","YR09","YR09","YR09"))
all_models_2more$OrigModelYr_ch <- as.character(all_models_2more$OrigModelYr)
all_models_2more$DataYear_ch <- as.character(all_models_2more$DataYear)

uniquebm <- unique(all_models_2more$BaseModel)
uniquedatayr <- unique(all_models_2more$DataYear)
#all_models_melt <- melt(all_models_2more, id.var = c("IndVars","OrigModelYr","p_withheld"))
#head(all_models_melt)

# ------ Calculate ratio of medyrmse at percent withheld 10%:50% for each ModelYr-DataYear combo -------------------
all_models_2more$medyrmseratio <- 0
for (eachbm in uniquebm){
  print(eachbm)
  for (eachdy in uniquedatayr) {
    print(eachdy)
    fitvar_pw10 <- all_models_2more[all_models_2more$p_withheld == 10 & all_models_2more$BaseModel == eachbm & all_models_2more$DataYear == eachdy,"medyrmse"]
    print(fitvar_pw10)
    fitvar_pw50 <- all_models_2more[all_models_2more$p_withheld == 50 & all_models_2more$BaseModel == eachbm & all_models_2more$DataYear == eachdy,"medyrmse"]
    print(fitvar_pw50)
    # calc medyrmse ratio
    medyrmseratio <- fitvar_pw10/fitvar_pw50
    medyrmseratio_ch <- paste(medyrmseratio, ":1",sep='')
    print(medyrmseratio_ch)
    all_models_2more$medyrmseratio[all_models_2more$BaseModel == eachbm & all_models_2more$DataYear == eachdy] <- medyrmseratio
    #all_models_2more$medyrmseratio <- ifelse(all_models_2more$p_withheld == 10 & all_models_2more$BaseModel == eachbm & all_models_2more$DataYear == eachdy, medyrmseratio, all_models_2more$medyrmseratio)
  }
}

# save all_models_2more
write.csv(all_models_2more, file = "all_models_2more.csv")
all_models_2more <- read.table("all_models_2more.csv", header=T,sep=',',stringsAsFactors = FALSE)

# ------ Heatmap for 10 - 50 median Y RMSE ratio -------------
#df2yyrmse<-unique(all_models_2more$DataYear)
#df2x<-(unique(all_models_2more$BaseModel))
z_medyrmseratio <-round(reshape2::acast(all_models_2more,OrigModelYr~DataYear,fun=mean,value.var="medyrmseratio"),2)
means <- as.matrix(rowMeans(z_medyrmseratio))
out <- data.frame(z_medyrmseratio,means)
pdf("rmse_ratio_plot.pdf")
rmse_ratio <- ggplot(data = all_models_2more, mapping = aes(y = all_models_2more$OrigModelYr_ch, x = all_models_2more$DataYear_ch,fill = medyrmseratio)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(medyrmseratio,2))) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  #geom_raster(aes(all_models_2more$fitvarexp),hjust = 0.5,vjust = 0.5) +
  #aes(stringr::str_wrap("Model performance as a function of training data withheld", 15)) + #xlab(NULL) + ylab(yaxis_label) +
  labs(x="Data Year", y="Base Model Year",title="Model performance as a function of training data withheld", subtitle=stringr::str_wrap("Cell values = Ratio of median rmse for 100 iterations of 10% training data withheld / 50% withheld",width = 65)) +
  #theme(axis.text.x=element_text(angle=-45, hjust=0, vjust=1)) 
  #theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  #theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  #theme(legend.title = element_text(size=10, hjust=0.5, face=1, colour="black"),legend.justification = "center") + labs(fill = "Med Y RMSE Ratio") +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  #annotate(geom="text", x=6.25, y=6.35, label=round(rowMeans(z_medyrmseratio)[1],3), color="black",size=3) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(z_medyrmseratio)[6],2), color="black",size=4) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(z_medyrmseratio)[5],2), color="black",size=4) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(z_medyrmseratio)[4],2), color="black",size=4) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(z_medyrmseratio)[3],2), color="black",size=4) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(z_medyrmseratio)[2],2), color="black",size=4) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(z_medyrmseratio)[1],2), color="black",size=4) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(rmse_ratio)
dev.off()

  
# ------ Heatmap  for fitvarexp-----------
# bm.set <- all_models_2more[which(all_models_2more$BaseModel == eachbm),]
# #df2y<-unique(bm.set$DataYear)
# #df2x<-(unique(bm.set$OrigModelYr)) # this will be single value for single base model
# #df2x<-bm.set$OrigModelYr
# z_mean_fitvarexp <-reshape2::acast(bm.set,DataYear~OrigModelYr,fun=mean,value.var="fitvarexp")
# z_mean_fitvarexp
# # base package heatmap() reorders axis based on dendrogram
# heatmap(z_mean_fitvarexp,scale="column", col = heat.colors(256))

z_fitvarexp <-round(reshape2::acast(all_models_2more,OrigModelYr~DataYear,fun=mean,value.var="fitvarexp"),2)
fitmeans <- as.matrix(rowMeans(z_fitvarexp))
fitout <- data.frame(z_fitvarexp,fitmeans)
pdf("fitvarexp_plot.pdf")
fitvarexp_plot <- ggplot(data = all_models_2more, mapping = aes(y = all_models_2more$OrigModelYr_ch, x = all_models_2more$DataYear_ch,fill = fitvarexp)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(fitvarexp,0))) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  #geom_raster(aes(all_models_2more$fitvarexp),hjust = 0.5,vjust = 0.5) +
  #aes(stringr::str_wrap("Model performance as a function of training data withheld", 15)) + #xlab(NULL) + ylab(yaxis_label) +
  labs(x="Data Year", y="Base Model Year",title="Model % Variance Explained") +
  #theme(axis.text.x=element_text(angle=-45, hjust=0, vjust=1)) 
  #theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  #theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.title = element_text(size=10, hjust=0, face=1, colour="black"),legend.justification = "center") + labs(fill = "% Variance Explained") +
  theme(legend.position = "bottom") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  #annotate(geom="text", x=6.25, y=6.35, label=round(rowMeans(z_medyrmseratio)[1],3), color="black",size=3) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.86, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(z_fitvarexp)[6],0), color="black",size=4) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(z_fitvarexp)[5],0), color="black",size=4) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(z_fitvarexp)[4],0), color="black",size=4) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(z_fitvarexp)[3],0), color="black",size=4) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(z_fitvarexp)[2],0), color="black",size=4) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(z_fitvarexp)[1],00), color="black",size=4) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1))
  #theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(fitvarexp_plot)
dev.off()


# ------ Plotting combination table as a surface using cast & persp3d -------

surf3d <- all_models_2more[,c("p_withheld", "OrigModelYr", "fitvarexp")]
head(surf3d)
rgl::plot3d(surf3d) # generates points, not a surface

# List of variables that we are measuring against year/iterations/percen withheld, etc
#"fitvarexp","fitmse","medyrmse","varyrmse","medymbe","medymae","medmodelmse","DataYear","var1","var2","var3","iter")

# Establish x & y data
df2x<-unique(all_models_2more$DataYear)
df2y<-sort(unique(all_models_2more$OrigModelYr))

# Start the RGL 3D interactive popup
open3d()

# Establish the z data matrix:
z_mean_fitvarexp <-reshape2::acast(all_models_2more,DataYear~OrigModelYr,fun=mean,value.var="fitvarexp")
z_mean_fitvarexp

# Use persp3d to create a rotate-able graph in the RGL 3Dpopup:
persp3d(df2x,df2y,z_mean_fitvarexp, theta=40, phi=30, r=2, shade=0.75, axes=TRUE,scale=TRUE, box=TRUE, 
      nticks=40, ticktype="simple", zlim = c(55,80) , col=rainbow(50), ltheta = 120, lphi = 60,  xlab="DataYear", 
      ylab="OrigModelYr", zlab="Mean Fit Variance Exp", main="Mean Fit Variance Exp x OrigModelYr & DataYear")


# Use persp (no 3d) to create graph in the default Plots window
persp(df2x,df2y,z_mean_fitvarexp, theta=40, phi=30, r=2, shade=0.75, axes=TRUE,scale=TRUE, box=TRUE, 
        nticks=40, ticktype="simple", zlim = c(55,80) , col=rainbow(50), ltheta = 120, lphi = 60,  xlab="DataYear", 
        ylab="OrigModelYr", zlab="Mean Fit Variance Exp", main="Mean Fit Variance Exp x OrigModelYr & DataYear")

# ------ More graphing re percent withheld ---------

z_mean_medymae <-reshape2::acast(all_models_2more,year~p_withheld,fun=mean,value.var="medymae")
z_mean_medymae
persp(df2x,df2y,z_mean_medymae, theta=40, phi=30, r=2, shade=0.75, axes=TRUE,scale=TRUE, box=TRUE, 
      nticks=40, ticktype="simple", zlim = c(3.2,5.8) , col=rainbow(50), ltheta = 120, lphi = 60,  xlab="Year", 
      ylab="Percent withheld", zlab="Mean Med Y MAE", main="Mean Med Y MAE x % withheld & Year")

z_mean_medymbe <-reshape2::acast(all_models_2more,year~p_withheld,fun=mean,value.var="medymbe")
z_mean_medymbe
persp(df2x,df2y,z_mean_medymbe, theta=40, phi=30, r=2, shade=0.75, axes=TRUE,scale=TRUE, box=TRUE, 
      nticks=40, ticktype="simple", zlim = c(-.6,-.10) , col=rainbow(50), ltheta = 120, lphi = 60,  xlab="Year", 
      ylab="Percent withheld", zlab="Mean Med Y MBE", main="Mean Med Y MBE x % withheld & Year")
#"fitvarexp","fitmse","medyrmse","varyrmse","medymbe","medymae","medmodelmse","year","var1","var2","var3","var4","iter")


#  Iterations vs Percent Withheld
df2x<-unique(all_models_2more$iter)
df2y<-sort(unique(all_models_2more$p_withheld))

z_mean_medymae <-reshape2::acast(all_models_2more,iter~p_withheld,fun=mean,value.var="medymae")
z_mean_medymae
persp(df2x,df2y,z_mean_medymae, theta=40, phi=30, r=2, shade=0.75, axes=TRUE,scale=TRUE, box=TRUE, 
        nticks=40, ticktype="simple", zlim = c(3.9,5.4) , col=rainbow(50), ltheta = 120, lphi = 60,  xlab="Iteration", 
        ylab="Percent withheld", zlab="Mean Med Y MAE", main="Mean Med Y MAE x % withheld & Iterations")

z_mean_medymbe <-reshape2::acast(all_models_2more_dt,iter~p_withheld,fun=mean,value.var="medymbe")
z_mean_medymbe
persp3d(df2x,df2y,z_mean_medymbe, theta=40, phi=30, r=2, shade=0.75, axes=TRUE,scale=TRUE, box=TRUE, 
        nticks=40, ticktype="simple", zlim = c(-.5,-.1) , col=rainbow(50), ltheta = 120, lphi = 60,  xlab="iterations", 
        ylab="Percent withheld", zlab="Mean Med Y MBE", main="Mean Med Y MBE x % withheld & iterations")

z_mean_varyrmse <-reshape2::acast(all_models_2more_dt,iter~p_withheld,fun=mean,value.var="varyrmse")
z_mean_varyrmse
persp3d(df2x,df2y,z_mean_varyrmse, theta=40, phi=30, r=2, shade=0.75, axes=TRUE,scale=TRUE, box=TRUE, 
      nticks=40, ticktype="simple", zlim = c(.9,3.1) , col=rainbow(50), ltheta = 120, lphi = 60,  xlab="iterations", 
      ylab="Percent withheld", zlab="Mean Var Y RMSE", main="Mean Var Y RMSE x % withheld & iterations")

z_mean_medyrmse <-reshape2::acast(all_models_2more_dt,iter~p_withheld,fun=mean,value.var="medyrmse")
z_mean_medyrmse
persp(df2x,df2y,z_mean_medyrmse, theta=40, phi=30, r=2, shade=0.75, axes=TRUE,scale=TRUE, box=TRUE, 
      nticks=40, ticktype="simple", zlim = c(6,8.5) , col=rainbow(50), ltheta = 120, lphi = 60,  xlab="iterations", 
      ylab="Percent withheld", zlab="Mean Y RMSE", main="Mean Y RMSE x % withheld & iterations")

rgl::persp3d(df2x,df2y,z_mean_medyrmse, theta=40, phi=30, r=2, shade=0.75, axes=TRUE,scale=TRUE, box=TRUE, 
      nticks=40, ticktype="simple", zlim = c(6,8.5) , col=rainbow(50), ltheta = 120, lphi = 60,  xlab="iterations", 
      ylab="Percent withheld", zlab="Mean Y RMSE", main="Mean Y RMSE x % withheld & iterations")

z_mean_fitmse <-reshape2::acast(all_models_2more_dt,iter~p_withheld,fun=mean,value.var="fitmse")
z_mean_fitmse
persp3d(df2x,df2y,z_mean_fitmse, theta=40, phi=30, r=2, shade=0.75, axes=TRUE,scale=TRUE, box=TRUE, 
      nticks=40, ticktype="simple", zlim = c(441.5,443.5) , col=rainbow(50), ltheta = 120, lphi = 60,  xlab="iterations", 
      ylab="Percent withheld", zlab="Mean Fit MSE", main="Mean Fit MSE x % withheld & iterations")

z_meanfitvarexp <-reshape2::acast(all_models_2more_dt,iter~p_withheld,fun=mean,value.var="fitvarexp")
z_meanfitvarexp
persp3d(df2x,df2y,z_meanfitvarexp, theta=40, phi=30, r=2, shade=0.75, axes=TRUE,scale=TRUE, box=TRUE, 
      nticks=40, ticktype="simple", zlim = c(71.1,71.3) , col=rainbow(50), ltheta = 120, lphi = 60,  xlab="iterations", 
      ylab="Percent withheld", zlab="Mean Fit Var Exp", main="Mean Fit Var Exp x % withheld & iterations")

z_meanmedmodelmse <-reshape2::acast(all_models_2more_dt,iter~p_withheld,fun=mean,value.var="medmodelmse")
z_meanmedmodelmse
persp(df2x,df2y,z_meanmedmodelmse, theta=70, phi=30, r=2, shade=0.75, axes=TRUE,scale=TRUE, box=TRUE, 
      nticks=40, ticktype="simple", zlim = c(40,80) , col=rainbow(50), ltheta = 120, lphi = 60,  xlab="iterations", 
      ylab="Percent withheld", zlab="Mean Med model MSE", main="Mean med.model.mse x % withheld & iterations")

z_summedmodelmse <-reshape2::acast(all_models_2more_dt,iter~p_withheld,fun=sum,value.var="medmodelmse")
z_summedmodelmse
persp(df2x,df2y,z_summedmodelmse, theta=80, phi=30, r=2, shade=0.75, axes=TRUE,scale=TRUE, box=TRUE, 
      nticks=40, ticktype="simple", zlim = c(1500,2600) , col=rainbow(50), ltheta = 120, lphi = 60,  xlab="iterations", 
      ylab="Percent withheld", zlab="SUM Med model MSE", main="Sum med.model.mse x % withheld & iterations")

# open3d()
# surface3d(df2x,df2y,df2z)
# persp(df2x,df2y,z_meanmedmodelmse, theta=30, phi=30, r=2, shade=0.75, axes=TRUE,scale=TRUE, box=TRUE, 
#       nticks=40, ticktype="simple", zlim = c(40,80) , col=rainbow(50), ltheta = 120, lphi = 60,  xlab="iterations", 
#       ylab="Percent withheld", zlab="Med model MSE", main="Mean med.model.mse x % withheld & iterations")



# ------ Plotting - scatterplot with labels  ----
library(ggplot2)
library(plotly)
pl2 <- ggplot(data=all_models_2more, aes(x=iter, y=medyrmse)) + 
  # geom_point() +
  xlab("Iterations") + ylab("Median Y RMSE") +
  # geom_text(aes(label=MMI_4$ROIID), size = 3) +
  geom_jitter(width = 0.01) +
  scale_colour_manual(values = c("Dark Green", "Red", "Green", "Orange"))
print(pl2)
#NOTE: do not do geom_point and geom_jitter in same plot, or points are duplicated

# need to use R & not RStudio for plot_ly plots
#plot_ly(data = all_models_2more,x = all_models_2more$iter,y = all_models_2more$p_withheld, z = all_models_2more$fit.var.exp,type = "scatter3d",showlegend = FALSE)

# ------ Manipulate package - interesting stuff for later ---------
manipulate(plot(1:x), x = slider(1, 100))
manipulate(
  plot(cars, xlim=c(0,x.max)),  
  x.max=slider(15,25))

# ------ Plotting - boxplots -------------
ggplot(data = all_models_2more,aes(x=iter,y=med.y.rmse) ) + geom_boxplot(aes(group=iter)) +
  labs(title="Grouped by Percent Withheld",x="p_withheld", y = "Median Y RMSE")
#xlab("p_withheld") + ylab("Median Y RMSE") + title("Iterations by RMSE")
ggplot(data = all_models_2more,aes(x=iter,y=med.y.rmse) ) + geom_boxplot(aes(group=OrigModelYr))
#ggplot(data = all_models_2more, aes(x=iter,y=med.y.rmse) ) + geom_boxplot(aes(fill=p_withheld)) + facet_wrap( ~ variable, scales="free")

# ------ Set of 7 boxplots -----------
summary(all_models_2more)
pdfname <- "boxplots.pdf"
pdf(pdfname) 
#x11(width=7, height=8)
plot.new()
par(mfrow=c(4,2),oma=c(7,3,2,2.5),mai=c(0,.1,0,0),tck=-0.025)

ii = 1
for(ii in 1:7) {
  # Set relevant parameters for each box
  if(ii == 1){
    ylab. <- "fit.var.exp"
    yfact <- "fit.var.exp"
    ylim. <- c(55,100)
    at. <- c(50,75,100)
    axis. <- 2  
    metric. <- all_models_2more$fit.var.exp
  }
  if(ii == 2){
    ylab. <- "fit.mse"
    yfact <- "fit.mse"
    ylim. <- c(300,700)      
    at. <- c(250,500,700)
    axis. <- 4          
    metric. <- all_models_2more$fit.mse
  }
  if(ii == 3){
    ylab. <- "med.y.rmse"
    yfact <- "med.y.rmse"
    ylim. <- c(3,12)      
    at. <- c(0,7.5,15)
    axis. <- 2
    metric. <- all_models_2more$med.y.rmse
  }
  
  if(ii == 4){
    ylab. <- "var.y.rmse"
    yfact <- "var.y.rmse"
    ylim. <- c(0,9)      
    at. <- c(0,2.5,5,9)
    axis. <- 4
    metric. <- all_models_2more$var.y.rmse  
  }
  
  if(ii == 5){
    ylab. <- "med.y.mbe"
    yfact <- "med.y.mbe"
    ylim. <- c(-2.5,2.5)      
    at. <- c(-2.5,0,2.5)
    axis. <- 2
    metric. <- all_models_2more$med.y.mbe  
  }
  if(ii == 6){
    ylab. <- "med.y.mae"
    yfact <- "med.y.mae"
    ylim. <- c(2,7)      
    at. <- c(2,5,7)
    axis. <- 4
    metric. <- all_models_2more$med.y.mae  
  }
  if(ii == 7){
    ylab. <- "med.model.mse"
    yfact <- "med.model.mse"
    ylim. <- c(30,100)      
    at. <- c(25,50,75,100)
    axis. <- 2
    metric. <- all_models_2more$med.model.mse  
  }
  # Boxplots
  bp <- boxplot(get(yfact) ~ p_withheld, data=all_models_2more, 
                ylim=ylim., xaxt='n', yaxt='n')
  
  legend(y=max(ylim.)+(0.075*max(ylim.)), x=0, ylab., box.lty=0, cex=1.4)
  axis(axis., at=at., labels=T, tick=T, outer=T, las=2)          
}
mtext(side=1, paste("Boxplots by Percent Withheld",sep=""), outer=T, font=1, cex=1.0, line=3)
dev.off()

# ------ Predict/Create rasters with Random Forest from RF, EBP, & BR models -------------------------
setwd("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\predictions20200203")
# U11 Predict herb_2009 with 2009 rasters & 2018 EBP model vars - ARGmea0906 + NDVmax0903 + ARGstd0903
rfmodel <- randomForest(heb_2009 ~ NDVstd0903 + ARGmea0906 + B3std0906, 
                   data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmea0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea0906_1_ARG_mean_u11.tif")
NDVmax0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmax0903_0_NDVI_max_u11.tif")
ARGstd0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGstd0903_0_ARG_stdDev_u11.tif")
outstack <- raster::stack(ARGmea0906,NDVmax0903,ARGstd0903); names(outstack) <- c("ARGmea0906","NDVmax0903","ARGstd0903")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2009_2018bm_u11.tif", format='GTiff',overwrite=TRUE)

# U10
ARGmea0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea0906_1_ARG_mean_u10.tif")
NDVmax0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmax0903_0_NDVI_max_u10.tif")
ARGstd0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGstd0903_0_ARG_stdDev_u10.tif")
outstack <- raster::stack(ARGmea0906,NDVmax0903,ARGstd0903); names(outstack) <- c("ARGmea0906","NDVmax0903","ARGstd0903")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2009_2018bm_u10.tif", format='GTiff',overwrite=TRUE)

# U11 Predict herb_2009 with 2009 rasters & 2009 RF model vars - ARGmax0903 + NBRstd0903 + ARGmed0906
rfmodel <- randomForest(heb_2009 ~ ARGmax0903 + NBRstd0903 + ARGmed0906, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax0903_0_ARG_max_u11.tif")
NBRstd0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd0903_0_NBR_stdDev_u11.tif")
ARGmed0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed0906_1_ARG_median_u11.tif")
outstack <- raster::stack(ARGmax0903,NBRstd0903,ARGmed0906); names(outstack) <- c("ARGmax0903","NBRstd0903","ARGmed0906")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2009_top3_u11.tif", format='GTiff',overwrite=TRUE)

# U10
ARGmax0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax0903_0_ARG_max_u10.tif")
NBRstd0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd0903_0_NBR_stdDev_u10.tif")
ARGmed0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed0906_1_ARG_median_u10.tif")
outstack <- raster::stack(ARGmax0903,NBRstd0903,ARGmed0906); names(outstack) <- c("ARGmax0903","NBRstd0903","ARGmed0906")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2009_top3_u10.tif", format='GTiff',overwrite=TRUE)

# U11 Predict herb_2009 with 2009 rasters & 2009 BR model vars -  B3std0909 + NBRstd0903 + ARGmed0906
rfmodel <- randomForest(heb_2009 ~  B3std0909 + NBRstd0903 + ARGmed0906, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
B3std0909 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std0909_2_B3_stdDev_u11.tif")
NBRstd0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd0903_0_NBR_stdDev_u11.tif")
ARGmed0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed0906_1_ARG_median_u11.tif")
outstack <- raster::stack(B3std0909,NBRstd0903,ARGmed0906); names(outstack) <- c("B3std0909","NBRstd0903","ARGmed0906")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "br2009_top3_u11.tif", format='GTiff',overwrite=TRUE)

# U10
B3std0909 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std0909_2_B3_stdDev_u10.tif")
NBRstd0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd0903_0_NBR_stdDev_u10.tif")
ARGmed0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed0906_1_ARG_median_u10.tif")
outstack <- raster::stack(B3std0909,NBRstd0903,ARGmed0906); names(outstack) <- c("B3std0909","NBRstd0903","ARGmed0906")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "br2009_top3_u10.tif", format='GTiff',overwrite=TRUE)


# Predict herb_2009 with 2009 rasters & 2009 RF top 4 model vars: ARGmax0903 + NBRstd0903 + ARGmed0906 + B3std0909
# U11 
rfmodel <- randomForest(heb_2009 ~  ARGmax0903 + NBRstd0903 + ARGmed0906 + B3std0909, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax0903_0_ARG_max_u11.tif")
B3std0909 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std0909_2_B3_stdDev_u11.tif")
NBRstd0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd0903_0_NBR_stdDev_u11.tif")
ARGmed0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed0906_1_ARG_median_u11.tif")
outstack <- raster::stack(ARGmax0903,B3std0909,NBRstd0903,ARGmed0906); names(outstack) <- c("ARGmax0903","B3std0909","NBRstd0903","ARGmed0906")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2009_top4_u11.tif", format='GTiff',overwrite=TRUE)

# U10
ARGmax0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax0903_0_ARG_max_u10.tif")
B3std0909 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std0909_2_B3_stdDev_u10.tif")
NBRstd0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd0903_0_NBR_stdDev_u10.tif")
ARGmed0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed0906_1_ARG_median_u10.tif")
outstack <- raster::stack(ARGmax0903,B3std0909,NBRstd0903,ARGmed0906); names(outstack) <- c("ARGmax0903","B3std0909","NBRstd0903","ARGmed0906")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2009_top4_u10.tif", format='GTiff',overwrite=TRUE)




# U11 Predict herb_2011 with 2011 rasters & 2018 EBP model vars - ARGmea1106 + NDVmax1103 + ARGstd1103
rfmodel <- randomForest(heb_2011 ~ ARGmea1106 + NDVmax1103 + ARGstd1103, 
                        data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmea1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1106_1_ARG_mean_u11.tif")
NDVmax1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmax1103_0_NDVI_max_u11.tif")
ARGstd1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGstd1103_0_ARG_stdDev_u11.tif")
outstack <- raster::stack(ARGmea1106,NDVmax1103,ARGstd1103); names(outstack) <- c("ARGmea1106","NDVmax1103","ARGstd1103")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2011_2018bm_u11.tif", format='GTiff',overwrite=TRUE)

# U10
ARGmea1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1106_1_ARG_mean_u10.tif")
NDVmax1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmax1103_0_NDVI_max_u10.tif")
ARGstd1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGstd1103_0_ARG_stdDev_u10.tif")
outstack <- raster::stack(ARGmea1106,NDVmax1103,ARGstd1103); names(outstack) <- c("ARGmea1106","NDVmax1103","ARGstd1103")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2011_2018bm_u10.tif", format='GTiff',overwrite=TRUE)


# U11 Predict herb_2011 with 2011 rasters & top 3 RF & BR model vars - ARGmea1106 + B3std1106 + NDVstd1103
rfmodel <- randomForest(heb_2011 ~ ARGmea1106 + B3std1106 + NDVstd1103, 
                        data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmea1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1106_1_ARG_mean_u11.tif")
B3std1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1106_1_B3_stdDev_u11.tif")
NDVstd1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1103_0_NDVI_stdDev_u11.tif")
outstack <- raster::stack(ARGmea1106,B3std1106,NDVstd1103); names(outstack) <- c("ARGmea1106","B3std1106","NDVstd1103")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2011_rfbrtop3_u11.tif", format='GTiff',overwrite=TRUE)

# U10
ARGmea1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1106_1_ARG_mean_u10.tif")
B3std1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1106_1_B3_stdDev_u10.tif")
NDVstd1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1103_0_NDVI_stdDev_u10.tif")
outstack <- raster::stack(ARGmea1106,B3std1106,NDVstd1103); names(outstack) <- c("ARGmea1106","B3std1106","NDVstd1103")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2011_rfbrtop3_u10.tif", format='GTiff',overwrite=TRUE)


# U11 Predict herb_2013 with 2013 rasters & 2018 EBP model vars - ARGmea1306 + NDVmax1303 + ARGstd1303
rfmodel <- randomForest(heb_2013 ~ ARGmea1306 + NDVmax1303 + ARGstd1303, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmea1306 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1306_1_ARG_mean_u11.tif")
NDVmax1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmax1303_0_NDVI_max_u11.tif")
ARGstd1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGstd1303_0_ARG_stdDev_u11.tif")
outstack <- raster::stack(ARGmea1306,NDVmax1303,ARGstd1303); names(outstack) <- c("ARGmea1306","NDVmax1303","ARGstd1303")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2013_2018bm_u11.tif", format='GTiff',overwrite=TRUE)

# U10
ARGmea1306 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1306_1_ARG_mean_u10.tif")
NDVmax1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmax1303_0_NDVI_max_u10.tif")
ARGstd1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGstd1303_0_ARG_stdDev_u10.tif")
outstack <- raster::stack(ARGmea1306,NDVmax1303,ARGstd1303); names(outstack) <- c("ARGmea1306","NDVmax1303","ARGstd1303")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2013_2018bm_u10.tif", format='GTiff',overwrite=TRUE)


# U11 Predict herb_2013 with 2013 rasters & top 3 RF & BR model vars - ARGmea1306 + NDVstd1303 + B3max1309
rfmodel <- randomForest(heb_2013 ~ ARGmea1306 + NDVstd1303 + B3max1309, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmea1306 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1306_1_ARG_mean_u11.tif")
B3max1309 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3max1309_2_B3_max_u11.tif")
NDVstd1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1303_0_NDVI_stdDev_u11.tif")
outstack <- raster::stack(ARGmea1306,B3max1309,NDVstd1303); names(outstack) <- c("ARGmea1306","B3max1309","NDVstd1303")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2013_rfbrtop3_u11.tif", format='GTiff',overwrite=TRUE)

# U10
ARGmea1306 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1306_1_ARG_mean_u10.tif")
B3max1309 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3max1309_2_B3_max_u10.tif")
NDVstd1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1303_0_NDVI_stdDev_u10.tif")
outstack <- raster::stack(ARGmea1306,B3max1309,NDVstd1303); names(outstack) <- c("ARGmea1306","B3max1309","NDVstd1303")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2013_rfbrtop3_u10.tif", format='GTiff',overwrite=TRUE)


# U11 Predict herb_2015 with 2015 rasters & 2018 EBP model vars - ARGmea1506 + NDVmax1503 + ARGstd1503
rfmodel <- randomForest(heb_2015 ~ ARGmea1506 + NDVmax1503 + ARGstd1503, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmea1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1506_1_ARG_mean_u11.tif")
NDVmax1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmax1503_0_NDVI_max_u11.tif")
ARGstd1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGstd1503_0_ARG_stdDev_u11.tif")
outstack <- raster::stack(ARGmea1506,NDVmax1503,ARGstd1503); names(outstack) <- c("ARGmea1506","NDVmax1503","ARGstd1503")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2015_2018bm_u11.tif", format='GTiff',overwrite=TRUE)

# U10
ARGmea1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1506_1_ARG_mean_u10.tif")
NDVmax1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmax1503_0_NDVI_max_u10.tif")
ARGstd1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGstd1503_0_ARG_stdDev_u10.tif")
outstack <- raster::stack(ARGmea1506,NDVmax1503,ARGstd1503); names(outstack) <- c("ARGmea1506","NDVmax1503","ARGstd1503")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2015_2018bm_u10.tif", format='GTiff',overwrite=TRUE)


# U11 Predict herb_2015 with 2015 rasters & top 3 RF model vars - ARGmax1503 + ARGmax1506 + NDVmin1509
rfmodel <- randomForest(heb_2015 ~ ARGmax1503 + ARGmax1506 + NDVmin1509, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1503_0_ARG_max_u11.tif")
ARGmax1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1506_1_ARG_max_u11.tif")
NDVmin1509 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmin1509_2_NDVI_min_u11.tif")
outstack <- raster::stack(ARGmax1503,ARGmax1506,NDVmin1509); names(outstack) <- c("ARGmax1503","ARGmax1506","NDVmin1509")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
head(outpredict)
values(outpredict)[values(outpredict) > 31.594 & values(outpredict) < 31.595] <- 0
writeRaster(outpredict, "rf2015_rftop3_u11.tif", format='GTiff',overwrite=TRUE)
# U10
ARGmax1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1503_0_ARG_max_u10.tif")
ARGmax1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1506_1_ARG_max_u10.tif")
NDVmin1509 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmin1509_2_NDVI_min_u10.tif")
outstack <- raster::stack(ARGmax1503,ARGmax1506,NDVmin1509); names(outstack) <- c("ARGmax1503","ARGmax1506","NDVmin1509")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
head(outpredict)
values(outpredict)[values(outpredict) > 31.594 & values(outpredict) < 31.595] <- 0
head(outpredict)
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2015_rftop3_u10.tif",overwrite=TRUE)



# Predict herb_2015 with 2015 rasters & top 4 RF model vars - ARGmax1503 + ARGmed1506 + ARGmax1506 + NBRstd1503
# U11 
rfmodel <- randomForest(heb_2015 ~ ARGmax1503 + ARGmed1506 + ARGmax1506 + NBRstd1503, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1503_0_ARG_max_u11.tif")
ARGmax1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1506_1_ARG_max_u11.tif")
ARGmed1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1506_1_ARG_med_u11.tif")
NBRstd1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1503_0_NBR_stdDev_u11.tif")
outstack <- raster::stack(ARGmax1503,ARGmax1506,ARGmed1506,NBRstd1503); names(outstack) <- c("ARGmax1503","ARGmax1506","ARGmed1506","NBRstd1503")

start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
head(outpredict)
#values(outpredict)[values(outpredict) > 31.594 & values(outpredict) < 31.595] <- 0
writeRaster(outpredict, "rf2015_rftop4_u11.tif", format='GTiff',overwrite=TRUE)

# U10
ARGmax1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1503_0_ARG_max_u10.tif")
ARGmax1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1506_1_ARG_max_u10.tif")
ARGmed1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1506_1_ARG_med_u10.tif")
NBRstd1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1503_0_NBR_stdDev_u10.tif")
outstack <- raster::stack(ARGmax1503,ARGmax1506,ARGmed1506,NBRstd1503); names(outstack) <- c("ARGmax1503","ARGmax1506","ARGmed1506","NBRstd1503")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
head(outpredict)
#values(outpredict)[values(outpredict) > 31.594 & values(outpredict) < 31.595] <- 0
#head(outpredict)
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2015_rftop4_u10.tif",overwrite=TRUE)

# 
# # U11 Predict herb_2015 with 2015 rasters & top 3 BR model vars - NDImin1503 + NDVmax1503 + B3mea1509
# rfmodel <- randomForest(heb_2015 ~ NDImin1503 + NDVmax1503 + B3mea1509, 
#                         data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
# NDImin1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDImin1503_0_NDII_min_u11.tif")
# NDVmax1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmax1503_1_NDVI_max_u11.tif")
# B3mea1509 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3mea1509_2_B3_mean_u11.tif")
# outstack <- raster::stack(NDImin1503,NDVmax1503,B3mea1509); names(outstack) <- c("NDImin1503","NDVmax1503","B3mea1509")
# start = Sys.time();beginCluster(36)
# outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
# endCluster();    finish = Sys.time(); finish-start
# writeRaster(outpredict, "rf2015_brtop3_u11.tif", format='GTiff',overwrite=TRUE)
# 
# # U10
# NDImin1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\u10\\NDImin1503_0_NDII_min_u10.tif")
# NDVmax1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\u10\\NDVmax1503_1_NDVI_max_u10.tif")
# B3mea1509 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\u10\\B3mea1509_2_B3_mean_u10.tif")
# outstack <- raster::stack(NDImin1503,NDVmax1503,B3mea1509); names(outstack) <- c("NDImin1503","NDVmax1503","B3mea1509")
# start = Sys.time();beginCluster(36)
# outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
# endCluster();    finish = Sys.time(); finish-start
# writeRaster(outpredict, "rf2015_brtop3_u10.tif", format='GTiff',overwrite=TRUE)


# Predict herb_2017 with 2017 rasters & 2018 EBP model vars - ARGmea1706 + NDVmax1703 + ARGstd1703
# U11 
rfmodel <- randomForest(heb_2017 ~ ARGmea1706 + NDVmax1703 + ARGstd1703, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmea1706 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1706_1_ARG_mean_u11.tif")
NDVmax1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmax1703_0_NDVI_max_u11.tif")
ARGstd1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGstd1703_0_ARG_stdDev_u11.tif")
outstack <- raster::stack(ARGmea1706,NDVmax1703,ARGstd1703); names(outstack) <- c("ARGmea1706","NDVmax1703","ARGstd1703")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2017_2018bm_u11.tif", format='GTiff',overwrite=TRUE)

# U10
ARGmea1706 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1706_1_ARG_mean_u10.tif")
NDVmax1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmax1703_0_NDVI_max_u10.tif")
ARGstd1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGstd1703_0_ARG_stdDev_u10.tif")
outstack <- raster::stack(ARGmea1706,NDVmax1703,ARGstd1703); names(outstack) <- c("ARGmea1706","NDVmax1703","ARGstd1703")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2017_2018bm_u10.tif", format='GTiff',overwrite=TRUE)


# U11 Predict herb_2017 with 2017 rasters & top 3 RF model vars - ARGmax1706 + NBRstd1703 + ARGmed1709
rfmodel <- randomForest(heb_2017 ~ ARGmax1706 + NBRstd1703 + ARGmed1709, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1706 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1706_1_ARG_max_u11.tif")
NBRstd1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1703_0_NBR_stdDev_u11.tif")
ARGmed1709 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1709_2_ARG_median_u11.tif")
outstack <- raster::stack(ARGmax1706,NBRstd1703,ARGmed1709); names(outstack) <- c("ARGmax1706","NBRstd1703","ARGmed1709")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2017_rfbrtop3_u11.tif", format='GTiff',overwrite=TRUE)

# U10
ARGmax1706 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1706_1_ARG_max_u10.tif")
NBRstd1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1703_0_NBR_stdDev_u10.tif")
ARGmed1709 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1709_2_ARG_median_u10.tif")
outstack <- raster::stack(ARGmax1706,NBRstd1703,ARGmed1709); names(outstack) <- c("ARGmax1706","NBRstd1703","ARGmed1709")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2017_rfbrtop3_u10.tif", format='GTiff',overwrite=TRUE)


# U11 Predict Herb_2018 with 2018 rasters & 2018 EBP model vars - ARGmea1806 + NDVmax1803 + ARGstd1803
rfmodel <- randomForest(Herb_2018 ~ ARGmea1806 + NDVmax1803 + ARGstd1803, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmea1806 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1806_1_ARG_mean_u11.tif")
NDVmax1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmax1803_0_NDVI_max_u11.tif")
ARGstd1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGstd1803_0_ARG_stdDev_u11.tif")
outstack <- raster::stack(ARGmea1806,NDVmax1803,ARGstd1803); names(outstack) <- c("ARGmea1806","NDVmax1803","ARGstd1803")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2018_2018bm_u11.tif", format='GTiff',overwrite=TRUE)

# U10
ARGmea1806 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1806_1_ARG_mean_u10.tif")
NDVmax1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmax1803_0_NDVI_max_u10.tif")
ARGstd1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGstd1803_0_ARG_stdDev_u10.tif")
outstack <- raster::stack(ARGmea1806,NDVmax1803,ARGstd1803); names(outstack) <- c("ARGmea1806","NDVmax1803","ARGstd1803")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2018_2018bm_u10.tif", format='GTiff',overwrite=TRUE)


# U11 Predict Herb_2018 with 2018 rasters & top 3 RF model vars - ARGmax1806 + NBRstd1803 +  B3std1803
rfmodel <- randomForest(Herb_2018 ~ ARGmax1806 + NBRstd1803 +  B3std1803, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1806 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1806_1_ARG_max_u11.tif")
NBRstd1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1803_0_NBR_stdDev_u11.tif")
B3std1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1803_0_B3_stdDev_u11.tif")
outstack <- raster::stack(ARGmax1806,NBRstd1803,B3std1803); names(outstack) <- c("ARGmax1806","NBRstd1803","B3std1803")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2018_rftop3_u11.tif", format='GTiff',overwrite=TRUE)

#U10
ARGmax1806 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1806_1_ARG_max_u10.tif")
NBRstd1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1803_0_NBR_stdDev_u10.tif")
B3std1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1803_0_B3_stdDev_u10.tif")
outstack <- raster::stack(ARGmax1806,NBRstd1803,B3std1803); names(outstack) <- c("ARGmax1806","NBRstd1803","B3std1803")
start = Sys.time();beginCluster(45)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2018_rftop3_u10.tif", format='GTiff',overwrite=TRUE)








# ------ Betareg   -----------------------------------------------------------------
# Can transform data into 0 & 1 with bernouli, but might not be as helpful
# Need to scale dependent to bw 0 & 1

str2useinmodel
detach()

# Create scaled field hpsc & populate based on herb pct
y09_nofires$hpsc <- y09_nofires$heb_2009 * .01
y09_nofires$hpsc[y09_nofires$hpsc == 1] <- .99
y09_nofires$hpsc[y09_nofires$hpsc == 0] <- .001
y11_nofires$hpsc <- y11_nofires$heb_2011 * .01
y11_nofires$hpsc[y11_nofires$hpsc == 1] <- .99
y11_nofires$hpsc[y11_nofires$hpsc == 0] <- .001
y13_nofires$hpsc <- y13_nofires$heb_2013 * .01
y13_nofires$hpsc[y13_nofires$hpsc == 1] <- .99
y13_nofires$hpsc[y13_nofires$hpsc == 0] <- .001
y15_nofires$hpsc <- y15_nofires$heb_2015 * .01
y15_nofires$hpsc[y15_nofires$hpsc == 1] <- .99
y15_nofires$hpsc[y15_nofires$hpsc == 0] <- .001
y17_nofires$hpsc <- y17_nofires$heb_2017 * .01
y17_nofires$hpsc[y17_nofires$hpsc == 1] <- .99
y17_nofires$hpsc[y17_nofires$hpsc == 0] <- .001
y18_nofires$hpsc <- y18_nofires$Herb_2018 * .01
y18_nofires$hpsc[y18_nofires$hpsc == 1] <- .99
y18_nofires$hpsc[y18_nofires$hpsc == 0] <- .001

beta_model <- betareg(hpsc ~ NDVstd0903 + ARGmea0906 + B3std0906, na.action=na.exclude, link = "probit",data=y09_nofires)
summary(beta_model)
bptest(beta_model, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test 
betareg_predict <- cbind(
  predict(beta_model, type = "response"),
  predict(beta_model, type = "link"),
  predict(beta_model, type = "precision"),
  predict(beta_model, type = "variance"),
  predict(beta_model, type = "quantile", at = c(0.25, 0.5, 0.75))
)
colnames(betareg_predict) <- c("brresponse","link","precision","variance","quantilep25","quantilep5","quantilep75")
aa <- cbind(betareg_predict,y09_nofires)
paste("rmse = ", round(rmse(aa$hpsc,aa$brresponse),4))
#  plot betareg output w ggplot
ggplot(aa, aes(x=hpsc, y=brresponse)) +
  geom_point() + geom_text(aes(label=pltnum)) +
  scale_fill_grey() +
  scale_colour_manual("", values = c("red", "blue")) +
  scale_linetype_manual("", values = c("solid", "dashed")) +
  geom_smooth(method = "lm", se = TRUE,colour = "red") +
  abline(0,1,col = "green") #abline not showing

detach()
attach(y18_nofires)
plot(y13_nofires$hpsc, fitted(beta_model), ylim = c(0,1), xlim = c(0,1),main = "betareg model")
c <- lm(hpsc ~ ARGmed1706 + NBRstd1703 + ARGmed1709)
summary(c); vif(c)



# ------ Compare & plot rf repsonse vs betareg response ---------------------------
# use aarf from rf output
bresp <- as.data.frame(cbind(aarf$rfresponse,aa))
#colnames(bresp[,1]) <- c("rfresponse")
ggplot(bresp, aes(x=aarf$rfresponse, y=brresponse)) +
  geom_point() + geom_text(aes(label=pltnum)) +
  scale_fill_grey() +
  scale_colour_manual("", values = c("red", "blue")) +
  scale_linetype_manual("", values = c("solid", "dashed")) +
  geom_smooth(method = "lm", se = TRUE,colour = "red")


paste("rmse betareg vs random forest= ", round(rmse(bresp$`aarf$rfresponse`,bresp$brresponse),4))

# Generalized Leverages plot
glev <- gleverage(beta_model);glev
aaglev <- cbind(aa,glev)
ggplot(aaglev, aes(x=glev, y=hpsc)) +
  geom_point() + geom_text(aes(label=pltnum)) +
  scale_fill_grey() +
  scale_colour_manual("", values = c("red", "blue")) +
  scale_linetype_manual("", values = c("solid", "dashed")) +
  geom_smooth(method = "lm", se = TRUE,colour = "red") +
  abline(0,1,col = "green") #abline not showing
theme_bw()




# ------ Predict to raster best models from year applied to other years ------------------------------------------
models_cast <-round(reshape2::acast(all_models_2more,IndVars~OrigModelYr,fun=mean,value.var="OrigModelYr",drop=TRUE,fill=NULL),2)
models_cast
setwd("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\predictions20200616")

# U11 Predict 2009 model on all years ------------
# U11 Predict 2009 model with 2009 data: heb_2009 ~ ARGmax0903 + NBRstd0903 + ARGmed0906
rfmodel <- randomForest(heb_2009 ~ ARGmax0903 + NBRstd0903 + ARGmed0906, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax0903_0_ARG_max_u11.tif")
NBRstd0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd0903_0_NBR_std_u11.tif")
#values(NBRmax0903)[values(NBRmax0903) == 0] = NA
ARGmed0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed0906_1_ARG_med_u11.tif")
outstack <- raster::stack(ARGmax0903,NBRstd0903,ARGmed0906); names(outstack) <- c("ARGmax0903","NBRstd0903","ARGmed0906")
start = Sys.time();beginCluster(8)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
start = Sys.time() # no clustering
outpredict <- raster::predict(outstack,rfmodel,fun=predict,na.rm=TRUE)
finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2009_2009bm_u11.tif", format='GTiff',overwrite=TRUE)


# U11 Predict 2009 model with 2011 data: heb_2011 ~ ARGmax1103 + NBRstd1103 + ARGmed1106
rfmodel <- randomForest(heb_2011 ~ ARGmax1103 + NBRstd1103 + ARGmed1106, 
                        data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1103_0_ARG_max_u11.tif")
NBRstd1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1103_0_NBR_std_u11.tif")
ARGmed1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1106_1_ARG_med_u11.tif")
outstack <- raster::stack(ARGmax1103,NBRstd1103,ARGmed1106); names(outstack) <- c("ARGmax1103","NBRstd1103","ARGmed1106")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
start = Sys.time() # no clustering
outpredict <- raster::predict(outstack,rfmodel,fun=predict,na.rm=TRUE)
finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2011_2009bm_u11.tif", format='GTiff',overwrite=TRUE)

# U11 Predict 2009 model with 2013 data: heb_2013 ~ ARGmax1303 + NBRstd1303 + ARGmed1306
rfmodel <- randomForest(heb_2013 ~ ARGmax1303 + NBRstd1303 + ARGmed1306, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1303_0_ARG_max_u11.tif")
NBRstd1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1303_0_NBR_std_u11.tif")
ARGmed1306 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1306_1_ARG_med_u11.tif")
outstack <- raster::stack(ARGmax1303,NBRstd1303,ARGmed1306); names(outstack) <- c("ARGmax1303","NBRstd1303","ARGmed1306")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2013_2009bm_u11.tif", format='GTiff',overwrite=TRUE)

# U11 Predict 2009 model with 2015 data: heb_2015 ~ ARGmax1503 + NBRstd1503 + ARGmed1506
rfmodel <- randomForest(heb_2015 ~ ARGmax1503 + NBRstd1503 + ARGmed1506, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1503_0_ARG_max_u11.tif")
NBRstd1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1503_0_NBR_std_u11.tif")
ARGmed1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1506_1_ARG_med_u11.tif")
outstack <- raster::stack(ARGmax1503,NBRstd1503,ARGmed1506); names(outstack) <- c("ARGmax1503","NBRstd1503","ARGmed1506")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2015_2009bm_u11.tif", format='GTiff',overwrite=TRUE)

# U11 Predict 2009 model with 2017 data: heb_2017 ~ ARGmax1703 + NBRstd1703 + ARGmed1706
rfmodel <- randomForest(heb_2017 ~ ARGmax1703 + NBRstd1703 + ARGmed1706, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1703_0_ARG_max_u11.tif")
NBRstd1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1703_0_NBR_std_u11.tif")
ARGmed1706 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1706_1_ARG_med_u11.tif")
outstack <- raster::stack(ARGmax1703,NBRstd1703,ARGmed1706); names(outstack) <- c("ARGmax1703","NBRstd1703","ARGmed1706")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2017_2009bm_u11.tif", format='GTiff',overwrite=TRUE)

# U11 Predict 2009 model with 2018 data: Herb_2018 ~ ARGmax1803 + NBRstd1803 + ARGmed1806
rfmodel <- randomForest(Herb_2018 ~ ARGmax1803 + NBRstd1803 + ARGmed1806, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1803_0_ARG_max_u11.tif")
NBRstd1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1803_0_NBR_std_u11.tif")
ARGmed1806 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1806_1_ARG_med_u11.tif")
outstack <- raster::stack(ARGmax1803,NBRstd1803,ARGmed1806); names(outstack) <- c("ARGmax1803","NBRstd1803","ARGmed1806")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2018_2009bm_u11.tif", format='GTiff',overwrite=TRUE)

# U10 Predict 2009 model on all years ------------
# U10 Predict 2009 model with 2009 data: heb_2009 ~ ARGmax0903 + NBRstd0903 + ARGmed0906
rfmodel <- randomForest(heb_2009 ~ ARGmax0903 + NBRstd0903 + ARGmed0906, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax0903_0_ARG_max_u10.tif")
NBRstd0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd0903_0_NBR_std_u10.tif")
ARGmed0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed0906_1_ARG_med_u10.tif")
outstack <- raster::stack(ARGmax0903,NBRstd0903,ARGmed0906); names(outstack) <- c("ARGmax0903","NBRstd0903","ARGmed0906")
start = Sys.time();beginCluster(8)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()  
writeRaster(outpredict, "rf2009_2009bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2009 model with 2011 data: heb_2011 ~ ARGmax1103 + NBRstd1103 + ARGmed1106
rfmodel <- randomForest(heb_2011 ~ ARGmax1103 + NBRstd1103 + ARGmed1106, 
                        data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1103_0_ARG_max_u10.tif")
NBRstd1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1103_0_NBR_std_u10.tif")
ARGmed1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1106_1_ARG_med_u10.tif")
outstack <- raster::stack(ARGmax1103,NBRstd1103,ARGmed1106); names(outstack) <- c("ARGmax1103","NBRstd1103","ARGmed1106")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2011_2009bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2009 model with 2013 data: heb_2013 ~ ARGmax1303 + NBRstd1303 + ARGmed1306
rfmodel <- randomForest(heb_2013 ~ ARGmax1303 + NBRstd1303 + ARGmed1306, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1303_0_ARG_max_u10.tif")
NBRstd1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1303_0_NBR_std_u10.tif")
ARGmed1306 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1306_1_ARG_med_u10.tif")
outstack <- raster::stack(ARGmax1303,NBRstd1303,ARGmed1306); names(outstack) <- c("ARGmax1303","NBRstd1303","ARGmed1306")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2013_2009bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2009 model with 2015 data: heb_2015 ~ ARGmax1503 + NBRstd1503 + ARGmed1506
rfmodel <- randomForest(heb_2015 ~ ARGmax1503 + NBRstd1503 + ARGmed1506, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1503_0_ARG_max_u10.tif")
NBRstd1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1503_0_NBR_std_u10.tif")
ARGmed1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1506_1_ARG_med_u10.tif")
outstack <- raster::stack(ARGmax1503,NBRstd1503,ARGmed1506); names(outstack) <- c("ARGmax1503","NBRstd1503","ARGmed1506")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2015_2009bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2009 model with 2017 data: heb_2017 ~ ARGmax1703 + NBRstd1703 + ARGmed1706
rfmodel <- randomForest(heb_2017 ~ ARGmax1703 + NBRstd1703 + ARGmed1706, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1703_0_ARG_max_u10.tif")
NBRstd1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1703_0_NBR_std_u10.tif")
ARGmed1706 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1706_1_ARG_med_u10.tif")
outstack <- raster::stack(ARGmax1703,NBRstd1703,ARGmed1706); names(outstack) <- c("ARGmax1703","NBRstd1703","ARGmed1706")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2017_2009bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2009 model with 2018 data: Herb_2018 ~ ARGmax1803 + NBRstd1803 + ARGmed1806
rfmodel <- randomForest(Herb_2018 ~ ARGmax1803 + NBRstd1803 + ARGmed1806, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1803_0_ARG_max_u10.tif")
NBRstd1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1803_0_NBR_std_u10.tif")
ARGmed1806 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1806_1_ARG_med_u10.tif")
outstack <- raster::stack(ARGmax1803,NBRstd1803,ARGmed1806); names(outstack) <- c("ARGmax1803","NBRstd1803","ARGmed1806")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()   
writeRaster(outpredict, "rf2018_2009bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start


# U11 Predict 2011 model on all years ------------
# U11 Predict 2011 model with 2009 data: heb_2009 ~ NDVstd0903 + ARGmea0906 + B3std0906
rfmodel <- randomForest(heb_2009 ~ NDVstd0903 + ARGmea0906 + B3std0906, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd0903_0_NDV_std_u11.tif")
ARGmea0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea0906_1_ARG_mean_u11.tif")
B3std0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std0906_1_B3_std_u11.tif")
outstack <- raster::stack(NDVstd0903,ARGmea0906,B3std0906); names(outstack) <- c("NDVstd0903","ARGmea0906","B3std0906")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2009_2011bm_u11.tif", format='GTiff',overwrite=TRUE)

# U11 Predict 2011 model with 2011 data: heb_2011 ~ NDVstd1103 + ARGmea1106 + B3std1106
rfmodel <- randomForest(heb_2011 ~ NDVstd1103 + ARGmea1106 + B3std1106, 
                        data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1103_0_NDV_std_u11.tif")
ARGmea1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1106_1_ARG_mean_u11.tif")
B3std1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1106_1_B3_std_u11.tif")
outstack <- raster::stack(NDVstd1103,ARGmea1106,B3std1106); names(outstack) <- c("NDVstd1103","ARGmea1106","B3std1106")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2011_2011bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start


# U11 Predict 2011 model with 2013 data: heb_2011 ~ NDVstd1303 + ARGmea1306 + B3std1306
rfmodel <- randomForest(heb_2013 ~ NDVstd1303 + ARGmea1306 + B3std1306, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1303_0_NDV_std_u11.tif")
ARGmea1306 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1306_1_ARG_mean_u11.tif")
B3std1306 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1306_1_B3_std_u11.tif")
outstack <- raster::stack(NDVstd1303,ARGmea1306,B3std1306); names(outstack) <- c("NDVstd1303","ARGmea1306","B3std1306")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2013_2011bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U11 Predict 2009 model with 2015 data: heb_2015 ~ NDVstd1503 + ARGmea1506 + B3std1506
rfmodel <- randomForest(heb_2015 ~ NDVstd1503 + ARGmea1506 + B3std1506, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1503_0_NDV_std_u11.tif")
ARGmea1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1506_1_ARG_mean_u11.tif")
B3std1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1506_1_B3_std_u11.tif")
outstack <- raster::stack(NDVstd1503,ARGmea1506,B3std1506); names(outstack) <- c("NDVstd1503","ARGmea1506","B3std1506")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2015_2011bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U11 Predict 2011 model with 2017 data: heb_2017 ~ NDVstd1703 + ARGmea1706 + B3std1706
rfmodel <- randomForest(heb_2017 ~ NDVstd1703 + ARGmea1706 + B3std1706, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1703_0_NDV_std_u11.tif")
ARGmea1706 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1706_1_ARG_mean_u11.tif")
B3std1706 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1706_1_B3_std_u11.tif")
outstack <- raster::stack(NDVstd1703,ARGmea1706,B3std1706); names(outstack) <- c("NDVstd1703","ARGmea1706","B3std1706")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2017_2011bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U11 Predict 2011 model with 2018 data: Herb_2018 ~ NDVstd1803 + ARGmea1806 + B3std1806
rfmodel <- randomForest(Herb_2018 ~ NDVstd1803 + ARGmea1806 + B3std1806, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1803_0_NDV_std_u11.tif")
ARGmea1806 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1806_1_ARG_mean_u11.tif")
B3std1806 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1806_1_B3_std_u11.tif")
outstack <- raster::stack(NDVstd1803,ARGmea1806,B3std1806); names(outstack) <- c("NDVstd1803","ARGmea1806","B3std1806")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2018_2011bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2011 model on all years ------------
# U10 Predict 2011 model with 2009 data: heb_2009 ~ NDVstd0903 + ARGmea0906 + B3std0906
rfmodel <- randomForest(heb_2009 ~ NDVstd0903 + ARGmea0906 + B3std0906, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd0903_0_NDV_std_u10.tif")
ARGmea0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea0906_1_ARG_mean_u10.tif")
B3std0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std0906_1_B3_std_u10.tif")
outstack <- raster::stack(NDVstd0903,ARGmea0906,B3std0906); names(outstack) <- c("NDVstd0903","ARGmea0906","B3std0906")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2009_2011bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2011 model with 2011 data: heb_2011 ~ NDVstd1103 + ARGmea1106 + B3std1106
rfmodel <- randomForest(heb_2011 ~ NDVstd1103 + ARGmea1106 + B3std1106, 
                        data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1103_0_NDV_std_u10.tif")
ARGmea1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1106_1_ARG_mean_u10.tif")
B3std1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1106_1_B3_std_u10.tif")
outstack <- raster::stack(NDVstd1103,ARGmea1106,B3std1106); names(outstack) <- c("NDVstd1103","ARGmea1106","B3std1106")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2011_2011bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start


# U10 Predict 2011 model with 2013 data: heb_2011 ~ NDVstd1303 + ARGmea1306 + B3std1306
rfmodel <- randomForest(heb_2013 ~ NDVstd1303 + ARGmea1306 + B3std1306, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1303_0_NDV_std_u10.tif")
ARGmea1306 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1306_1_ARG_mean_u10.tif")
B3std1306 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1306_1_B3_std_u10.tif")
outstack <- raster::stack(NDVstd1303,ARGmea1306,B3std1306); names(outstack) <- c("NDVstd1303","ARGmea1306","B3std1306")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2013_2011bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2011 model with 2015 data: heb_2015 ~ NDVstd1503 + ARGmea1506 + B3std1506
rfmodel <- randomForest(heb_2015 ~ NDVstd1503 + ARGmea1506 + B3std1506, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1503_0_NDV_std_u10.tif")
ARGmea1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1506_1_ARG_mean_u10.tif")
B3std1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1506_1_B3_std_u10.tif")
outstack <- raster::stack(NDVstd1503,ARGmea1506,B3std1506); names(outstack) <- c("NDVstd1503","ARGmea1506","B3std1506")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2015_2011bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2011 model with 2017 data: heb_2017 ~ NDVstd1703 + ARGmea1706 + B3std1706
rfmodel <- randomForest(heb_2017 ~ NDVstd1703 + ARGmea1706 + B3std1706, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1703_0_NDV_std_u10.tif")
ARGmea1706 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1706_1_ARG_mean_u10.tif")
B3std1706 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1706_1_B3_std_u10.tif")
outstack <- raster::stack(NDVstd1703,ARGmea1706,B3std1706); names(outstack) <- c("NDVstd1703","ARGmea1706","B3std1706")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2017_2010bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2011 model with 2018 data: Herb_2018 ~ NDVstd1803 + ARGmea1806 + B3std1806
rfmodel <- randomForest(Herb_2018 ~ NDVstd1803 + ARGmea1806 + B3std1806, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1803_0_NDV_std_u10.tif")
ARGmea1806 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1806_1_ARG_mean_u10.tif")
B3std1806 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1806_1_B3_std_u10.tif")
outstack <- raster::stack(NDVstd1803,ARGmea1806,B3std1806); names(outstack) <- c("NDVstd1803","ARGmea1806","B3std1806")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2018_2011bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start



# U11 Predict 2013 model on all years ------------
# U11 Predict 2013 model with 2009 data: heb_2009 ~ ARGmea0906 + NDVstd0903 + B3max0909
rfmodel <- randomForest(heb_2009 ~ NDVstd0903 + ARGmea0906 + B3max0909, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmea0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea0906_1_ARG_mean_u11.tif")
NDVstd0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd0903_0_NDV_std_u11.tif")
B3max0909 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3max0909_2_B3_max_u11.tif")
outstack <- raster::stack(NDVstd0903,ARGmea0906,B3max0909); names(outstack) <- c("NDVstd0903","ARGmea0906","B3max0909")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2009_2013bm_u11.tif", format='GTiff',overwrite=TRUE)

# U11 Predict 2013 model with 2011 data: heb_2011 ~ ARGmea0906 + NDVstd0903 + B3max0909
rfmodel <- randomForest(heb_2011 ~ ARGmea1106 + NDVstd1103 + B3max1109, 
                        data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmea1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1106_1_ARG_mean_u11.tif")
NDVstd1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1103_0_NDV_std_u11.tif")
B3max1109 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3max1109_2_B3_max_u11.tif")
outstack <- raster::stack(ARGmea1106,NDVstd1103,B3max1109); names(outstack) <- c("ARGmea1106","NDVstd1103","B3max1109")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2011_2013bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start


# U11 Predict 2013 model with 2013 data: heb_2013 ~ ARGmea1306 + NDVstd1303 + B3max1309
rfmodel <- randomForest(heb_2013 ~ ARGmea1306 + NDVstd1303 + B3max1309, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1303_0_NDV_std_u11.tif")
ARGmea1306 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1306_1_ARG_mean_u11.tif")
B3max1309 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3max1309_2_B3_max_u11.tif")
outstack <- raster::stack(NDVstd1303,ARGmea1306,B3max1309); names(outstack) <- c("NDVstd1303","ARGmea1306","B3max1309")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2013_2013bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U11 Predict 2013 model with 2015 data: heb_2015 ~ ARGmea1506 + NDVstd1503 + B3max1509
rfmodel <- randomForest(heb_2015 ~ ARGmea1506 + NDVstd1503 + B3max1509, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1503_0_NDV_std_u11.tif")
ARGmea1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1506_1_ARG_mean_u11.tif")
B3max1509 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3max1509_2_B3_max_u11.tif")
outstack <- raster::stack(NDVstd1503,ARGmea1506,B3max1509); names(outstack) <- c("NDVstd1503","ARGmea1506","B3max1509")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2015_2013bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U11 Predict 2013 model with 2017 data: heb_2017 ~ ARGmea1706 + NDVstd1703 + B3max1709
rfmodel <- randomForest(heb_2017 ~ ARGmea1706 + NDVstd1703 + B3max1709, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1703_0_NDV_std_u11.tif")
ARGmea1706 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1706_1_ARG_mean_u11.tif")
B3max1709 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3max1709_2_B3_max_u11.tif")
outstack <- raster::stack(NDVstd1703,ARGmea1706,B3max1709); names(outstack) <- c("NDVstd1703","ARGmea1706","B3max1709")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2017_2013bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U11 Predict 2013 model with 2018 data: Herb_2018 ~ ARGmea0906 + NDVstd0903 + B3max1809
rfmodel <- randomForest(Herb_2018 ~ NDVstd1803 + ARGmea1806 + B3max1809, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1803_0_NDV_std_u11.tif")
ARGmea1806 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1806_1_ARG_mean_u11.tif")
B3max1809 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3max0909_2_B3_max_u11.tif")
outstack <- raster::stack(NDVstd1803,ARGmea1806,B3max1809); names(outstack) <- c("NDVstd1803","ARGmea1806","B3max1809")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2018_2013bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2013 model on all years ------------
# U10 Predict 2013 model with 2009 data: heb_2009 ~ ARGmea0906 + NDVstd0903 + B3max0909
rfmodel <- randomForest(heb_2009 ~ NDVstd0903 + ARGmea0906 + B3max0909, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmea0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea0906_1_ARG_mean_u10.tif")
NDVstd0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd0903_0_NDV_std_u10.tif")
B3max0909 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3max0909_2_B3_max_u10.tif")
outstack <- raster::stack(NDVstd0903,ARGmea0906,B3max0909); names(outstack) <- c("NDVstd0903","ARGmea0906","B3max0909")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2009_2013bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2013 model with 2011 data: heb_2011 ~ ARGmea0906 + NDVstd0903 + B3max0909
rfmodel <- randomForest(heb_2011 ~ ARGmea1106 + NDVstd1103 + B3max1109, 
                        data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmea1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1106_1_ARG_mean_u10.tif")
NDVstd1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1103_0_NDV_std_u10.tif")
B3max1109 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3max1109_2_B3_max_u10.tif")
outstack <- raster::stack(ARGmea1106,NDVstd1103,B3max1109); names(outstack) <- c("ARGmea1106","NDVstd1103","B3max1109")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2011_2013bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start


# U10 Predict 2013 model with 2013 data: heb_2013 ~ ARGmea1306 + NDVstd1303 + B3max1309
rfmodel <- randomForest(heb_2013 ~ ARGmea1306 + NDVstd1303 + B3max1309, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1303_0_NDV_std_u10.tif")
ARGmea1306 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1306_1_ARG_mean_u10.tif")
B3max1309 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3max1309_2_B3_max_u10.tif")
outstack <- raster::stack(NDVstd1303,ARGmea1306,B3max1309); names(outstack) <- c("NDVstd1303","ARGmea1306","B3max1309")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2013_2013bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2013 model with 2015 data: heb_2015 ~ ARGmea1506 + NDVstd1503 + B3max1509
rfmodel <- randomForest(heb_2015 ~ ARGmea1506 + NDVstd1503 + B3max1509, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1503_0_NDV_std_u10.tif")
ARGmea1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1506_1_ARG_mean_u10.tif")
B3max1509 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3max1509_2_B3_max_u10.tif")
outstack <- raster::stack(NDVstd1503,ARGmea1506,B3max1509); names(outstack) <- c("NDVstd1503","ARGmea1506","B3max1509")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2015_2013bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2013 model with 2017 data: heb_2017 ~ ARGmea1706 + NDVstd1703 + B3max1709
rfmodel <- randomForest(heb_2017 ~ ARGmea1706 + NDVstd1703 + B3max1709, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1703_0_NDV_std_u10.tif")
ARGmea1706 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1706_1_ARG_mean_u10.tif")
B3max1709 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3max1709_2_B3_max_u10.tif")
outstack <- raster::stack(NDVstd1703,ARGmea1706,B3max1709); names(outstack) <- c("NDVstd1703","ARGmea1706","B3max1709")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2017_2013bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2013 model with 2018 data: Herb_2018 ~ ARGmea0906 + NDVstd0903 + B3max1809
rfmodel <- randomForest(Herb_2018 ~ NDVstd1803 + ARGmea1806 + B3max1809, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
NDVstd1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVstd1803_0_NDV_std_u10.tif")
ARGmea1806 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmea1806_1_ARG_mean_u10.tif")
B3max1809 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3max0909_2_B3_max_u10.tif")
outstack <- raster::stack(NDVstd1803,ARGmea1806,B3max1809); names(outstack) <- c("NDVstd1803","ARGmea1806","B3max1809")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2018_2013bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start


# U11 Predict 2015 model on all years ------------
# U11 Predict 2015 model with 2009 data: heb_2009 ~  ARGmax0903 + ARGmax0906 + NDVmin0909
rfmodel <- randomForest(heb_2009 ~ ARGmax0903 + ARGmax0906 + NDVmin0909, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax0903_0_ARG_max_u11.tif")
ARGmax0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax0906_1_ARG_max_u11.tif")
NDVmin0909 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmin0909_2_NDV_min_u11.tif")
outstack <- raster::stack(ARGmax0903,ARGmax0906,NDVmin0909); names(outstack) <- c("ARGmax0903","ARGmax0906","NDVmin0909")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2009_2015bm_u11.tif", format='GTiff',overwrite=TRUE)

# U11 Predict 2015 model with 2011 data: heb_2011 ~ ARGmax1103 + ARGmax1106 + NDVmin1109
rfmodel <- randomForest(heb_2011 ~ ARGmax1103 + ARGmax1106 + NDVmin1109, 
                        data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1103_0_ARG_max_u11.tif")
ARGmax1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1106_1_ARG_max_u11.tif")
NDVmin1109 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmin1109_2_NDV_min_u11.tif")
outstack <- raster::stack(ARGmax1103,ARGmax1106,NDVmin1109); names(outstack) <- c("ARGmax1103","ARGmax1106","NDVmin1109")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2011_2015bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start


# U11 Predict 2015 model with 2013 data: heb_2013 ~ ARGmax1303 + ARGmax1306 + NDVmin1309
rfmodel <- randomForest(heb_2013 ~ ARGmax1303 + ARGmax1306 + NDVmin1309, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1303_0_ARG_max_u11.tif")
ARGmax1306 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1306_1_ARG_max_u11.tif")
NDVmin1309 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmin1309_2_NDV_min_u11.tif")
outstack <- raster::stack(ARGmax1303,ARGmax1306,NDVmin1309); names(outstack) <- c("ARGmax1303","ARGmax1306","NDVmin1309")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2013_2015bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U11 Predict 2015 model with 2015 data: heb_2015 ~ ARGmax1503 + ARGmax1506 + NDVmin1509
rfmodel <- randomForest(heb_2015 ~ ARGmax1503 + ARGmax1506 + NDVmin1509, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1503_0_ARG_max_u11.tif")
ARGmax1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1506_1_ARG_max_u11.tif")
NDVmin1509 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmin1509_2_NDV_min_u11.tif")
outstack <- raster::stack(ARGmax1503,ARGmax1506,NDVmin1509); names(outstack) <- c("ARGmax1503","ARGmax1506","NDVmin1509")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2015_2015bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U11 Predict 2015 model with 2017 data: heb_2017 ~  ARGmax1703 + ARGmax1706 + NDVmin1709
rfmodel <- randomForest(heb_2017 ~ ARGmax1703 + ARGmax1706 + NDVmin1709, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1703_0_ARG_max_u11.tif")
ARGmax1706 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1706_1_ARG_max_u11.tif")
NDVmin1709 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmin1709_2_NDV_min_u11.tif")
outstack <- raster::stack(ARGmax1703,ARGmax1706,NDVmin1709); names(outstack) <- c("ARGmax1703","ARGmax1706","NDVmin1709")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2017_2015bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U11 Predict 2013 model with 2018 data: Herb_2018 ~ ARGmax1803 + ARGmax1806 + NDVmin1809
rfmodel <- randomForest(Herb_2018 ~ ARGmax1803 + ARGmax1806 + NDVmin1809, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1803_0_ARG_max_u11.tif")
ARGmax1806 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1806_1_ARG_max_u11.tif")
NDVmin1809 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmin1809_2_NDV_min_u11.tif")
outstack <- raster::stack(ARGmax1803,ARGmax1806,NDVmin1809); names(outstack) <- c("ARGmax1803","ARGmax1806","NDVmin1809")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2018_2015bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start




# U10 Predict 2015 model on all years ------------
# U10 Predict 2015 model with 2009 data: heb_2009 ~  ARGmax0903 + ARGmax0906 + NDVmin0909
rfmodel <- randomForest(heb_2009 ~ ARGmax0903 + ARGmax0906 + NDVmin0909, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax0903_0_ARG_max_u10.tif")
ARGmax0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax0906_1_ARG_max_u10.tif")
NDVmin0909 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmin0909_2_NDV_min_u10.tif")
outstack <- raster::stack(ARGmax0903,ARGmax0906,NDVmin0909); names(outstack) <- c("ARGmax0903","ARGmax0906","NDVmin0909")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2009_2015bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2015 model with 2011 data: heb_2011 ~ ARGmax1103 + ARGmax1106 + NDVmin1109
rfmodel <- randomForest(heb_2011 ~ ARGmax1103 + ARGmax1106 + NDVmin1109, 
                        data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1103_0_ARG_max_u10.tif")
ARGmax1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1106_1_ARG_max_u10.tif")
NDVmin1109 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmin1109_2_NDV_min_u10.tif")
outstack <- raster::stack(ARGmax1103,ARGmax1106,NDVmin1109); names(outstack) <- c("ARGmax1103","ARGmax1106","NDVmin1109")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2011_2015bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start


# U10 Predict 2015 model with 2013 data: heb_2013 ~ ARGmax1303 + ARGmax1306 + NDVmin1309
rfmodel <- randomForest(heb_2013 ~ ARGmax1303 + ARGmax1306 + NDVmin1309, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1303_0_ARG_max_u10.tif")
ARGmax1306 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1306_1_ARG_max_u10.tif")
NDVmin1309 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmin1309_2_NDV_min_u10.tif")
outstack <- raster::stack(ARGmax1303,ARGmax1306,NDVmin1309); names(outstack) <- c("ARGmax1303","ARGmax1306","NDVmin1309")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2013_2015bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2015 model with 2015 data: heb_2015 ~ ARGmax1503 + ARGmax1506 + NDVmin1509
rfmodel <- randomForest(heb_2015 ~ ARGmax1503 + ARGmax1506 + NDVmin1509, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1503_0_ARG_max_u10.tif")
ARGmax1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1506_1_ARG_max_u10.tif")
NDVmin1509 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmin1509_2_NDV_min_u10.tif")
outstack <- raster::stack(ARGmax1503,ARGmax1506,NDVmin1509); names(outstack) <- c("ARGmax1503","ARGmax1506","NDVmin1509")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2015_2015bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2015 model with 2017 data: heb_2017 ~  ARGmax1703 + ARGmax1706 + NDVmin1709
rfmodel <- randomForest(heb_2017 ~ ARGmax1703 + ARGmax1706 + NDVmin1709, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1703_0_ARG_max_u10.tif")
ARGmax1706 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1706_1_ARG_max_u10.tif")
NDVmin1709 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmin1709_2_NDV_min_u10.tif")
outstack <- raster::stack(ARGmax1703,ARGmax1706,NDVmin1709); names(outstack) <- c("ARGmax1703","ARGmax1706","NDVmin1709")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2017_2015bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2013 model with 2018 data: Herb_2018 ~ ARGmax1803 + ARGmax1806 + NDVmin1809
rfmodel <- randomForest(Herb_2018 ~ ARGmax1803 + ARGmax1806 + NDVmin1809, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1803_0_ARG_max_u10.tif")
ARGmax1806 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1806_1_ARG_max_u10.tif")
NDVmin1809 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NDVmin1809_2_NDV_min_u10.tif")
outstack <- raster::stack(ARGmax1803,ARGmax1806,NDVmin1809); names(outstack) <- c("ARGmax1803","ARGmax1806","NDVmin1809")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2018_2015bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start




# U11 Predict 2017 model on all years ------------
# U11 Predict 2017 model with 2009 data: heb_2009 ~   ARGmax0906 + NBRstd0903 + ARGmed0909
rfmodel <- randomForest(heb_2009 ~ ARGmax0906 + NBRstd0903 + ARGmed0909, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax0906_1_ARG_max_u11.tif")
NBRstd0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd0903_0_NBR_std_u11.tif")
ARGmed0909 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed0909_2_ARG_med_u11.tif")
outstack <- raster::stack(NBRstd0903,ARGmax0906,ARGmed0909); names(outstack) <- c("NBRstd0903","ARGmax0906","ARGmed0909")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2009_2017bm_u11.tif", format='GTiff',overwrite=TRUE)

# U11 Predict 2017 model with 2011 data: heb_2011 ~ ARGmax1106 + NBRstd1103 + ARGmed1109
rfmodel <- randomForest(heb_2011 ~ ARGmax1106 + NBRstd1103 + ARGmed1109, 
                        data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1106_1_ARG_max_u11.tif")
NBRstd1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1103_0_NBR_std_u11.tif")
ARGmed1109 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1109_2_ARG_med_u11.tif")
outstack <- raster::stack(NBRstd1103,ARGmax1106,ARGmed1109); names(outstack) <- c("NBRstd1103","ARGmax1106","ARGmed1109")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2011_2017bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start


# U11 Predict 2017 model with 2013 data: heb_2013 ~ ARGmax1306 + NBRstd1303 + ARGmed1309
rfmodel <- randomForest(heb_2013 ~ ARGmax1306 + NBRstd1303 + ARGmed1309, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1306 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1306_1_ARG_max_u11.tif")
NBRstd1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1303_0_NBR_std_u11.tif")
ARGmed1309 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1309_2_ARG_med_u11.tif")
outstack <- raster::stack(NBRstd1303,ARGmax1306,ARGmed1309); names(outstack) <- c("NBRstd1303","ARGmax1306","ARGmed1309")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2013_2017bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start


# U11 Predict 2017 model with 2015 data: heb_2015 ~ ARGmax1506 + NBRstd1503 + ARGmed1509
rfmodel <- randomForest(heb_2015 ~ ARGmax1506 + NBRstd1503 + ARGmed1509, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1506_1_ARG_max_u11.tif")
NBRstd1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1503_0_NBR_std_u11.tif")
ARGmed1509 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1509_2_ARG_med_u11.tif")
outstack <- raster::stack(NBRstd1503,ARGmax1506,ARGmed1509); names(outstack) <- c("NBRstd1503","ARGmax1506","ARGmed1509")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2015_2017bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U11 Predict 2017 model with 2017 data: heb_2017 ~  ARGmax1706 + NBRstd1703 + ARGmed1709
rfmodel <- randomForest(heb_2017 ~ ARGmax1706 + NBRstd1703 + ARGmed1709, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1706 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1706_1_ARG_max_u11.tif")
NBRstd1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1703_0_NBR_std_u11.tif")
ARGmed1709 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1709_2_ARG_med_u11.tif")
outstack <- raster::stack(NBRstd1703,ARGmax1706,ARGmed1709); names(outstack) <- c("NBRstd1703","ARGmax1706","ARGmed1709")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2017_2017bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U11 Predict 2017 model with 2018 data: Herb_2018 ~ ARGmax1806 + NBRstd1803 + ARGmed1809
rfmodel <- randomForest(Herb_2018 ~ ARGmax1806 + NBRstd1803 + ARGmed1809, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1806 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1806_1_ARG_max_u11.tif")
NBRstd1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1803_0_NBR_std_u11.tif")
ARGmed1809 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1809_2_ARG_med_u11.tif")
outstack <- raster::stack(NBRstd1803,ARGmax1806,ARGmed1809); names(outstack) <- c("NBRstd1803","ARGmax1806","ARGmed1809")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2018_2017bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start





# U10 Predict 2017 model on all years ------------
# U10 Predict 2017 model with 2009 data: heb_2009 ~   ARGmax0906 + NBRstd0903 + ARGmed0909
rfmodel <- randomForest(heb_2009 ~ ARGmax0906 + NBRstd0903 + ARGmed0909, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax0906_1_ARG_max_u10.tif")
NBRstd0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd0903_0_NBR_std_u10.tif")
ARGmed0909 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed0909_2_ARG_med_u10.tif")
outstack <- raster::stack(NBRstd0903,ARGmax0906,ARGmed0909); names(outstack) <- c("NBRstd0903","ARGmax0906","ARGmed0909")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster();    finish = Sys.time(); finish-start
writeRaster(outpredict, "rf2009_2017bm_u10.tif", format='GTiff',overwrite=TRUE)

# U10 Predict 2017 model with 2011 data: heb_2011 ~ ARGmax1106 + NBRstd1103 + ARGmed1109
rfmodel <- randomForest(heb_2011 ~ ARGmax1106 + NBRstd1103 + ARGmed1109, 
                        data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1106_1_ARG_max_u10.tif")
NBRstd1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1103_0_NBR_std_u10.tif")
ARGmed1109 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1109_2_ARG_med_u10.tif")
outstack <- raster::stack(NBRstd1103,ARGmax1106,ARGmed1109); names(outstack) <- c("NBRstd1103","ARGmax1106","ARGmed1109")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2011_2017bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2017 model with 2013 data: heb_2013 ~ ARGmax1306 + NBRstd1303 + ARGmed1309
rfmodel <- randomForest(heb_2013 ~ ARGmax1306 + NBRstd1303 + ARGmed1309, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1306 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1306_1_ARG_max_u10.tif")
NBRstd1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1303_0_NBR_std_u10.tif")
ARGmed1309 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1309_2_ARG_med_u10.tif")
outstack <- raster::stack(NBRstd1303,ARGmax1306,ARGmed1309); names(outstack) <- c("NBRstd1303","ARGmax1306","ARGmed1309")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2013_2017bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U11 Predict 2017 model with 2015 data: heb_2015 ~ ARGmax1506 + NBRstd1503 + ARGmed1509
rfmodel <- randomForest(heb_2015 ~ ARGmax1506 + NBRstd1503 + ARGmed1509, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1506_1_ARG_max_u10.tif")
NBRstd1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1503_0_NBR_std_u10.tif")
ARGmed1509 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1509_2_ARG_med_u10.tif")
outstack <- raster::stack(NBRstd1503,ARGmax1506,ARGmed1509); names(outstack) <- c("NBRstd1503","ARGmax1506","ARGmed1509")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2015_2017bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U11 Predict 2017 model with 2017 data: heb_2017 ~  ARGmax1706 + NBRstd1703 + ARGmed1709
rfmodel <- randomForest(heb_2017 ~ ARGmax1706 + NBRstd1703 + ARGmed1709, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1706 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1706_1_ARG_max_u10.tif")
NBRstd1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1703_0_NBR_std_u10.tif")
ARGmed1709 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1709_2_ARG_med_u10.tif")
outstack <- raster::stack(NBRstd1703,ARGmax1706,ARGmed1709); names(outstack) <- c("NBRstd1703","ARGmax1706","ARGmed1709")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2017_2017bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U11 Predict 2017 model with 2018 data: Herb_2018 ~ ARGmax1806 + NBRstd1803 + ARGmed1809
rfmodel <- randomForest(Herb_2018 ~ ARGmax1806 + NBRstd1803 + ARGmed1809, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1806 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1806_1_ARG_max_u10.tif")
NBRstd1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1803_0_NBR_std_u10.tif")
ARGmed1809 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmed1809_2_ARG_med_u10.tif")
outstack <- raster::stack(NBRstd1803,ARGmax1806,ARGmed1809); names(outstack) <- c("NBRstd1803","ARGmax1806","ARGmed1809")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2018_2017bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start


# U11 Predict 2018 model on all years ------------
# U11 Predict 2018 model with 2009 data: heb_2009 ~  ARGmax0906 + NBRstd1803 +  B3std1803
rfmodel <- randomForest(heb_2009 ~ ARGmax0906 + NBRstd0903 +  B3std0903, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax0906_1_ARG_max_u11.tif")
NBRstd0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd0903_0_NBR_std_u11.tif")
B3std0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std0903_0_B3_std_u11.tif")
outstack <- raster::stack(NBRstd0903,ARGmax0906,B3std0903); names(outstack) <- c("NBRstd0903","ARGmax0906","B3std0903")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2009_2018bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U11 Predict 2018 model with 2011 data: heb_2011 ~ ARGmax1106 + NBRstd1103 + B3std1103
rfmodel <- randomForest(heb_2011 ~ ARGmax1106 + NBRstd1103 + B3std1103, 
                        data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1106_1_ARG_max_u11.tif")
NBRstd1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1103_0_NBR_std_u11.tif")
B3std1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1103_0_B3_std_u11.tif")
outstack <- raster::stack(NBRstd1103,ARGmax1106,B3std1103); names(outstack) <- c("NBRstd1103","ARGmax1106","B3std1103")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2011_2018bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start


# U11 Predict 2018 model with 2013 data: heb_2013 ~ ARGmax1306 + NBRstd1303 + B3std1303
rfmodel <- randomForest(heb_2013 ~ ARGmax1306 + NBRstd1303 + B3std1303, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1306 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1306_1_ARG_max_u11.tif")
NBRstd1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1303_0_NBR_std_u11.tif")
B3std1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1303_0_B3_std_u11.tif")
outstack <- raster::stack(NBRstd1303,ARGmax1306,B3std1303); names(outstack) <- c("NBRstd1303","ARGmax1306","B3std1303")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2013_2018bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U11 Predict 2018 model with 2015 data: heb_2015 ~ ARGmax1506 + NBRstd1503 + B3std1503
rfmodel <- randomForest(heb_2015 ~ ARGmax1506 + NBRstd1503 + ARGmed1509, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1506_1_ARG_max_u11.tif")
NBRstd1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1503_0_NBR_std_u11.tif")
B3std1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1503_0_B3_std_u11.tif")
outstack <- raster::stack(NBRstd1503,ARGmax1506,B3std1503); names(outstack) <- c("NBRstd1503","ARGmax1506","B3std1503")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2015_2018bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U11 Predict 2018 model with 2017 data: heb_2017 ~  ARGmax1706 + NBRstd1703 + B3std1703
rfmodel <- randomForest(heb_2017 ~ ARGmax1706 + NBRstd1703 + B3std1703, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1706 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1706_1_ARG_max_u11.tif")
NBRstd1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1703_0_NBR_std_u11.tif")
B3std1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1703_0_B3_std_u11.tif")
outstack <- raster::stack(NBRstd1703,ARGmax1706,B3std1703); names(outstack) <- c("NBRstd1703","ARGmax1706","B3std1703")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2017_2018bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U11 Predict 2018 model with 2018 data: Herb_2018 ~ ARGmax1806 + NBRstd1803 + B3std1803
rfmodel <- randomForest(Herb_2018 ~ ARGmax1806 + NBRstd1803 + ARGmed1809, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1806 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1806_1_ARG_max_u11.tif")
NBRstd1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1803_0_NBR_std_u11.tif")
B3std1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1803_0_B3_std_u11.tif")
outstack <- raster::stack(NBRstd1803,ARGmax1806,B3std1803); names(outstack) <- c("NBRstd1803","ARGmax1806","B3std1803")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2018_2018bm_u11.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2018 model on all years ------------
# U10 Predict 2018 model with 2009 data: heb_2009 ~  ARGmax0906 + NBRstd1803 +  B3std1803
rfmodel <- randomForest(heb_2009 ~ ARGmax0906 + NBRstd0903 +  B3std0903, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax0906 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax0906_1_ARG_max_u10.tif")
NBRstd0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd0903_0_NBR_std_u10.tif")
B3std0903 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std0903_0_B3_std_u10.tif")
outstack <- raster::stack(NBRstd0903,ARGmax0906,B3std0903); names(outstack) <- c("NBRstd0903","ARGmax0906","B3std0903")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2009_2018bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2018 model with 2011 data: heb_2011 ~ ARGmax1106 + NBRstd1103 + B3std1103
rfmodel <- randomForest(heb_2011 ~ ARGmax1106 + NBRstd1103 + B3std1103, 
                        data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1106 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1106_1_ARG_max_u10.tif")
NBRstd1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1103_0_NBR_std_u10.tif")
B3std1103 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1103_0_B3_std_u10.tif")
outstack <- raster::stack(NBRstd1103,ARGmax1106,B3std1103); names(outstack) <- c("NBRstd1103","ARGmax1106","B3std1103")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2011_2018bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2018 model with 2013 data: heb_2013 ~ ARGmax1306 + NBRstd1303 + B3std1303
rfmodel <- randomForest(heb_2013 ~ ARGmax1306 + NBRstd1303 + B3std1303, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1306 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1306_1_ARG_max_u10.tif")
NBRstd1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1303_0_NBR_std_u10.tif")
B3std1303 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1303_0_B3_std_u10.tif")
outstack <- raster::stack(NBRstd1303,ARGmax1306,B3std1303); names(outstack) <- c("NBRstd1303","ARGmax1306","B3std1303")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2013_2018bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2018 model with 2015 data: heb_2015 ~ ARGmax1506 + NBRstd1503 + B3std1503
rfmodel <- randomForest(heb_2015 ~ ARGmax1506 + NBRstd1503 + B3std1503, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1506 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1506_1_ARG_max_u10.tif")
NBRstd1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1503_0_NBR_std_u10.tif")
B3std1503 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1503_0_B3_std_u10.tif")
outstack <- raster::stack(NBRstd1503,ARGmax1506,B3std1503); names(outstack) <- c("NBRstd1503","ARGmax1506","B3std1503")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2015_2018bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2018 model with 2017 data: heb_2017 ~  ARGmax1706 + NBRstd1703 + B3std1703
rfmodel <- randomForest(heb_2017 ~ ARGmax1706 + NBRstd1703 + B3std1703, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1706 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1706_1_ARG_max_u10.tif")
NBRstd1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1703_0_NBR_std_u10.tif")
B3std1703 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1703_0_B3_std_u10.tif")
outstack <- raster::stack(NBRstd1703,ARGmax1706,B3std1703); names(outstack) <- c("NBRstd1703","ARGmax1706","B3std1703")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2017_2018bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# U10 Predict 2018 model with 2018 data: Herb_2018 ~ ARGmax1806 + NBRstd1803 + B3std1803
rfmodel <- randomForest(Herb_2018 ~ ARGmax1806 + NBRstd1803 + B3std1803, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
ARGmax1806 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\ARGmax1806_1_ARG_max_u10.tif")
NBRstd1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\NBRstd1803_0_NBR_std_u10.tif")
B3std1803 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\fromEE\\B3std1803_0_B3_std_u10.tif")
outstack <- raster::stack(NBRstd1803,ARGmax1806,B3std1803); names(outstack) <- c("NBRstd1803","ARGmax1806","B3std1803")
start = Sys.time();beginCluster(12)
outpredict <- clusterR(outstack, raster::predict, args=list(rfmodel), progress='text', type='prob')
endCluster()
writeRaster(outpredict, "rf2018_2018bm_u10.tif", format='GTiff',overwrite=TRUE)
finish = Sys.time(); finish-start

# ------------------------------------------------------------------------------------------------------------
# OLDFIRES ##### Generate residuals for model yr 2009 for all years & save tables in cvtables_oldfires #####-------
require(broom)
# OLDFIRES - Generate residuals for model yr 2009 with 2009 data
rfts <- randomForest(heb_2009 ~ ARGmax0903 + NBRstd0903 + ARGmed0906, 
                     data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y09_nofires$pred_m09d09 <- predict(rfts)
y09_nofires$residm09d09 <- y09_nofires$pred_m09d09 - y09_nofires$heb_2009
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2009 # data year
cvtabof[1,2] <- 2009 # model year
cvtabof[1,3] <- mean(y09_nofires$residm09d09[y09$timesince > 199])
cvtabof[1,4] <- mean(y09_nofires$residm09d09[y09$timesince < 200])
# lm is based on all records & ts for that year, rather than the _nofires version.
d <- lm(y09_nofires$residm09d09 ~ y09_nofires$ts2009) # lm for all
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y09_nofires$residm09d09[y09_nofires$timesince < 200] ~ y09_nofires$timesince[y09_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m09d09_mresid_oldfires.csv"
write.csv(cvtabof,tname)


# OLDFIRES - Generate residuals for model year 2009 with 2011 data
rfts <- randomForest(heb_2011 ~ ARGmax1103 + NBRstd1103 + ARGmed1106, 
                     data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
# Generate residuals
y11_nofires$pred_m09d11 <- predict(rfts)
y11_nofires$residm09d11 <- y11_nofires$pred_m09d11 - y11_nofires$heb_2011
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2011 # data year
cvtabof[1,2] <- 2009 # model year
cvtabof[1,3] <- mean(y11_nofires$residm09d11[y11_nofires$timesince > 199])
cvtabof[1,4] <- mean(y11_nofires$residm09d11[y11_nofires$timesince < 200])
d <- lm(y11_nofires$residm09d11 ~ y11_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y11_nofires$residm09d11[y11_nofires$timesince < 200] ~ y11_nofires$timesince[y11_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m09d11_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate residuals for model year 2009 with 2013 data
rfts <- randomForest(heb_2013 ~ ARGmax1303 + NBRstd1303 + ARGmed1306, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y13_nofires$pred_m09d13 <- predict(rfts)
y13_nofires$residm09d13 <- y13_nofires$pred_m09d13 - y13_nofires$heb_2013
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 999, old fires < 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2013 # data year
cvtabof[1,2] <- 2009 # model year
cvtabof[1,3] <- mean(y13_nofires$residm09d13[y13_nofires$timesince > 199])
cvtabof[1,4] <- mean(y13_nofires$residm09d13[y13_nofires$timesince < 200])
d <- lm(y13_nofires$residm09d13 ~ y13_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y13_nofires$residm09d13[y13_nofires$timesince < 200] ~ y13_nofires$timesince[y13_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m09d13_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate residuals for model year 2009 with 2015 data
rfts <- randomForest(heb_2015 ~ ARGmax1503 + NBRstd1503 + ARGmed1506, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y15_nofires$pred_m09d15 <- predict(rfts)
y15_nofires$residm09d15 <- y15_nofires$pred_m09d15 - y15_nofires$heb_2015
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 999, old fires < 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2015 # data year
cvtabof[1,2] <- 2009 # model year
cvtabof[1,3] <- mean(y15_nofires$residm09d15[y15_nofires$timesince > 199])
cvtabof[1,4] <- mean(y15_nofires$residm09d15[y15_nofires$timesince < 200])
d <- lm(y15_nofires$residm09d15 ~ y15_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y15_nofires$residm09d15[y15_nofires$timesince < 200] ~ y15_nofires$timesince[y15_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m09d15_mresid_oldfires.csv"
write.csv(cvtabof,tname)


#  OLDFIRES - Generate residuals for model year 2009 with 2017 data
rfts <- randomForest(heb_2017 ~ ARGmax1703 + NBRstd1703 + ARGmed1706, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y17_nofires$pred_m09d17 <- predict(rfts)
y17_nofires$residm09d17 <- y17_nofires$pred_m09d17 - y17_nofires$heb_2017
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 999, old fires < 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2017 # data year
cvtabof[1,2] <- 2009 # model year
cvtabof[1,3] <- mean(y17_nofires$residm09d17[y17_nofires$timesince > 199]) # Select where we have no fire history
cvtabof[1,4] <- mean(y17_nofires$residm09d17[y17_nofires$timesince < 200]) # select where there are previous fires
d <- lm(y17_nofires$residm09d17 ~ y17_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y17_nofires$residm09d17[y17_nofires$timesince < 200] ~ y17_nofires$timesince[y17_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m09d17_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate residuals for model year 2009 with 2018 data
rfts <- randomForest(Herb_2018 ~ ARGmax1803 + NBRstd1803 + ARGmed1806, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y18_nofires$pred_m09d18 <- predict(rfts)
y18_nofires$residm09d18 <- y18_nofires$pred_m09d18 - y18_nofires$Herb_2018
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 999, old fires < 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2018 # data year
cvtabof[1,2] <- 2009 # model year
cvtabof[1,3] <- mean(y18_nofires$residm09d18[y18_nofires$timesince > 199])
cvtabof[1,4] <- mean(y18_nofires$residm09d18[y18_nofires$timesince < 200])
d <- lm(y18_nofires$residm09d18 ~ y18_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y18_nofires$residm09d18[y18_nofires$timesince < 200] ~ y18_nofires$timesince[y18_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m09d18_mresid_oldfires.csv"
write.csv(cvtabof,tname)


# OLDFIRES ##### Generate residuals for model yr 2011 for all years & save tables in cvtables_oldfires #####-------
# OLDFIRES - Generate 2011 model with 2009 data
rfts <- randomForest(heb_2009 ~ NDVstd0903 + ARGmea0906 + B3std0906, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y09_nofires$pred_m11d09 <- predict(rfts)
y09_nofires$residm11d09 <- y09_nofires$pred_m11d09 - y09_nofires$heb_2009
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2009 # data year
cvtabof[1,2] <- 2011 # model year
cvtabof[1,3] <- mean(y09_nofires$residm11d09[y09_nofires$timesince > 199])
cvtabof[1,4] <- mean(y09_nofires$residm11d09[y09_nofires$timesince < 200])
#plot(y09_nofires$residm11d09,y09_nofires$timesince)
d <- lm(y09_nofires$residm11d09 ~ y09_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y09_nofires$residm11d09[y09_nofires$timesince < 200] ~ y09_nofires$timesince[y09_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m11d09_mresid_oldfires.csv"
write.csv(cvtabof,tname)


# OLDFIRES - Generate 2011 model with 2011 data
rfts <- randomForest(heb_2011 ~ NDVstd1103 + ARGmea1106 + B3std1106, 
                        data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y11_nofires$pred_m11d11 <- predict(rfts)
y11_nofires$residm11d11 <- y11_nofires$pred_m11d11 - y11_nofires$heb_2011
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2011 # data year
cvtabof[1,2] <- 2011 # model year
cvtabof[1,3] <- mean(y11_nofires$residm11d11[y11_nofires$timesince > 199])
cvtabof[1,4] <- mean(y11_nofires$residm11d11[y11_nofires$timesince < 200])
d <- lm(y11_nofires$residm11d11 ~ y11_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y11_nofires$residm11d11[y11_nofires$timesince < 200] ~ y11_nofires$timesince[y11_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m11d11_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate 2011 model with 2013 data
rfts <- randomForest(heb_2013 ~ NDVstd1303 + ARGmea1306 + B3std1306, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y13_nofires$pred_m11d13 <- predict(rfts)
y13_nofires$residm11d13 <- y13_nofires$pred_m11d13 - y13_nofires$heb_2013
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2013 # data year
cvtabof[1,2] <- 2011 # model year
cvtabof[1,3] <- mean(y13_nofires$residm11d13[y13_nofires$timesince > 199]) 
cvtabof[1,4] <- mean(y13_nofires$residm11d13[y13_nofires$timesince < 200])
d <- lm(y13_nofires$residm11d13 ~ y13_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y13_nofires$residm11d13[y13_nofires$timesince < 200] ~ y13_nofires$timesince[y13_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m11d13_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate 2011 model with 2015 data
rfts <- randomForest(heb_2015 ~ NDVstd1503 + ARGmea1506 + B3std1506, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y15_nofires$pred_m11d15 <- predict(rfts)
y15_nofires$residm11d15 <- y15_nofires$pred_m11d15 - y15_nofires$heb_2015
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2015 # data year
cvtabof[1,2] <- 2011 # model year
cvtabof[1,3] <- mean(y15_nofires$residm11d15[y15_nofires$timesince > 199]) 
cvtabof[1,4] <- mean(y15_nofires$residm11d15[y15_nofires$timesince < 200])
d <- lm(y15_nofires$residm11d15 ~ y15_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y15_nofires$residm11d15[y15_nofires$timesince < 200] ~ y15_nofires$timesince[y15_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m11d15_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate 2011 model with 2017 data
rfts <- randomForest(heb_2017 ~ NDVstd1703 + ARGmea1706 + B3std1706, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y17_nofires$pred_m11d17 <- predict(rfts)
y17_nofires$residm11d17 <- y17_nofires$pred_m11d17 - y17_nofires$heb_2017
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2017 # data year
cvtabof[1,2] <- 2011 # model year
cvtabof[1,3] <- mean(y17_nofires$residm11d17[y17_nofires$timesince > 199]) 
cvtabof[1,4] <- mean(y17_nofires$residm11d17[y17_nofires$timesince < 200])
d <- lm(y17_nofires$residm11d17 ~ y17_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y17_nofires$residm11d17[y17_nofires$timesince < 200] ~ y17_nofires$timesince[y17_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m11d17_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate 2011 model with 2018 data
rfts <- randomForest(Herb_2018 ~ NDVstd1803 + ARGmea1806 + B3std1806, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y18_nofires$pred_m11d18 <- predict(rfts)
y18_nofires$residm11d18 <- y18_nofires$pred_m11d18 - y18_nofires$Herb_2018
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2018 # data year
cvtabof[1,2] <- 2011 # model year
cvtabof[1,3] <- mean(y18_nofires$residm11d18[y18_nofires$timesince > 199]) 
cvtabof[1,4] <- mean(y18_nofires$residm11d18[y18_nofires$timesince < 200])
d <- lm(y18_nofires$residm11d18 ~ y18_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y18_nofires$residm11d18[y18_nofires$timesince < 200] ~ y18_nofires$timesince[y18_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m11d18_mresid_oldfires.csv"
write.csv(cvtabof,tname)


# OLDFIRES ##### Generate residuals for model yr 2013 for all years & save tables in cvtables_oldfires -----------
# OLDFIRES - Generate 2013 model with 2009 data
rfts <- randomForest(heb_2009 ~ NDVstd0903 + ARGmea0906 + B3max0909, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y09_nofires$pred_m13d09 <- predict(rfts)
y09_nofires$residm13d09 <- y09_nofires$pred_m13d09 - y09_nofires$heb_2009
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2009 # data year
cvtabof[1,2] <- 2013 # model year
cvtabof[1,3] <- mean(y09_nofires$residm13d09[y09_nofires$timesince > 199])
cvtabof[1,4] <- mean(y09_nofires$residm13d09[y09_nofires$timesince < 200])
#plot(y09_nofires$residm11d09,y09_nofires$timesince)
d <- lm(y09_nofires$residm13d09 ~ y09_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y09_nofires$residm13d09[y09_nofires$timesince < 200] ~ y09_nofires$timesince[y09_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m13d09_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate 2013 model with 2011 data
rfts <- randomForest(heb_2011 ~ ARGmea1106 + NDVstd1103 + B3max1109, 
                        data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y11_nofires$pred_m13d11 <- predict(rfts)
y11_nofires$residm13d11 <- y11_nofires$pred_m13d11 - y11_nofires$heb_2011
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2011 # data year
cvtabof[1,2] <- 2013 # model year
cvtabof[1,3] <- mean(y11_nofires$residm13d11[y11_nofires$timesince > 199])
cvtabof[1,4] <- mean(y11_nofires$residm13d11[y11_nofires$timesince < 200])
#plot(y09_nofires$residm11d09,y09_nofires$timesince)
d <- lm(y11_nofires$residm13d11 ~ y11_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y11_nofires$residm13d11[y11_nofires$timesince < 200] ~ y11_nofires$timesince[y11_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m13d11_mresid_oldfires.csv"
write.csv(cvtabof,tname)


# OLDFIRES - Generate 2013 model with 2013 data: heb_2013 ~ ARGmea1306 + NDVstd1303 + B3max1309
rfts <- randomForest(heb_2013 ~ ARGmea1306 + NDVstd1303 + B3max1309, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y13_nofires$pred_m13d13 <- predict(rfts)
y13_nofires$residm13d13 <- y13_nofires$pred_m13d13 - y13_nofires$heb_2013
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2013 # data year
cvtabof[1,2] <- 2013 # model year
cvtabof[1,3] <- mean(y13_nofires$residm13d13[y13_nofires$timesince > 199])
cvtabof[1,4] <- mean(y13_nofires$residm13d13[y13_nofires$timesince < 200])
#plot(y09_nofires$residm11d09,y09_nofires$timesince)
d <- lm(y13_nofires$residm13d13 ~ y13_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y13_nofires$residm13d13[y13_nofires$timesince < 200] ~ y13_nofires$timesince[y13_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m13d13_mresid_oldfires.csv"
write.csv(cvtabof,tname)


# OLDFIRES - Generate 2013 model with 2015 data: heb_2015 ~ ARGmea1506 + NDVstd1503 + B3max1509
rfts <- randomForest(heb_2015 ~ ARGmea1506 + NDVstd1503 + B3max1509, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y15_nofires$pred_m13d15 <- predict(rfts)
y15_nofires$residm13d15 <- y15_nofires$pred_m13d15 - y15_nofires$heb_2015
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2015 # data year
cvtabof[1,2] <- 2013 # model year
cvtabof[1,3] <- mean(y15_nofires$residm13d15[y15_nofires$timesince > 199])
cvtabof[1,4] <- mean(y15_nofires$residm13d15[y15_nofires$timesince < 200])
#plot(y09_nofires$residm11d09,y09_nofires$timesince)
d <- lm(y15_nofires$residm13d15 ~ y15_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y15_nofires$residm13d15[y15_nofires$timesince < 200] ~ y15_nofires$timesince[y15_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m13d15_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate 2013 model with 2017 data: heb_2017 ~ ARGmea1706 + NDVstd1703 + B3max1709
rfts <- randomForest(heb_2017 ~ ARGmea1706 + NDVstd1703 + B3max1709, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y17_nofires$pred_m13d17 <- predict(rfts)
y17_nofires$residm13d17 <- y17_nofires$pred_m13d17 - y17_nofires$heb_2017
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2017 # data year
cvtabof[1,2] <- 2013 # model year
cvtabof[1,3] <- mean(y17_nofires$residm13d17[y17_nofires$timesince > 199])
cvtabof[1,4] <- mean(y17_nofires$residm13d17[y17_nofires$timesince < 200])
#plot(y09_nofires$residm11d09,y09_nofires$timesince)
d <- lm(y17_nofires$residm13d17 ~ y17_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y17_nofires$residm13d17[y17_nofires$timesince < 200] ~ y17_nofires$timesince[y17_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m13d17_mresid_oldfires.csv"
write.csv(cvtabof,tname)


# OLDFIRES - Generate 2013 model with 2018 data: Herb_2018 ~ ARGmea0906 + NDVstd0903 + B3max1809
rfts <- randomForest(Herb_2018 ~ NDVstd1803 + ARGmea1806 + B3max1809, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y18_nofires$pred_m13d18 <- predict(rfts)
y18_nofires$residm13d18 <- y18_nofires$pred_m13d18 - y18_nofires$Herb_2018
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2018 # data year
cvtabof[1,2] <- 2013 # model year
cvtabof[1,3] <- mean(y18_nofires$residm13d18[y18_nofires$timesince > 199])
cvtabof[1,4] <- mean(y18_nofires$residm13d18[y18_nofires$timesince < 200])
#plot(y09_nofires$residm11d09,y09_nofires$timesince)
d <- lm(y18_nofires$residm13d18 ~ y18_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y18_nofires$residm13d18[y18_nofires$timesince < 200] ~ y18_nofires$timesince[y18_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m13d18_mresid_oldfires.csv"
write.csv(cvtabof,tname)


# OLDFIRES ##### Generate residuals for model yr 2015 for all years & save tables in cvtables_oldfires ------------
# OLDFIRES - Generate 2015 model with 2009 data: heb_2009 ~  ARGmax0903 + ARGmax0906 + NDVmin0909
rfts <- randomForest(heb_2009 ~ ARGmax0903 + ARGmax0906 + NDVmin0909, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y09_nofires$pred_m15d09 <- predict(rfts)
y09_nofires$residm15d09 <- y09_nofires$pred_m15d09 - y09_nofires$heb_2009
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2009 # data year
cvtabof[1,2] <- 2015 # model year
cvtabof[1,3] <- mean(y09_nofires$residm15d09[y09_nofires$timesince > 199])
cvtabof[1,4] <- mean(y09_nofires$residm15d09[y09_nofires$timesince < 200])
d <- lm(y09_nofires$residm15d09 ~ y09_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y09_nofires$residm15d09[y09_nofires$timesince < 200] ~ y09_nofires$timesince[y09_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m15d09_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate 2015 model with 2011 data: heb_2011 ~ ARGmax1103 + ARGmax1106 + NDVmin1109
rfts <- randomForest(heb_2011 ~ ARGmax1103 + ARGmax1106 + NDVmin1109, 
                        data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y11_nofires$pred_m15d11 <- predict(rfts)
y11_nofires$residm15d11 <- y11_nofires$pred_m15d11 - y11_nofires$heb_2011
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2011 # data year
cvtabof[1,2] <- 2015 # model year
cvtabof[1,3] <- mean(y11_nofires$residm15d11[y11_nofires$timesince > 199])
cvtabof[1,4] <- mean(y11_nofires$residm15d11[y11_nofires$timesince < 200])
# Make linear model for residuals ~ timesince
d <- lm(y11_nofires$residm15d11 ~ y11_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y11_nofires$residm15d11[y11_nofires$timesince < 200] ~ y11_nofires$timesince[y11_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m15d11_mresid_oldfires.csv"
write.csv(cvtabof,tname)


# OLDFIRES - Generate 2015 model with 2013 data: heb_2013 ~ ARGmax1303 + ARGmax1306 + NDVmin1309
rfts <- randomForest(heb_2013 ~ ARGmax1303 + ARGmax1306 + NDVmin1309, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y13_nofires$pred_m15d13 <- predict(rfts)
y13_nofires$residm15d13 <- y13_nofires$pred_m15d13 - y13_nofires$heb_2013
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2013 # data year
cvtabof[1,2] <- 2015 # model year
cvtabof[1,3] <- mean(y13_nofires$residm15d13[y13_nofires$timesince > 199])
cvtabof[1,4] <- mean(y13_nofires$residm15d13[y13_nofires$timesince < 200])
# Make linear model for residuals ~ timesince
d <- lm(y13_nofires$residm15d13 ~ y13_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y13_nofires$residm15d13[y13_nofires$timesince < 200] ~ y13_nofires$timesince[y13_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m15d13_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate 2015 model with 2015 data: heb_2015 ~ ARGmax1503 + ARGmax1506 + NDVmin1509
rfts <- randomForest(heb_2015 ~ ARGmax1503 + ARGmax1506 + NDVmin1509, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y15_nofires$pred_m15d15 <- predict(rfts)
y15_nofires$residm15d15 <- y15_nofires$pred_m15d15 - y15_nofires$heb_2015
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2015 # data year
cvtabof[1,2] <- 2015 # model year
cvtabof[1,3] <- mean(y15_nofires$residm15d15[y15_nofires$timesince > 199])
cvtabof[1,4] <- mean(y15_nofires$residm15d15[y15_nofires$timesince < 200])
# Make linear model for residuals ~ timesince
d <- lm(y15_nofires$residm15d15 ~ y15_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y15_nofires$residm15d15[y15_nofires$timesince < 200] ~ y15_nofires$timesince[y15_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m15d15_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate 2015 model with 2017 data: heb_2017 ~  ARGmax1703 + ARGmax1706 + NDVmin1709
rfts <- randomForest(heb_2017 ~ ARGmax1703 + ARGmax1706 + NDVmin1709, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y17_nofires$pred_m15d17 <- predict(rfts)
y17_nofires$residm15d17 <- y17_nofires$pred_m15d17 - y17_nofires$heb_2017
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2017 # data year
cvtabof[1,2] <- 2015 # model year
cvtabof[1,3] <- mean(y17_nofires$residm15d17[y17_nofires$timesince > 199])
cvtabof[1,4] <- mean(y17_nofires$residm15d17[y17_nofires$timesince < 200])
# Make linear model for residuals ~ timesince
d <- lm(y17_nofires$residm15d17 ~ y17_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y17_nofires$residm15d17[y17_nofires$timesince < 200] ~ y17_nofires$timesince[y17_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m15d17_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate 2015 model with 2018 data: Herb_2018 ~ ARGmax1803 + ARGmax1806 + NDVmin1809
rfts <- randomForest(Herb_2018 ~ ARGmax1803 + ARGmax1806 + NDVmin1809, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y18_nofires$pred_m15d18 <- predict(rfts)
y18_nofires$residm15d18 <- y18_nofires$pred_m15d18 - y18_nofires$Herb_2018
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2018 # data year
cvtabof[1,2] <- 2015 # model year
cvtabof[1,3] <- mean(y18_nofires$residm15d18[y18_nofires$timesince > 199])
cvtabof[1,4] <- mean(y18_nofires$residm15d18[y18_nofires$timesince < 200])
# Make linear model for residuals ~ timesince
d <- lm(y18_nofires$residm15d18 ~ y18_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y18_nofires$residm15d18[y18_nofires$timesince < 200] ~ y18_nofires$timesince[y18_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m15d18_mresid_oldfires.csv"
write.csv(cvtabof,tname)


# OLDFIRES ##### Generate residuals for model yr 2017 for all years & save tables in cvtables_oldfires ----------
# OLDFIRES - Generate 2017 model with 2009 data: heb_2009 ~   ARGmax0906 + NBRstd0903 + ARGmed0909
rfts <- randomForest(heb_2009 ~ ARGmax0906 + NBRstd0903 + ARGmed0909, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y09_nofires$pred_m17d09 <- predict(rfts)
y09_nofires$residm17d09 <- y09_nofires$pred_m17d09 - y09_nofires$heb_2009
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2009 # data year
cvtabof[1,2] <- 2017 # model year
cvtabof[1,3] <- mean(y09_nofires$residm17d09[y09_nofires$timesince > 199])
cvtabof[1,4] <- mean(y09_nofires$residm17d09[y09_nofires$timesince < 200])
d <- lm(y09_nofires$residm17d09 ~ y09_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y09_nofires$residm17d09[y09_nofires$timesince < 200] ~ y09_nofires$timesince[y09_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m17d09_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate 2017 model with 2011 data: heb_2011 ~ ARGmax1106 + NBRstd1103 + ARGmed1109
rfts <- randomForest(heb_2011 ~ ARGmax1106 + NBRstd1103 + ARGmed1109, 
                        data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y11_nofires$pred_m17d11 <- predict(rfts)
y11_nofires$residm17d11 <- y11_nofires$pred_m17d11 - y11_nofires$heb_2011
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2011 # data year
cvtabof[1,2] <- 2017 # model year
cvtabof[1,3] <- mean(y11_nofires$residm17d11[y11_nofires$timesince > 199])
cvtabof[1,4] <- mean(y11_nofires$residm17d11[y11_nofires$timesince < 200])
d <- lm(y11_nofires$residm17d11 ~ y11_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y11_nofires$residm17d11[y11_nofires$timesince < 200] ~ y11_nofires$timesince[y11_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m17d11_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate 2017 model with 2013 data: heb_2013 ~ ARGmax1306 + NBRstd1303 + ARGmed1309
rfts <- randomForest(heb_2013 ~ ARGmax1306 + NBRstd1303 + ARGmed1309, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y13_nofires$pred_m17d13 <- predict(rfts)
y13_nofires$residm17d13 <- y13_nofires$pred_m17d13 - y13_nofires$heb_2013
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2013 # data year
cvtabof[1,2] <- 2017 # model year
cvtabof[1,3] <- mean(y13_nofires$residm17d13[y13_nofires$timesince > 199])
cvtabof[1,4] <- mean(y13_nofires$residm17d13[y13_nofires$timesince < 200])
d <- lm(y13_nofires$residm17d13 ~ y13_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y13_nofires$residm17d13[y13_nofires$timesince < 200] ~ y13_nofires$timesince[y13_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m17d13_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate 2017 model with 2015 data: heb_2015 ~ ARGmax1506 + NBRstd1503 + ARGmed1509
rfts <- randomForest(heb_2015 ~ ARGmax1506 + NBRstd1503 + ARGmed1509, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y15_nofires$pred_m17d15 <- predict(rfts)
y15_nofires$residm17d15 <- y15_nofires$pred_m17d15 - y15_nofires$heb_2015
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2015 # data year
cvtabof[1,2] <- 2017 # model year
cvtabof[1,3] <- mean(y15_nofires$residm17d15[y15_nofires$timesince > 199])
cvtabof[1,4] <- mean(y15_nofires$residm17d15[y15_nofires$timesince < 200])
d <- lm(y15_nofires$residm17d15 ~ y15_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y15_nofires$residm17d15[y15_nofires$timesince < 200] ~ y15_nofires$timesince[y15_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m17d15_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate 2017 model with 2017 data: heb_2017 ~  ARGmax1706 + NBRstd1703 + ARGmed1709
rfts <- randomForest(heb_2017 ~ ARGmax1706 + NBRstd1703 + ARGmed1709, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y17_nofires$pred_m17d17 <- predict(rfts)
y17_nofires$residm17d17 <- y17_nofires$pred_m17d17 - y17_nofires$heb_2017
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2017 # data year
cvtabof[1,2] <- 2017 # model year
cvtabof[1,3] <- mean(y17_nofires$residm17d17[y17_nofires$timesince > 199])
cvtabof[1,4] <- mean(y17_nofires$residm17d17[y17_nofires$timesince < 200])
d <- lm(y17_nofires$residm17d17 ~ y17_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y17_nofires$residm17d17[y17_nofires$timesince < 200] ~ y17_nofires$timesince[y17_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m17d17_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate 2017 model with 2018 data: Herb_2018 ~ ARGmax1806 + NBRstd1803 + ARGmed1809
rfts <- randomForest(Herb_2018 ~ ARGmax1806 + NBRstd1803 + ARGmed1809, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y18_nofires$pred_m17d18 <- predict(rfts)
y18_nofires$residm17d18 <- y18_nofires$pred_m17d18 - y18_nofires$Herb_2018
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2018 # data year
cvtabof[1,2] <- 2017 # model year
cvtabof[1,3] <- mean(y18_nofires$residm17d18[y18_nofires$timesince > 199])
cvtabof[1,4] <- mean(y18_nofires$residm17d18[y18_nofires$timesince < 200])
d <- lm(y18_nofires$residm17d18 ~ y18_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y18_nofires$residm17d18[y18_nofires$timesince < 200] ~ y18_nofires$timesince[y18_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m17d18_mresid_oldfires.csv"
write.csv(cvtabof,tname)


# OLDFIRES ##### Generate residuals for model yr 2018 for all years & save tables in cvtables_oldfires ----------
# OLDFIRES - Generate 2018 model with 2009 data: heb_2009 ~  ARGmax0906 + NBRstd1803 +  B3std1803
rfts <- randomForest(heb_2009 ~ ARGmax0906 + NBRstd0903 +  B3std0903, 
                        data = y09_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y09_nofires$pred_m18d09 <- predict(rfts)
y09_nofires$residm18d09 <- y09_nofires$pred_m18d09 - y09_nofires$heb_2009
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2009 # data year
cvtabof[1,2] <- 2018 # model year
cvtabof[1,3] <- mean(y09_nofires$residm18d09[y09_nofires$timesince > 199])
cvtabof[1,4] <- mean(y09_nofires$residm18d09[y09_nofires$timesince < 200])
d <- lm(y09_nofires$residm18d09 ~ y09_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y09_nofires$residm18d09[y09_nofires$timesince < 200] ~ y09_nofires$timesince[y09_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m18d09_mresid_oldfires.csv"
write.csv(cvtabof,tname)

#OLDFIRES - Generate 2018 model with 2011 data: heb_2011 ~ ARGmax1106 + NBRstd1103 + B3std1103
rfts <- randomForest(heb_2011 ~ ARGmax1106 + NBRstd1103 + B3std1103, 
                        data = y11_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y11_nofires$pred_m18d11 <- predict(rfts)
y11_nofires$residm18d11 <- y11_nofires$pred_m18d11 - y11_nofires$heb_2011
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2011 # data year
cvtabof[1,2] <- 2018 # model year
cvtabof[1,3] <- mean(y11_nofires$residm18d11[y11_nofires$timesince > 199])
cvtabof[1,4] <- mean(y11_nofires$residm18d11[y11_nofires$timesince < 200])
d <- lm(y11_nofires$residm18d11 ~ y11_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y11_nofires$residm18d11[y11_nofires$timesince < 200] ~ y11_nofires$timesince[y11_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m18d11_mresid_oldfires.csv"
write.csv(cvtabof,tname)


# OLDFIRES - Generate 2018 model with 2013 data: heb_2013 ~ ARGmax1306 + NBRstd1303 + B3std1303
rfts <- randomForest(heb_2013 ~ ARGmax1306 + NBRstd1303 + B3std1303, 
                        data = y13_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y13_nofires$pred_m18d13 <- predict(rfts)
y13_nofires$residm18d13 <- y13_nofires$pred_m18d13 - y13_nofires$heb_2013
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2013 # data year
cvtabof[1,2] <- 2018 # model year
cvtabof[1,3] <- mean(y13_nofires$residm18d13[y13_nofires$timesince > 199])
cvtabof[1,4] <- mean(y13_nofires$residm18d13[y13_nofires$timesince < 200])
d <- lm(y13_nofires$residm18d13 ~ y13_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y13_nofires$residm18d13[y13_nofires$timesince < 200] ~ y13_nofires$timesince[y13_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m18d13_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate 2018 model with 2015 data: heb_2015 ~ ARGmax1506 + NBRstd1503 + B3std1503
rfts <- randomForest(heb_2015 ~ ARGmax1506 + NBRstd1503 + ARGmed1509, 
                        data = y15_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y15_nofires$pred_m18d15 <- predict(rfts)
y15_nofires$residm18d15 <- y15_nofires$pred_m18d15 - y15_nofires$heb_2015
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2015 # data year
cvtabof[1,2] <- 2018 # model year
cvtabof[1,3] <- mean(y15_nofires$residm18d15[y15_nofires$timesince > 199])
cvtabof[1,4] <- mean(y15_nofires$residm18d15[y15_nofires$timesince < 200])
d <- lm(y15_nofires$residm18d15 ~ y15_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y15_nofires$residm18d15[y15_nofires$timesince < 200] ~ y15_nofires$timesince[y15_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m18d15_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate 2018 model with 2017 data: heb_2017 ~  ARGmax1706 + NBRstd1703 + B3std1703
rfts <- randomForest(heb_2017 ~ ARGmax1706 + NBRstd1703 + B3std1703, 
                        data = y17_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y17_nofires$pred_m18d17 <- predict(rfts)
y17_nofires$residm18d17 <- y17_nofires$pred_m18d17 - y17_nofires$heb_2017
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2017 # data year
cvtabof[1,2] <- 2018 # model year
cvtabof[1,3] <- mean(y17_nofires$residm18d17[y17_nofires$timesince > 199])
cvtabof[1,4] <- mean(y17_nofires$residm18d17[y17_nofires$timesince < 200])
d <- lm(y17_nofires$residm18d17 ~ y17_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y17_nofires$residm18d17[y17_nofires$timesince < 200] ~ y17_nofires$timesince[y17_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m18d17_mresid_oldfires.csv"
write.csv(cvtabof,tname)

# OLDFIRES - Generate 2018 model with 2018 data: Herb_2018 ~ ARGmax1806 + NBRstd1803 + B3std1803
rfts <- randomForest(Herb_2018 ~ ARGmax1806 + NBRstd1803 + ARGmed1809, 
                        data = y18_nofires, importance = TRUE, type = "regression", na.action = na.exclude)
y18_nofires$pred_m18d18 <- predict(rfts)
y18_nofires$residm18d18 <- y18_nofires$pred_m18d18 - y18_nofires$Herb_2018
# Separate by presence or absence of old fire. 'timesince' for NO old fires = 200
# Add means to table & write to folder N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires
cvtabof <- data.frame()
cvtabof[1,1] <- 2018 # data year
cvtabof[1,2] <- 2018 # model year
cvtabof[1,3] <- mean(y18_nofires$residm18d18[y18_nofires$timesince > 199])
cvtabof[1,4] <- mean(y18_nofires$residm18d18[y18_nofires$timesince < 200])
d <- lm(y18_nofires$residm18d18 ~ y18_nofires$timesince)
# add timesince Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y18_nofires$residm18d18[y18_nofires$timesince < 200] ~ y18_nofires$timesince[y18_nofires$timesince < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES old fires
cvtabof[1,8] <- glance(dYof)$p.value # YES old fires
colnames(cvtabof) <- c("datayr","modelyr","mresid_Nof","mresid_Yof","Aoflm_estimate","Aofpvalue","Yoflm_estimate","Yofpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_oldfires\\rfcv_m18d18_mresid_oldfires.csv"
write.csv(cvtabof,tname)






# OLDFIRES ##### Append residual model results into a single file ------------------------------
library(purrr)
library(readr)
library(mgsub)
# substrRight <- function(x, n){
#   substr(x, nchar(x)-n+1, nchar(x))
# }
setwd("N:/project/monitoring_vol/SoCalShrubs/EE_Extractions/cvtables_oldfires")
list_of_files <- list.files(path = "N:/project/monitoring_vol/SoCalShrubs/EE_Extractions/cvtables_oldfires",
                            pattern = "rfcv",
                            full.names = FALSE)


all_models_tempint_of <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read_csv(.x, col_types = cols(), col_names = TRUE), .id = "file_name") # apply read_csv to all elements 


#colnames(all_models_tempint) <- c("file_name","X1", "datayr","modelyr","mresid_Nof","mresid_Yof","lmestimate","pvalue")
all_models_tempint_of$datayr_ch <-as.character(all_models_tempint_of$datayr)
all_models_tempint_of$modelyr_ch <-as.character(all_models_tempint_of$modelyr)

write.csv(all_models_tempint_of, file = "N:/project/monitoring_vol/SoCalShrubs/EE_Extractions/cvtables_oldfires/oldfires_sum_20201015.csv")


# OLDFIRES ##### Heatmaps for for yes & no old fires -------------
# Mean residuals for NO old fires
round_mresid_nof <-round(reshape2::acast(all_models_tempint_of,modelyr~datayr,fun=mean,value.var="mresid_Nof"),2)
means <- as.matrix(rowMeans(round_mresid_nof))
out <- data.frame(round_mresid_nof,means)
pdf("NO_of_mean_resid_plot.pdf")
mresid_no_old_fires <- ggplot(data = all_models_tempint_of, mapping = aes(y = all_models_tempint_of$modelyr_ch, x = all_models_tempint_of$datayr_ch,fill = mresid_Nof)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(mresid_Nof,2))) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title="NO old fires", subtitle=stringr::str_wrap("Cell values = mean of residuals",width = 65)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(round_mresid_nof)[6],2), color="black",size=4) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(round_mresid_nof)[5],2), color="black",size=4) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(round_mresid_nof)[4],2), color="black",size=4) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(round_mresid_nof)[3],2), color="black",size=4) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(round_mresid_nof)[2],2), color="black",size=4) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(round_mresid_nof)[1],2), color="black",size=4) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(mresid_no_old_fires)
dev.off()

# Mean residuals for YES old fires
round_mresid_yof <-round(reshape2::acast(all_models_tempint_of,modelyr~datayr,fun=mean,value.var="mresid_Yof"),2)
means <- as.matrix(rowMeans(round_mresid_yof))
out <- data.frame(round_mresid_yof,means)
pdf("ONLY_of_mean_resid__plot.pdf")
mresid_yes_old_fires <- ggplot(data = all_models_tempint_of, mapping = aes(y = all_models_tempint_of$modelyr_ch, x = all_models_tempint_of$datayr_ch,fill = mresid_Yof)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(mresid_Yof,2))) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title="ONLY old fires", subtitle=stringr::str_wrap("Cell values = mean of residuals",width = 65)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(round_mresid_yof)[6],2), color="black",size=4) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(round_mresid_yof)[5],2), color="black",size=4) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(round_mresid_yof)[4],2), color="black",size=4) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(round_mresid_yof)[3],2), color="black",size=4) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(round_mresid_yof)[2],2), color="black",size=4) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(round_mresid_yof)[1],2), color="black",size=4) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(mresid_yes_old_fires)
dev.off()

# All (old fires / no old fires) lm esimate
rounded <-round(reshape2::acast(all_models_tempint_of,modelyr~datayr,fun=mean,value.var="Aoflm_estimate"),3)
means <- as.matrix(rowMeans(rounded))
#out <- data.frame(round_Aoflm_estimate,means)
pdf("All_of_lm_estimate_plot.pdf")
Aoflm_estimate <- ggplot(data = all_models_tempint_of, mapping = aes(y = all_models_tempint_of$modelyr_ch, x = all_models_tempint_of$datayr_ch,fill = Aoflm_estimate)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(Aoflm_estimate,3)),size=3.5) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title="All - both old fires & no old fires", subtitle=stringr::str_wrap("Cell values = lm estimate",width = 65)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(rounded)[6],3), color="black",size=3.3) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(rounded)[5],3), color="black",size=3.3) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(rounded)[4],3), color="black",size=3.3) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(rounded)[3],3), color="black",size=3.3) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(rounded)[2],3), color="black",size=3.3) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(rounded)[1],3), color="black",size=3.3) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(Aoflm_estimate)
dev.off()

# All (old fires / no old fires) p-value
rounded <-round(reshape2::acast(all_models_tempint_of,modelyr~datayr,fun=mean,value.var="Aofpvalue"),3)
#means <- as.matrix(rowMeans(rounded))
#out <- data.frame(round_Aoflm_estimate,means)
pdf("All_of_p_value_plot.pdf")
Aof_p_value <- ggplot(data = all_models_tempint_of, mapping = aes(y = all_models_tempint_of$modelyr_ch, x = all_models_tempint_of$datayr_ch,fill = Aofpvalue)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(Aofpvalue,3)),size=3.5) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title="All - both old fires & no old fires", subtitle=stringr::str_wrap("Cell values = p-value",width = 65)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(rounded)[6],2), color="black",size=3.3) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(rounded)[5],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(rounded)[4],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(rounded)[3],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(rounded)[2],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(rounded)[1],2), color="black",size=3.3) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(Aof_p_value)
dev.off()

# ONLY old fires lm estimate
rounded <-round(reshape2::acast(all_models_tempint_of,modelyr~datayr,fun=mean,value.var="Yoflm_estimate"),3)
#means <- as.matrix(rowMeans(rounded))
#out <- data.frame(round_Aoflm_estimate,means)
pdf("ONLY_of_lm_estimate_plot.pdf")
Yoflm_estimate <- ggplot(data = all_models_tempint_of, mapping = aes(y = all_models_tempint_of$modelyr_ch, x = all_models_tempint_of$datayr_ch,fill = Yoflm_estimate)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(Yoflm_estimate,3)),size=3.5) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title= "ONLY old fires) ", subtitle=stringr::str_wrap("Cell values = lm estimate",width = 65)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(rounded)[6],2), color="black",size=3.3) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(rounded)[5],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(rounded)[4],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(rounded)[3],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(rounded)[2],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(rounded)[1],2), color="black",size=3.3) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(Yoflm_estimate)
dev.off()

# ONLY old fires Yofpvalue
rounded <-round(reshape2::acast(all_models_tempint_of,modelyr~datayr,fun=mean,value.var="Yofpvalue"),3)
#means <- as.matrix(rowMeans(rounded))
#out <- data.frame(round_Aoflm_estimate,means)
pdf("ONLY_of_p_value_plot.pdf")
Yofpvalue <- ggplot(data = all_models_tempint_of, mapping = aes(y = all_models_tempint_of$modelyr_ch, x = all_models_tempint_of$datayr_ch,fill = Yofpvalue)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(Yofpvalue,3)),size=3.5) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title= "ONLY old fires) ", subtitle=stringr::str_wrap("Cell values = p value",width = 65)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(rounded)[6],2), color="black",size=3.3) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(rounded)[5],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(rounded)[4],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(rounded)[3],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(rounded)[2],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(rounded)[1],2), color="black",size=3.3) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(Yofpvalue)
dev.off()

# ONLY old fires Yoflm_estimate
rounded <-round(reshape2::acast(all_models_tempint_of,modelyr~datayr,fun=mean,value.var="Yoflm_estimate"),3)
#means <- as.matrix(rowMeans(rounded))
#out <- data.frame(round_Aoflm_estimate,means)
pdf("ONLY_of_lm_estimate_plot.pdf")
Yoflm_estimate <- ggplot(data = all_models_tempint_of, mapping = aes(y = all_models_tempint_of$modelyr_ch, x = all_models_tempint_of$datayr_ch,fill = Yoflm_estimate)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(Yoflm_estimate,3)),size=3.5) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title= "ONLY old fires) ", subtitle=stringr::str_wrap("Cell values = lm_estimate",width = 65)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(rounded)[6],2), color="black",size=3.3) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(rounded)[5],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(rounded)[4],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(rounded)[3],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(rounded)[2],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(rounded)[1],2), color="black",size=3.3) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(Yoflm_estimate)
dev.off()





# ------------------------------------------------------------------------------------------------------------
# OLDFIRES INCLUDE ALL ##### Generate rf residuals for model yr 2009 for all years & save tables as rfcv_**_mresid_all.csv#####-------
# These input data include all records, regardless of whether there was a fire or not
# 20201117 Rerun with the time since variables in each model & use output folder cvtables_20201104
# 20201118 Add BP test & bp test p value & generation of rmse of residuals. "rmseresid_A","BPresult","BPpvalue"
require(broom)
# OLDFIRES INCLUDE ALL - Generate residuals for model yr 2009 with 2009 data
rfts <- randomForest(heb_2009 ~ ARGmax0903 + NBRstd0903 + ARGmed0906 + ts2009, 
                     data = y09, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y09$pred_m09d09 <- predict(rfts)
y09$residm09d09 <- y09$pred_m09d09 - y09$heb_2009

sqrt(mean((y09$heb_2009 - y09$pred_m09d09)^2))
# Separate by presence or absence of past fires. 'timesince' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2009 # data year
cvtabof[1,2] <- 2009 # model year
cvtabof[1,3] <- mean(y09$residm09d09[y09$ts2009 > 199])
cvtabof[1,4] <- mean(y09$residm09d09[y09$ts2009 < 200])
d <- lm(y09$residm09d09 ~ y09$ts2009) # lm for ALL recs
# add ts2009 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2] #lm estimate
cvtabof[1,6] <- glance(d)$p.value # p.value
dYof <- lm(y09$residm09d09[y09$ts2009 < 200] ~ y09$ts2009[y09$ts2009 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y09$residm09d09) # mresid_All
cvtabof[1,10] <- rmse(y09$heb_2009,y09$pred_m09d09) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_20201104\\rfcv_m09d09_mresid_all.csv"
write.csv(cvtabof,tname)


# OLDFIRES INCLUDE ALL - Generate residuals for model year 2009 with 2011 data
rfts <- randomForest(heb_2011 ~ ARGmax1103 + NBRstd1103 + ARGmed1106 + ts2011, 
                     data = y11, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
# Generate residuals
y11$pred_m09d11 <- predict(rfts)
y11$residm09d11 <- y11$pred_m09d11 - y11$heb_2011
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2011 # data year
cvtabof[1,2] <- 2009 # model year
cvtabof[1,3] <- mean(y11$residm09d11[y11$ts2011 > 199]) # no past fires
cvtabof[1,4] <- mean(y11$residm09d11[y11$ts2011 < 200]) # yes past fires
d <- lm(y11$residm09d11 ~ y11$ts2011)
# add ts2011 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y11$residm09d11[y11$ts2011 < 200] ~ y11$ts2011[y11$ts2011 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y11$residm09d11)
cvtabof[1,10] <- rmse(y11$heb_2011,y11$pred_m09d11) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_20201104\\rfcv_m09d11_mresid_all.csv"
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate residuals for model year 2009 with 2013 data
rfts <- randomForest(heb_2013 ~ ARGmax1303 + NBRstd1303 + ARGmed1306 + ts2013, 
                     data = y13, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y13$pred_m09d13 <- predict(rfts)
y13$residm09d13 <- y13$pred_m09d13 - y13$heb_2013
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2013 # data year
cvtabof[1,2] <- 2009 # model year
cvtabof[1,3] <- mean(y13$residm09d13[y13$ts2013 > 199])
cvtabof[1,4] <- mean(y13$residm09d13[y13$ts2013 < 200])
d <- lm(y13$residm09d13 ~ y13$ts2013)
# add ts2013 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y13$residm09d13[y13$ts2013 < 200] ~ y13$ts2013[y13$ts2013 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y13$residm09d13)
cvtabof[1,10] <- rmse(y13$heb_2013,y13$pred_m09d13) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_20201104\\rfcv_m09d13_mresid_all.csv"
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate residuals for model year 2009 with 2015 data
rfts <- randomForest(heb_2015 ~ ARGmax1503 + NBRstd1503 + ARGmed1506 + ts2015, 
                     data = y15, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y15$pred_m09d15 <- predict(rfts)
y15$residm09d15 <- y15$pred_m09d15 - y15$heb_2015
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2015 # data year
cvtabof[1,2] <- 2009 # model year
cvtabof[1,3] <- mean(y15$residm09d15[y15$ts2015 > 199])
cvtabof[1,4] <- mean(y15$residm09d15[y15$ts2015 < 200])
d <- lm(y15$residm09d15 ~ y15$ts2015)
# add ts2015 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y15$residm09d15[y15$ts2015 < 200] ~ y15$ts2015[y15$ts2015 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y15$residm09d15)
cvtabof[1,10] <- rmse(y15$heb_2015,y15$pred_m09d15) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_20201104\\rfcv_m09d15_mresid_all.csv"
write.csv(cvtabof,tname)

#  OLDFIRES INCLUDE ALL - Generate residuals for model year 2009 with 2017 data
rfts <- randomForest(heb_2017 ~ ARGmax1703 + NBRstd1703 + ARGmed1706 + ts2017, 
                     data = y17, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y17$pred_m09d17 <- predict(rfts)
y17$residm09d17 <- y17$pred_m09d17 - y17$heb_2017
# Separate by presence or absence of past fires. 'ts2017' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2017 # data year
cvtabof[1,2] <- 2009 # model year
cvtabof[1,3] <- mean(y17$residm09d17[y17$ts2017 > 199]) # Select where we have no fire history
cvtabof[1,4] <- mean(y17$residm09d17[y17$ts2017 < 200]) # select where there are previous fires
d <- lm(y17$residm09d17 ~ y17$ts2017)
# add ts2017 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y17$residm09d17[y17$ts2017 < 200] ~ y17$ts2017[y17$ts2017 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y17$residm09d17)
cvtabof[1,10] <- rmse(y17$pred_m09d17,y17$pred_m09d17) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_20201104\\rfcv_m09d17_mresid_all.csv"
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate residuals for model year 2009 with 2018 data
rfts <- randomForest(Herb_2018 ~ ARGmax1803 + NBRstd1803 + ARGmed1806 + ts2018, 
                     data = y18, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y18$pred_m09d18 <- predict(rfts)
y18$residm09d18 <- y18$pred_m09d18 - y18$Herb_2018
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2018 # data year
cvtabof[1,2] <- 2009 # model year
cvtabof[1,3] <- mean(y18$residm09d18[y18$ts2018 > 199])
cvtabof[1,4] <- mean(y18$residm09d18[y18$ts2018 < 200])
d <- lm(y18$residm09d18 ~ y18$ts2018)
# add ts2018 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y18$residm09d18[y18$ts2018 < 200] ~ y18$ts2018[y18$ts2018 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y18$residm09d18)
cvtabof[1,10] <- rmse(y18$Herb_2018,y18$pred_m09d18) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_20201104\\rfcv_m09d18_mresid_all.csv"
write.csv(cvtabof,tname)


# OLDFIRES INCLUDE ALL ##### Generate rf residuals for model yr 2011 for all years & save tables as rfcv_**_mresid_all.csv #####-------
# OLDFIRES  INCLUDE ALL- Generate 2011 model with 2009 data
rfts <- randomForest(heb_2009 ~ NDVstd0903 + ARGmea0906 + B3std0906 + ts2009, 
                     data = y09, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y09$pred_m11d09 <- predict(rfts)
y09$residm11d09 <- y09$pred_m11d09 - y09$heb_2009
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2009 # data year
cvtabof[1,2] <- 2011 # model year
cvtabof[1,3] <- mean(y09$residm11d09[y09$ts2009 > 199])
cvtabof[1,4] <- mean(y09$residm11d09[y09$ts2009 < 200])
#plot(y09$residm11d09,y09$ts2009)
d <- lm(y09$residm11d09 ~ y09$ts2009)
# add ts2009 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y09$residm11d09[y09$ts2009 < 200] ~ y09$ts2009[y09$ts2009 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y09$residm11d09)
cvtabof[1,10] <- rmse(y09$heb_2009,y09$pred_m11d09) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_20201104\\rfcv_m11d09_mresid_all.csv"
write.csv(cvtabof,tname)


# OLDFIRES INCLUDE ALL - Generate 2011 model with 2011 data
rfts <- randomForest(heb_2011 ~ NDVstd1103 + ARGmea1106 + B3std1106 + ts2011, 
                     data = y11, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y11$pred_m11d11 <- predict(rfts)
y11$residm11d11 <- y11$pred_m11d11 - y11$heb_2011
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2011 # data year
cvtabof[1,2] <- 2011 # model year
cvtabof[1,3] <- mean(y11$residm11d11[y11$ts2011 > 199])
cvtabof[1,4] <- mean(y11$residm11d11[y11$ts2011 < 200])
d <- lm(y11$residm11d11 ~ y11$ts2011)
# add ts2011 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y11$residm11d11[y11$ts2011 < 200] ~ y11$ts2011[y11$ts2011 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y11$residm11d11)
cvtabof[1,10] <- rmse(y11$heb_2011,y11$pred_m11d11) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_20201104\\rfcv_m11d11_mresid_all.csv"
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate 2011 model with 2013 data
rfts <- randomForest(heb_2013 ~ NDVstd1303 + ARGmea1306 + B3std1306 + ts2013, 
                     data = y13, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y13$pred_m11d13 <- predict(rfts)
y13$residm11d13 <- y13$pred_m11d13 - y13$heb_2013
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2013 # data year
cvtabof[1,2] <- 2011 # model year
cvtabof[1,3] <- mean(y13$residm11d13[y13$ts2013 > 199]) 
cvtabof[1,4] <- mean(y13$residm11d13[y13$ts2013 < 200])
d <- lm(y13$residm11d13 ~ y13$ts2013)
# add ts2013 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y13$residm11d13[y13$ts2013 < 200] ~ y13$ts2013[y13$ts2013 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y13$residm11d13)
cvtabof[1,10] <- rmse(y13$heb_2013,y13$pred_m11d13) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_20201104\\rfcv_m11d13_mresid_all.csv"
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate 2011 model with 2015 data
rfts <- randomForest(heb_2015 ~ NDVstd1503 + ARGmea1506 + B3std1506 + ts2015, 
                     data = y15, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y15$pred_m11d15 <- predict(rfts)
y15$residm11d15 <- y15$pred_m11d15 - y15$heb_2015
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2015 # data year
cvtabof[1,2] <- 2011 # model year
cvtabof[1,3] <- mean(y15$residm11d15[y15$ts2015 > 199]) 
cvtabof[1,4] <- mean(y15$residm11d15[y15$ts2015 < 200])
d <- lm(y15$residm11d15 ~ y15$ts2015)
# add ts2015 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y15$residm11d15[y15$ts2015 < 200] ~ y15$ts2015[y15$ts2015 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y15$residm11d15)
cvtabof[1,10] <- rmse(y15$heb_2015,y15$pred_m11d15) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_20201104\\rfcv_m11d15_mresid_all.csv"
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate 2011 model with 2017 data
rfts <- randomForest(heb_2017 ~ NDVstd1703 + ARGmea1706 + B3std1706 + ts2017, 
                     data = y17, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y17$pred_m11d17 <- predict(rfts)
y17$residm11d17 <- y17$pred_m11d17 - y17$heb_2017
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2017 # data year
cvtabof[1,2] <- 2011 # model year
cvtabof[1,3] <- mean(y17$residm11d17[y17$ts2017 > 199]) 
cvtabof[1,4] <- mean(y17$residm11d17[y17$ts2017 < 200])
d <- lm(y17$residm11d17 ~ y17$ts2017)
# add ts2017 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y17$residm11d17[y17$ts2017 < 200] ~ y17$ts2017[y17$ts2017 < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y17$residm11d17)
cvtabof[1,10] <- rmse(y17$heb_2017,y17$pred_m11d17) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_20201104\\rfcv_m11d17_mresid_all.csv"
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate 2011 model with 2018 data
rfts <- randomForest(Herb_2018 ~ NDVstd1803 + ARGmea1806 + B3std1806 + ts2018, 
                     data = y18, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y18$pred_m11d18 <- predict(rfts)
y18$residm11d18 <- y18$pred_m11d18 - y18$Herb_2018
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2018 # data year
cvtabof[1,2] <- 2011 # model year
cvtabof[1,3] <- mean(y18$residm11d18[y18$ts2018 > 199]) 
cvtabof[1,4] <- mean(y18$residm11d18[y18$ts2018 < 200])
d <- lm(y18$residm11d18 ~ y18$ts2018)
# add ts2018 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y18$residm11d18[y18$ts2018 < 200] ~ y18$ts2018[y18$ts2018 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y18$residm11d18)
cvtabof[1,10] <- rmse(y18$Herb_2018,y18$pred_m11d18) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_20201104\\rfcv_m11d18_mresid_all.csv"
write.csv(cvtabof,tname)


# OLDFIRES INCLUDE ALL ##### Generate rf residuals for model yr 2013 for all years & save tables as rfcv_**_mresid_all.csv -----------
# OLDFIRES INCLUDE ALL - Generate 2013 model with 2009 data
folder <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_20201104"
rfts <- randomForest(heb_2009 ~ NDVstd0903 + ARGmea0906 + B3max0909 + ts2009, 
                     data = y09, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y09$pred_m13d09 <- predict(rfts)
y09$residm13d09 <- y09$pred_m13d09 - y09$heb_2009
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2009 # data year
cvtabof[1,2] <- 2013 # model year
cvtabof[1,3] <- mean(y09$residm13d09[y09$ts2009 > 199])
cvtabof[1,4] <- mean(y09$residm13d09[y09$ts2009 < 200])
#plot(y09$residm11d09,y09$ts2009)
d <- lm(y09$residm13d09 ~ y09$ts2009)
# add ts2009 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y09$residm13d09[y09$ts2009 < 200] ~ y09$ts2009[y09$ts2009 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y09$residm13d09)
cvtabof[1,10] <- rmse(y09$heb_2009,y09$pred_m13d09) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m13d09_mresid_all.csv")
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate 2013 model with 2011 data
rfts <- randomForest(heb_2011 ~ ARGmea1106 + NDVstd1103 + B3max1109 + ts2011, 
                     data = y11, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y11$pred_m13d11 <- predict(rfts)
y11$residm13d11 <- y11$pred_m13d11 - y11$heb_2011
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2011 # data year
cvtabof[1,2] <- 2013 # model year
cvtabof[1,3] <- mean(y11$residm13d11[y11$ts2011 > 199])
cvtabof[1,4] <- mean(y11$residm13d11[y11$ts2011 < 200])
#plot(y09$residm11d09,y09$ts2011)
d <- lm(y11$residm13d11 ~ y11$ts2011)
# add ts2011 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y11$residm13d11[y11$ts2011 < 200] ~ y11$ts2011[y11$ts2011 < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y11$residm13d11)
cvtabof[1,10] <- rmse(y11$heb_2011,y11$pred_m13d11) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m13d11_mresid_all.csv")
write.csv(cvtabof,tname)


# OLDFIRES INCLUDE ALL - Generate 2013 model with 2013 data: heb_2013 ~ ARGmea1306 + NDVstd1303 + B3max1309
rfts <- randomForest(heb_2013 ~ ARGmea1306 + NDVstd1303 + B3max1309 + ts2013, 
                     data = y13, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y13$pred_m13d13 <- predict(rfts)
y13$residm13d13 <- y13$pred_m13d13 - y13$heb_2013
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2013 # data year
cvtabof[1,2] <- 2013 # model year
cvtabof[1,3] <- mean(y13$residm13d13[y13$ts2013 > 199])
cvtabof[1,4] <- mean(y13$residm13d13[y13$ts2013 < 200])
#plot(y09$residm11d09,y09$ts2013)
d <- lm(y13$residm13d13 ~ y13$ts2013)
# add ts2013 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y13$residm13d13[y13$ts2013 < 200] ~ y13$ts2013[y13$ts2013 < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y13$residm13d13)
cvtabof[1,10] <- rmse(y13$heb_2013,y13$pred_m13d13) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m13d13_mresid_all.csv")
write.csv(cvtabof,tname)


# OLDFIRES INCLUDE ALL - Generate 2013 model with 2015 data: heb_2015 ~ ARGmea1506 + NDVstd1503 + B3max1509
rfts <- randomForest(heb_2015 ~ ARGmea1506 + NDVstd1503 + B3max1509 + ts2015, 
                     data = y15, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y15$pred_m13d15 <- predict(rfts)
y15$residm13d15 <- y15$pred_m13d15 - y15$heb_2015
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2015 # data year
cvtabof[1,2] <- 2013 # model year
cvtabof[1,3] <- mean(y15$residm13d15[y15$ts2015 > 199])
cvtabof[1,4] <- mean(y15$residm13d15[y15$ts2015 < 200])
#plot(y09$residm11d09,y09$ts2015)
d <- lm(y15$residm13d15 ~ y15$ts2015)
# add ts2015 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y15$residm13d15[y15$ts2015 < 200] ~ y15$ts2015[y15$ts2015 < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y15$residm13d15)
cvtabof[1,10] <- rmse(y15$heb_2015,y15$pred_m13d15) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m13d15_mresid_all.csv")
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate 2013 model with 2017 data: heb_2017 ~ ARGmea1706 + NDVstd1703 + B3max1709
rfts <- randomForest(heb_2017 ~ ARGmea1706 + NDVstd1703 + B3max1709 + ts2017, 
                     data = y17, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y17$pred_m13d17 <- predict(rfts)
y17$residm13d17 <- y17$pred_m13d17 - y17$heb_2017
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2017 # data year
cvtabof[1,2] <- 2013 # model year
cvtabof[1,3] <- mean(y17$residm13d17[y17$ts2017 > 199])
cvtabof[1,4] <- mean(y17$residm13d17[y17$ts2017 < 200])
#plot(y09$residm11d09,y09$ts2017)
d <- lm(y17$residm13d17 ~ y17$ts2017)
# add ts2017 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y17$residm13d17[y17$ts2017 < 200] ~ y17$ts2017[y17$ts2017 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y17$residm13d17)
cvtabof[1,10] <- rmse(y17$heb_2017,y17$pred_m13d17) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m13d17_mresid_all.csv")
write.csv(cvtabof,tname)


# OLDFIRES INCLUDE ALL - Generate 2013 model with 2018 data: Herb_2018 ~ ARGmea0906 + NDVstd0903 + B3max1809
rfts <- randomForest(Herb_2018 ~ NDVstd1803 + ARGmea1806 + B3max1809 + ts2018, 
                     data = y18, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y18$pred_m13d18 <- predict(rfts)
y18$residm13d18 <- y18$pred_m13d18 - y18$Herb_2018
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2018 # data year
cvtabof[1,2] <- 2013 # model year
cvtabof[1,3] <- mean(y18$residm13d18[y18$ts2018 > 199])
cvtabof[1,4] <- mean(y18$residm13d18[y18$ts2018 < 200])
#plot(y09$residm11d09,y09$ts2018)
d <- lm(y18$residm13d18 ~ y18$ts2018)
# add ts2018 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y18$residm13d18[y18$ts2018 < 200] ~ y18$ts2018[y18$ts2018 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y18$residm13d18)
cvtabof[1,10] <- rmse(y18$Herb_2018,y18$pred_m13d18) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m13d18_mresid_all.csv")
write.csv(cvtabof,tname)


# OLDFIRES INCLUDE ALL ##### Generate rf residuals for model yr 2015 for all years & save tables as rfcv_**_mresid_all.csv ------------
# OLDFIRES INCLUDE ALL - Generate 2015 model with 2009 data: heb_2009 ~  ARGmax0903 + ARGmax0906 + NDVmin0909
folder <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_20201104"
rfts <- randomForest(heb_2009 ~ ARGmax0903 + ARGmax0906 + NDVmin0909 + ts2009, 
                     data = y09, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y09$pred_m15d09 <- predict(rfts)
y09$residm15d09 <- y09$pred_m15d09 - y09$heb_2009
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2009 # data year
cvtabof[1,2] <- 2015 # model year
cvtabof[1,3] <- mean(y09$residm15d09[y09$ts2009 > 199])
cvtabof[1,4] <- mean(y09$residm15d09[y09$ts2009 < 200])
d <- lm(y09$residm15d09 ~ y09$ts2009)
# add ts2009 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y09$residm15d09[y09$ts2009 < 200] ~ y09$ts2009[y09$ts2009 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y09$residm15d09)
cvtabof[1,10] <- rmse(y09$heb_2009,y09$pred_m15d09) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m15d09_mresid_all.csv")
write.csv(cvtabof,tname)


# OLDFIRES INCLUDE ALL - Generate 2015 model with 2011 data: heb_2011 ~ ARGmax1103 + ARGmax1106 + NDVmin1109
rfts <- randomForest(heb_2011 ~ ARGmax1103 + ARGmax1106 + NDVmin1109 + ts2011, 
                     data = y11, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y11$pred_m15d11 <- predict(rfts)
y11$residm15d11 <- y11$pred_m15d11 - y11$heb_2011
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2011 # data year
cvtabof[1,2] <- 2015 # model year
cvtabof[1,3] <- mean(y11$residm15d11[y11$ts2011 > 199])
cvtabof[1,4] <- mean(y11$residm15d11[y11$ts2011 < 200])
# Make linear model for residuals ~ ts2011
d <- lm(y11$residm15d11 ~ y11$ts2011)
# add ts2011 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y11$residm15d11[y11$ts2011 < 200] ~ y11$ts2011[y11$ts2011 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y11$residm15d11)
cvtabof[1,10] <- rmse(y11$heb_2011,y11$pred_m15d11) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m15d11_mresid_all.csv")
write.csv(cvtabof,tname)


# OLDFIRES INCLUDE ALL - Generate 2015 model with 2013 data: heb_2013 ~ ARGmax1303 + ARGmax1306 + NDVmin1309
rfts <- randomForest(heb_2013 ~ ARGmax1303 + ARGmax1306 + NDVmin1309 + ts2013, 
                     data = y13, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y13$pred_m15d13 <- predict(rfts)
y13$residm15d13 <- y13$pred_m15d13 - y13$heb_2013
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2013 # data year
cvtabof[1,2] <- 2015 # model year
cvtabof[1,3] <- mean(y13$residm15d13[y13$ts2013 > 199])
cvtabof[1,4] <- mean(y13$residm15d13[y13$ts2013 < 200])
# Make linear model for residuals ~ ts2013
d <- lm(y13$residm15d13 ~ y13$ts2013)
# add ts2013 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y13$residm15d13[y13$ts2013 < 200] ~ y13$ts2013[y13$ts2013 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y13$residm15d13)
cvtabof[1,10] <- rmse(y13$heb_2013,y13$pred_m15d13) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m15d13_mresid_all.csv")
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate 2015 model with 2015 data: heb_2015 ~ ARGmax1503 + ARGmax1506 + NDVmin1509
rfts <- randomForest(heb_2015 ~ ARGmax1503 + ARGmax1506 + NDVmin1509 + ts2015, 
                     data = y15, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y15$pred_m15d15 <- predict(rfts)
y15$residm15d15 <- y15$pred_m15d15 - y15$heb_2015
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2015 # data year
cvtabof[1,2] <- 2015 # model year
cvtabof[1,3] <- mean(y15$residm15d15[y15$ts2015 > 199])
cvtabof[1,4] <- mean(y15$residm15d15[y15$ts2015 < 200])
# Make linear model for residuals ~ ts2015
d <- lm(y15$residm15d15 ~ y15$ts2015)
# add ts2015 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y15$residm15d15[y15$ts2015 < 200] ~ y15$ts2015[y15$ts2015 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y15$residm15d15)
cvtabof[1,10] <- rmse(y15$heb_2015,y15$pred_m15d15) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m15d15_mresid_all.csv")
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate 2015 model with 2017 data: heb_2017 ~  ARGmax1703 + ARGmax1706 + NDVmin1709
rfts <- randomForest(heb_2017 ~ ARGmax1703 + ARGmax1706 + NDVmin1709 + ts2017, 
                     data = y17, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y17$pred_m15d17 <- predict(rfts)
y17$residm15d17 <- y17$pred_m15d17 - y17$heb_2017
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2017 # data year
cvtabof[1,2] <- 2015 # model year
cvtabof[1,3] <- mean(y17$residm15d17[y17$ts2017 > 199])
cvtabof[1,4] <- mean(y17$residm15d17[y17$ts2017 < 200])
# Make linear model for residuals ~ ts2017
d <- lm(y17$residm15d17 ~ y17$ts2017)
# add ts2017 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y17$residm15d17[y17$ts2017 < 200] ~ y17$ts2017[y17$ts2017 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y17$residm15d17)
cvtabof[1,10] <- rmse(y17$heb_2017,y17$pred_m15d17) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m15d17_mresid_all.csv")
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate 2015 model with 2018 data: Herb_2018 ~ ARGmax1803 + ARGmax1806 + NDVmin1809
rfts <- randomForest(Herb_2018 ~ ARGmax1803 + ARGmax1806 + NDVmin1809 + ts2018, 
                     data = y18, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y18$pred_m15d18 <- predict(rfts)
y18$residm15d18 <- y18$pred_m15d18 - y18$Herb_2018
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2018 # data year
cvtabof[1,2] <- 2015 # model year
cvtabof[1,3] <- mean(y18$residm15d18[y18$ts2018 > 199])
cvtabof[1,4] <- mean(y18$residm15d18[y18$ts2018 < 200])
# Make linear model for residuals ~ ts2018
d <- lm(y18$residm15d18 ~ y18$ts2018)
# add ts2018 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y18$residm15d18[y18$ts2018 < 200] ~ y18$ts2018[y18$ts2018 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y18$residm15d18)
cvtabof[1,10] <- rmse(y18$Herb_2018,y18$pred_m15d18) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m15d18_mresid_all.csv")
write.csv(cvtabof,tname)


# OLDFIRES INCLUDE ALL ##### Generate rf residuals for model yr 2017 for all years & save tables as rfcv_**_mresid_all.csv ----------
# OLDFIRES INCLUDE ALL - Generate 2017 model with 2009 data: heb_2009 ~   ARGmax0906 + NBRstd0903 + ARGmed0909
folder <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_20201104"
rfts <- randomForest(heb_2009 ~ ARGmax0906 + NBRstd0903 + ARGmed0909 + ts2009, 
                     data = y09, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y09$pred_m17d09 <- predict(rfts)
y09$residm17d09 <- y09$pred_m17d09 - y09$heb_2009
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2009 # data year
cvtabof[1,2] <- 2017 # model year
cvtabof[1,3] <- mean(y09$residm17d09[y09$ts2009 > 199])
cvtabof[1,4] <- mean(y09$residm17d09[y09$ts2009 < 200])
d <- lm(y09$residm17d09 ~ y09$ts2009)
# add ts2009 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y09$residm17d09[y09$ts2009 < 200] ~ y09$ts2009[y09$ts2009 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y09$residm17d09)
cvtabof[1,10] <- rmse(y09$heb_2009,y09$pred_m17d09) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m17d09_mresid_all.csv")
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate 2017 model with 2011 data: heb_2011 ~ ARGmax1106 + NBRstd1103 + ARGmed1109
rfts <- randomForest(heb_2011 ~ ARGmax1106 + NBRstd1103 + ARGmed1109 + ts2011, 
                     data = y11, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y11$pred_m17d11 <- predict(rfts)
y11$residm17d11 <- y11$pred_m17d11 - y11$heb_2011
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2011 # data year
cvtabof[1,2] <- 2017 # model year
cvtabof[1,3] <- mean(y11$residm17d11[y11$ts2011 > 199])
cvtabof[1,4] <- mean(y11$residm17d11[y11$ts2011 < 200])
d <- lm(y11$residm17d11 ~ y11$ts2011)
# add ts2011 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y11$residm17d11[y11$ts2011 < 200] ~ y11$ts2011[y11$ts2011 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y11$residm17d11)
cvtabof[1,10] <- rmse(y11$heb_2011,y11$pred_m17d11) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m17d11_mresid_all.csv")
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate 2017 model with 2013 data: heb_2013 ~ ARGmax1306 + NBRstd1303 + ARGmed1309
rfts <- randomForest(heb_2013 ~ ARGmax1306 + NBRstd1303 + ARGmed1309 + ts2013, 
                     data = y13, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y13$pred_m17d13 <- predict(rfts)
y13$residm17d13 <- y13$pred_m17d13 - y13$heb_2013
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2013 # data year
cvtabof[1,2] <- 2017 # model year
cvtabof[1,3] <- mean(y13$residm17d13[y13$ts2013 > 199])
cvtabof[1,4] <- mean(y13$residm17d13[y13$ts2013 < 200])
d <- lm(y13$residm17d13 ~ y13$ts2013)
# add ts2013 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y13$residm17d13[y13$ts2013 < 200] ~ y13$ts2013[y13$ts2013 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y13$residm17d13)
cvtabof[1,10] <- rmse(y13$heb_2013,y13$pred_m17d13) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m17d13_mresid_all.csv")
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate 2017 model with 2015 data: heb_2015 ~ ARGmax1506 + NBRstd1503 + ARGmed1509
rfts <- randomForest(heb_2015 ~ ARGmax1506 + NBRstd1503 + ARGmed1509 + ts2015, 
                     data = y15, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y15$pred_m17d15 <- predict(rfts)
y15$residm17d15 <- y15$pred_m17d15 - y15$heb_2015
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2015 # data year
cvtabof[1,2] <- 2017 # model year
cvtabof[1,3] <- mean(y15$residm17d15[y15$ts2015 > 199])
cvtabof[1,4] <- mean(y15$residm17d15[y15$ts2015 < 200])
d <- lm(y15$residm17d15 ~ y15$ts2015)
# add ts2015 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y15$residm17d15[y15$ts2015 < 200] ~ y15$ts2015[y15$ts2015 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y15$residm17d15)
cvtabof[1,10] <- rmse(y13$heb_2013,y13$pred_m17d13) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m17d15_mresid_all.csv")
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate 2017 model with 2017 data: heb_2017 ~  ARGmax1706 + NBRstd1703 + ARGmed1709
rfts <- randomForest(heb_2017 ~ ARGmax1706 + NBRstd1703 + ARGmed1709 + ts2017, 
                     data = y17, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y17$pred_m17d17 <- predict(rfts)
y17$residm17d17 <- y17$pred_m17d17 - y17$heb_2017
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2017 # data year
cvtabof[1,2] <- 2017 # model year
cvtabof[1,3] <- mean(y17$residm17d17[y17$ts2017 > 199])
cvtabof[1,4] <- mean(y17$residm17d17[y17$ts2017 < 200])
d <- lm(y17$residm17d17 ~ y17$ts2017)
# add ts2017 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y17$residm17d17[y17$ts2017 < 200] ~ y17$ts2017[y17$ts2017 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y17$residm17d17)
cvtabof[1,10] <- rmse(y17$heb_2017,y17$pred_m17d17) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m17d17_mresid_all.csv")
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate 2017 model with 2018 data: Herb_2018 ~ ARGmax1806 + NBRstd1803 + ARGmed1809
rfts <- randomForest(Herb_2018 ~ ARGmax1806 + NBRstd1803 + ARGmed1809 + ts2018, 
                     data = y18, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y18$pred_m17d18 <- predict(rfts)
y18$residm17d18 <- y18$pred_m17d18 - y18$Herb_2018
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2018 # data year
cvtabof[1,2] <- 2017 # model year
cvtabof[1,3] <- mean(y18$residm17d18[y18$ts2018 > 199])
cvtabof[1,4] <- mean(y18$residm17d18[y18$ts2018 < 200])
d <- lm(y18$residm17d18 ~ y18$ts2018)
# add ts2018 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y18$residm17d18[y18$ts2018 < 200] ~ y18$ts2018[y18$ts2018 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y18$residm17d18)
cvtabof[1,10] <- rmse(y18$Herb_2018,y18$pred_m17d18) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m17d18_mresid_all.csv")
write.csv(cvtabof,tname)


# OLDFIRES INCLUDE ALL ##### Generate rf residuals for model yr 2018 for all years & save tables as rfcv_**_mresid_all.csv ----------
# OLDFIRES INCLUDE ALL - Generate 2018 model with 2009 data: heb_2009 ~  ARGmax0906 + NBRstd1803 +  B3std1803
folder <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_20201104"
rfts <- randomForest(heb_2009 ~ ARGmax0906 + NBRstd0903 +  B3std0903 + ts2009, 
                     data = y09, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y09$pred_m18d09 <- predict(rfts)
y09$residm18d09 <- y09$pred_m18d09 - y09$heb_2009
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2009 # data year
cvtabof[1,2] <- 2018 # model year
cvtabof[1,3] <- mean(y09$residm18d09[y09$ts2009 > 199])
cvtabof[1,4] <- mean(y09$residm18d09[y09$ts2009 < 200])
d <- lm(y09$residm18d09 ~ y09$ts2009)
# add ts2009 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y09$residm18d09[y09$ts2009 < 200] ~ y09$ts2009[y09$ts2009 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y09$residm18d09)
cvtabof[1,10] <- rmse(y09$heb_2009,y09$pred_m18d09) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m18d09_mresid_all.csv")
write.csv(cvtabof,tname)

#OLDFIRES INCLUDE ALL - Generate 2018 model with 2011 data: heb_2011 ~ ARGmax1106 + NBRstd1103 + B3std1103
rfts <- randomForest(heb_2011 ~ ARGmax1106 + NBRstd1103 + B3std1103 + ts2011, 
                     data = y11, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y11$pred_m18d11 <- predict(rfts)
y11$residm18d11 <- y11$pred_m18d11 - y11$heb_2011
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2011 # data year
cvtabof[1,2] <- 2018 # model year
cvtabof[1,3] <- mean(y11$residm18d11[y11$ts2011 > 199])
cvtabof[1,4] <- mean(y11$residm18d11[y11$ts2011 < 200])
d <- lm(y11$residm18d11 ~ y11$ts2011)
# add ts2011 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y11$residm18d11[y11$ts2011 < 200] ~ y11$ts2011[y11$ts2011 < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y11$residm18d11)
cvtabof[1,10] <- rmse(y11$heb_2011,y11$pred_m18d11) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m18d11_mresid_all.csv")
write.csv(cvtabof,tname)


# OLDFIRES INCLUDE ALL - Generate 2018 model with 2013 data: heb_2013 ~ ARGmax1306 + NBRstd1303 + B3std1303
rfts <- randomForest(heb_2013 ~ ARGmax1306 + NBRstd1303 + B3std1303 + ts2013, 
                     data = y13, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y13$pred_m18d13 <- predict(rfts)
y13$residm18d13 <- y13$pred_m18d13 - y13$heb_2013
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2013 # data year
cvtabof[1,2] <- 2018 # model year
cvtabof[1,3] <- mean(y13$residm18d13[y13$ts2013 > 199])
cvtabof[1,4] <- mean(y13$residm18d13[y13$ts2013 < 200])
d <- lm(y13$residm18d13 ~ y13$ts2013)
# add ts2013 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y13$residm18d13[y13$ts2013 < 200] ~ y13$ts2013[y13$ts2013 < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y13$residm18d13)
cvtabof[1,10] <- rmse(y13$heb_2013,y13$pred_m18d13) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m18d13_mresid_all.csv")
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate 2018 model with 2015 data: heb_2015 ~ ARGmax1506 + NBRstd1503 + B3std1503
rfts <- randomForest(heb_2015 ~ ARGmax1506 + NBRstd1503 + ARGmed1509 + ts2015, 
                     data = y15, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y15$pred_m18d15 <- predict(rfts)
y15$residm18d15 <- y15$pred_m18d15 - y15$heb_2015
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2015 # data year
cvtabof[1,2] <- 2018 # model year
cvtabof[1,3] <- mean(y15$residm18d15[y15$ts2015 > 199])
cvtabof[1,4] <- mean(y15$residm18d15[y15$ts2015 < 200])
d <- lm(y15$residm18d15 ~ y15$ts2015)
# add ts2015 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y15$residm18d15[y15$ts2015 < 200] ~ y15$ts2015[y15$ts2015 < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y15$residm18d15)
cvtabof[1,10] <- rmse(y15$heb_2015,y15$pred_m18d15) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m18d15_mresid_all.csv")
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate 2018 model with 2017 data: heb_2017 ~  ARGmax1706 + NBRstd1703 + B3std1703
rfts <- randomForest(heb_2017 ~ ARGmax1706 + NBRstd1703 + B3std1703 + ts2017, 
                     data = y17, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y17$pred_m18d17 <- predict(rfts)
y17$residm18d17 <- y17$pred_m18d17 - y17$heb_2017
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2017 # data year
cvtabof[1,2] <- 2018 # model year
cvtabof[1,3] <- mean(y17$residm18d17[y17$ts2017 > 199])
cvtabof[1,4] <- mean(y17$residm18d17[y17$ts2017 < 200])
d <- lm(y17$residm18d17 ~ y17$ts2017)
# add ts2017 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y17$residm18d17[y17$ts2017 < 200] ~ y17$ts2017[y17$ts2017 < 200]) #lm for ONLY old fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y17$residm18d17)
cvtabof[1,10] <- rmse(y17$heb_2017,y17$pred_m18d17) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m18d17_mresid_all.csv")
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL - Generate 2018 model with 2018 data: Herb_2018 ~ ARGmax1806 + NBRstd1803 + B3std1803
rfts <- randomForest(Herb_2018 ~ ARGmax1806 + NBRstd1803 + ARGmed1809 + ts2018, 
                     data = y18, importance = TRUE, type = "regression", na.action = na.exclude)
a <- bptest(rfts, varformula = NULL, studentize = TRUE, data = list())  #------- Breusch Pagan Test for BP results & p-value
y18$pred_m18d18 <- predict(rfts)
y18$residm18d18 <- y18$pred_m18d18 - y18$Herb_2018
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2018 # data year
cvtabof[1,2] <- 2018 # model year
cvtabof[1,3] <- mean(y18$residm18d18[y18$ts2018 > 199]) # no past fires
cvtabof[1,4] <- mean(y18$residm18d18[y18$ts2018 < 200]) # yes past fires
d <- lm(y18$residm18d18 ~ y18$ts2018)
# add ts2018 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y18$residm18d18[y18$ts2018 < 200] ~ y18$ts2018[y18$ts2018 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y18$residm18d18)
cvtabof[1,10] <- rmse(y18$Herb_2018,y18$pred_m18d18) # rmse of residuals - all recs
cvtabof[1,11] <- a$statistic # BP test result
cvtabof[1,12] <- a$p.value # BP test p value
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All","rmseresid_A","BPresult","BPpvalue")
tname <- paste0(folder,"\\rfcv_m18d18_mresid_all.csv")
write.csv(cvtabof,tname)

# OLDFIRES INCLUDE ALL ##### Save data year files to csvs ----------
folder <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_20201104"
y09_l <- y09[,c(1,78,79,82,84,86,88,90,92)]
write.csv(y09_l,paste0(folder,"\\PastFires_pointsWithResids_dy2009_l.csv"))
y11_l <- y11[,c(1,78,79,82,84,86,88,90,92)]
write.csv(y11_l,paste0(folder,"\\PastFires_pointsWithResids_dy2011_l.csv"))
y13_l <- y13[,c(1,78,79,82,84,86,88,90,92)]
write.csv(y13_l,paste0(folder,"\\PastFires_pointsWithResids_dy2013_l.csv"))
y15_l <- y15[,c(1,78,79,82,84,86,88,90,92)]
write.csv(y15_l,paste0(folder,"\\PastFires_pointsWithResids_dy2015_l.csv"))
y17_l <- y17[,c(1,78,79,82,84,86,88,90,92)]
write.csv(y17_l,paste0(folder,"\\PastFires_pointsWithResids_dy2017_l.csv"))
y18_l <- y18[,c(1,78,79,82,84,86,88,90,92)]
write.csv(y18_l,paste0(folder,"\\PastFires_pointsWithResids_dy2018_l.csv"))

# OLDFIRES INCLUDE ALL ##### Append residual model results into a single file ------------------------------
library(purrr)
library(readr)
library(mgsub)
# substrRight <- function(x, n){
#   substr(x, nchar(x)-n+1, nchar(x))
# }
setwd("N:/project/monitoring_vol/SoCalShrubs/EE_Extractions/cvtables_20201104")
list_of_files <- list.files(path = ".",
                            pattern = "mresid_all",
                            full.names = FALSE)


all_models_tempint_all <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read_csv(.x, col_types = cols(), col_names = TRUE), .id = "file_name") # apply read_csv to all elements 


#c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All")
all_models_tempint_all$datayr_ch <-as.character(all_models_tempint_all$datayr)
all_models_tempint_all$modelyr_ch <-as.character(all_models_tempint_all$modelyr)

write.csv(all_models_tempint_all, file = "N:/project/monitoring_vol/SoCalShrubs/EE_Extractions/cvtables_20201104/all_sum_20201119.csv")


# OLDFIRES INCLUDE ALL ##### Heatmaps for for yes & no past fires -------------
setwd("N:/project/monitoring_vol/SoCalShrubs/EE_Extractions/cvtables_20201104")
rounded <-round(reshape2::acast(all_models_tempint_all,modelyr~datayr,fun=mean,value.var="mresid_All"),2)
means <- as.matrix(rowMeans(rounded))
#out <- data.frame(round_mresid_NoPF,means)
#pdf("mresid_All.pdf")
tiff("OLDFIRES_INCLUDE_ALL_mresid_All.tiff", units="in", width=6, height=6, res=300)
mresid_All_plot <- ggplot(data = all_models_tempint_all, mapping = aes(y = modelyr_ch, x = datayr_ch,fill = mresid_All)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(mresid_All,2))) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title="All recs regardless of fires prior to data year", subtitle=stringr::str_wrap("Cell values = mean of residuals",width = 65)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(rounded)[6],2), color="black",size=4) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(rounded)[5],2), color="black",size=4) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(rounded)[4],2), color="black",size=4) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(rounded)[3],2), color="black",size=4) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(rounded)[2],2), color="black",size=4) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(rounded)[1],2), color="black",size=4) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(mresid_All_plot)
dev.off()


# Mean residuals for NO past fires
round_mresid_NoPF <-round(reshape2::acast(all_models_tempint_all,modelyr~datayr,fun=mean,value.var="mresid_NoPF"),2)
#means <- as.matrix(rowMeans(round_mresid_NoPF))
#out <- data.frame(round_mresid_NoPF,means)
tiff("OLDFIRES_INCLUDE_ALL_mresid_NoPastFires.tiff", units="in", width=6, height=6, res=300)
#pdf("mresid_NoPF.pdf")
mresid_NoPF <- ggplot(data = all_models_tempint_all, mapping = aes(y = modelyr_ch, x = datayr_ch,fill = mresid_NoPF)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(mresid_NoPF,2))) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title="No fires prior to data year", subtitle=stringr::str_wrap("Cell values = mean of residuals",width = 65)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(round_mresid_NoPF)[6],2), color="black",size=4) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(round_mresid_NoPF)[5],2), color="black",size=4) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(round_mresid_NoPF)[4],2), color="black",size=4) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(round_mresid_NoPF)[3],2), color="black",size=4) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(round_mresid_NoPF)[2],2), color="black",size=4) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(round_mresid_NoPF)[1],2), color="black",size=4) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(mresid_NoPF)
dev.off()

# Mean residuals for Recs having fire prior to data year
#c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue")
round_mresid_YPF <-round(reshape2::acast(all_models_tempint_all,modelyr~datayr,fun=mean,value.var="mresid_YPF"),2)
#means <- as.matrix(rowMeans(round_mresid_YPF))
#out <- data.frame(round_mresid_YPF,means)
tiff("OLDFIRES_INCLUDE_ALL_mresid_YPF.tiff", units="in", width=6, height=6, res=300)
#pdf("mresid_YPF.pdf")
mresid_YPF <- ggplot(data = all_models_tempint_all, mapping = aes(y = modelyr_ch, x = datayr_ch,fill = mresid_YPF)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(mresid_YPF,2))) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title="Records with past fire prior to data year", subtitle=stringr::str_wrap("Cell values = mean of residuals",width = 65)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(round_mresid_YPF)[6],2), color="black",size=4) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(round_mresid_YPF)[5],2), color="black",size=4) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(round_mresid_YPF)[4],2), color="black",size=4) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(round_mresid_YPF)[3],2), color="black",size=4) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(round_mresid_YPF)[2],2), color="black",size=4) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(round_mresid_YPF)[1],2), color="black",size=4) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(mresid_YPF)
dev.off()

# All - both past fires & no past fires lm esimate
#c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue")
round_Alllm_estimate <-round(reshape2::acast(all_models_tempint_all,modelyr~datayr,fun=mean,value.var="Alllm_estimate"),3)
#means <- as.matrix(rowMeans(round_Alllm_estimate))
#out <- data.frame(round_Alllm_estimate,means)
tiff("OLDFIRES_INCLUDE_ALL_Alllm_estimate.tiff", units="in", width=6, height=6, res=300)
#pdf("Alllm_estimate.pdf")
Alllm_estimate <- ggplot(data = all_models_tempint_all, mapping = aes(y = modelyr_ch, x = datayr_ch,fill = Alllm_estimate)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(Alllm_estimate,3)),size=3.5) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title="All - both past fires & no past fires", subtitle=stringr::str_wrap("Cell values = lm estimate",width = 65)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(round_Alllm_estimate)[6],3), color="black",size=3.3) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(round_Alllm_estimate)[5],3), color="black",size=3.3) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(round_Alllm_estimate)[4],3), color="black",size=3.3) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(round_Alllm_estimate)[3],3), color="black",size=3.3) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(round_Alllm_estimate)[2],3), color="black",size=3.3) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(round_Alllm_estimate)[1],3), color="black",size=3.3) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(Alllm_estimate)
dev.off()

# All - both past fires & no past fires p-value
#c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue")
rounded <-round(reshape2::acast(all_models_tempint_all,modelyr~datayr,fun=mean,value.var="Allpvalue"),3)
#pdf("Allpvalue.pdf")
tiff("OLDFIRES_INCLUDE_ALL_Allpvalue.tiff", units="in", width=6, height=6, res=300)
Allpvalue <- ggplot(data = all_models_tempint_all, mapping = aes(y = modelyr_ch, x = datayr_ch,fill = Allpvalue)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(Allpvalue,3)),size=3.5) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title="All - regardless of past fires", subtitle=stringr::str_wrap("Cell values = p-value",width = 65)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(rounded)[6],2), color="black",size=3.3) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(rounded)[5],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(rounded)[4],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(rounded)[3],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(rounded)[2],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(rounded)[1],2), color="black",size=3.3) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(Allpvalue)
dev.off()


# YPFlm_estimate
#c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue")
rounded <-round(reshape2::acast(all_models_tempint_all,modelyr~datayr,fun=mean,value.var="YPFlm_estimate"),3)
#means <- as.matrix(rowMeans(rounded))
#out <- data.frame(round_Aoflm_estimate,means)
tiff("OLDFIRES_INCLUDE_ALL_YPFlm_estimate.tiff", units="in", width=6, height=6, res=300)
#pdf("YPFlm_estimate.pdf")
YPFlm_estimateplot <- ggplot(data = all_models_tempint_all, mapping = aes(y = modelyr_ch, x = datayr_ch,fill = YPFlm_estimate)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(YPFlm_estimate,3)),size=3.5) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title= "Records with past fire prior to data year", subtitle=stringr::str_wrap("Cell values = lm estimate",width = 65)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(rounded)[6],3), color="black",size=3.3) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(rounded)[5],3), color="black",size=3.3) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(rounded)[4],3), color="black",size=3.3) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(rounded)[3],3), color="black",size=3.3) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(rounded)[2],3), color="black",size=3.3) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(rounded)[1],3), color="black",size=3.3) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(YPFlm_estimateplot)
dev.off()

# ONLY past fires Yofpvalue
#c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue")
rounded <-round(reshape2::acast(all_models_tempint_all,modelyr~datayr,fun=mean,value.var="YPFpvalue"),3)
#means <- as.matrix(rowMeans(rounded))
#out <- data.frame(round_Aoflm_estimate,means)
tiff("OLDFIRES_INCLUDE_ALL_YPFpvalue.tiff", units="in", width=6, height=6, res=300)
#pdf("YPFpvalue.pdf")
YPFpvalue <- ggplot(data = all_models_tempint_all, mapping = aes(y = modelyr_ch, x = datayr_ch,fill = YPFpvalue)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(YPFpvalue,3)),size=3.5) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title= "Records with past fire prior to data year", subtitle=stringr::str_wrap("Cell values = p value",width = 65)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(rounded)[6],2), color="black",size=3.3) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(rounded)[5],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(rounded)[4],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(rounded)[3],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(rounded)[2],2), color="black",size=3.3) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(rounded)[1],2), color="black",size=3.3) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(YPFpvalue)
dev.off()

# OLDFIRES INCLUDE ALL ##### Generate a cross-validation table for each rf run - rerun 20201117 to include time since -------------------
# ### use y09,y11, etc - data frames with no past fires removed
# ### create files called: paste0("stats2_", outstr,"_iter100.csv")

setwd("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_20201104")

intable <- y09
listofstrings <- c("heb_2009 ~ ARGmax0903 + NBRstd0903 + ARGmed0906 + ts2009",
                   "heb_2009 ~ NDVstd0903 + ARGmea0906 + B3std0906 + ts2009",
                   "heb_2009 ~ ARGmea0906 + NDVstd0903 + B3max0909 + ts2009",
                   "heb_2009 ~ ARGmax0903 + ARGmax0906 + NDVmin0909 + ts2009",
                   "heb_2009 ~ ARGmax0906 + NBRstd0903 + ARGmed0909 + ts2009",
                   "heb_2009 ~ ARGmax0906 + NBRstd0903 + B3std0903 + ts2009")

intable <- y11
listofstrings <- c("heb_2011 ~ ARGmax1103 + NBRstd1103 + ARGmed1106 + ts2011",
                   "heb_2011 ~ NDVstd1103 + ARGmea1106 + B3std1106 + ts2011",
                   "heb_2011 ~ ARGmea1106 + NDVstd1103 + B3max1109 + ts2011",
                   "heb_2011 ~ ARGmax1103 + ARGmax1106 + NDVmin1109 + ts2011",
                   "heb_2011 ~ ARGmax1106 + NBRstd1103 + ARGmed1109 + ts2011",
                   "heb_2011 ~ ARGmax1106 + NBRstd1103 + B3std1103 + ts2011")

intable <- y13
listofstrings <- c("heb_2013 ~ ARGmax1303 + NBRstd1303 + ARGmed1306 + ts2013",
                   "heb_2013 ~ NDVstd1303 + ARGmea1306 + B3std1306 + ts2013",
                   "heb_2013 ~ ARGmea1306 + NDVstd1303 + B3max1309 + ts2013",
                   "heb_2013 ~ ARGmax1303 + ARGmax1306 + NDVmin1309 + ts2013",
                   "heb_2013 ~ ARGmax1306 + NBRstd1303 + ARGmed1309 + ts2013",
                   "heb_2013 ~ ARGmax1306 + NBRstd1303 + B3std1303 + ts2013")

intable <- y15
listofstrings <- c("heb_2015 ~ ARGmax1503 + NBRstd1503 + ARGmed1506 + ts2015",
                   "heb_2015 ~ NDVstd1503 + ARGmea1506 + B3std1506 + ts2015",
                   "heb_2015 ~ ARGmea1506 + NDVstd1503 + B3max1509 + ts2015",
                   "heb_2015 ~ ARGmax1503 + ARGmax1506 + NDVmin1509 + ts2015",
                   "heb_2015 ~ ARGmax1506 + NBRstd1503 + ARGmed1509 + ts2015",
                   "heb_2015 ~ ARGmax1506 + NBRstd1503 + B3std1503 + ts2015")

intable <- y17
listofstrings <- c("heb_2017 ~ ARGmax1703 + NBRstd1703 + ARGmed1706 + ts2017",
                   "heb_2017 ~ NDVstd1703 + ARGmea1706 + B3std1706 + ts2017",
                   "heb_2017 ~ ARGmea1706 + NDVstd1703 + B3max1709 + ts2017",
                   "heb_2017 ~ ARGmax1703 + ARGmax1706 + NDVmin1709 + ts2017",
                   "heb_2017 ~ ARGmax1706 + NBRstd1703 + ARGmed1709 + ts2017",
                   "heb_2017 ~ ARGmax1706 + NBRstd1703 + B3std1703 + ts2017") 

intable <- y18
listofstrings <- c("Herb_2018 ~ ARGmax1803 + NBRstd1803 + ARGmed1806 + ts2018",
                   "Herb_2018 ~ NDVstd1803 + ARGmea1806 + B3std1806 + ts2018",
                   "Herb_2018 ~ ARGmea1806 + NDVstd1803 + B3max1809 + ts2018",
                   "Herb_2018 ~ ARGmax1803 + ARGmax1806 + NDVmin1809 + ts2018",
                   "Herb_2018 ~ ARGmax1806 + NBRstd1803 + ARGmed1809 + ts2018",
                   "Herb_2018 ~ ARGmax1806 + NBRstd1803 + B3std1803 + ts2018")

# Next time running this add the model year to the cvtab
setwd("N:/project/monitoring_vol/SoCalShrubs/EE_Extractions/cvtables_20201104")
for (instring in listofstrings){
  rf <- randomForest(as.formula(instring),data = intable, importance = TRUE, type = "regression", na.action = na.exclude)
  outstr <- mgsub::mgsub(instring, c("\\+","\\,", "\\~","\\ ","\\(","\\)"), c("_","","_","","","")); outstr
  cvtab <- data.frame()
  for (i in 1:5){ # i is p_withheld
    cvout <- rf.crossValidation(rf, intable, p=i/10, n=100, ntree=501)
    cvout
    cvtab[i,1] <- cvout$fit.var.exp
    cvtab[i,2] <- cvout$fit.mse
    cvtab[i,3] <- median(cvout$y.rmse)
    cvtab[i,4] <- var(cvout$y.rmse)
    cvtab[i,5] <- median(cvout$y.mbe)
    cvtab[i,6] <- median(cvout$y.mae)
    cvtab[i,7] <- median(cvout$model.mse)
    cvtab[i,8] <- mean(cvout$y.rmse)
    colnames(cvtab) <- c("fit.var.exp","fit.mse","med.y.rmse","var.y.rmse","med.y.mbe","med.y.mae","med.model.mse","mean.y.rmse")
    
    tname <- paste0("stats2_", outstr,"_iter100.csv")
    write.csv(cvtab,tname,col.names=TRUE)
  }
}


# OLDFIRES INCLUDE ALL ##### Append cross validation model results into one file and set up for graphing the combination table --------------------
require(reshape2)
library(reshape)
require(ggplot2)
require(scales)
require(plot3D)
require(akima) #Interpolation of Irregularly and Regularly Spaced Data
require(rgl) # 3D Visualization Using OpenGL
library(manipulate)
library(data.table)
library(purrr)
library(readr)
library(mgsub)
setwd("N:/project/monitoring_vol/SoCalShrubs/EE_Extractions/cvtables_20201104")
list_of_files <- list.files(path = ".",
                            pattern = "stats2",
                            full.names = FALSE)
all_models_cv <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read_csv(.x, col_types = cols(), col_names = TRUE), .id = "file_name")
all_models_cv$year <- gsub("_", "",substring(as.character(all_models_cv$file_name), 12, 16))
all_models_cv$var1 <- gsub("_", "",substring(as.character(all_models_cv$file_name), 17, 27))
all_models_cv$var2 <- gsub("_", "",substring(as.character(all_models_cv$file_name), 28, 38))
all_models_cv$var3 <- gsub("_", "",substring(as.character(all_models_cv$file_name), 39, 49))
all_models_cv$iter <- gsub("\\.", "",stringr::str_sub(as.character(all_models_cv$file_name), -7,-4)) # all 100 at this point

# Create BaseModel field with no data year info
all_models_cv$BaseModel <- mgsub::mgsub(all_models_cv$file_name,c("stats2_Herb_2018_","stats2_heb_2017_","stats2_heb_2015_","stats2_heb_2013_","stats2_heb_2011_","stats2_heb_2009_","_iter100.csv"),c("","","","","","",""))
all_models_cv$BaseModel <- mgsub::mgsub(all_models_cv$BaseModel, c("0903","1103","1303","1503","1703","1803"),c("YR03","YR03","YR03","YR03","YR03","YR03"))
all_models_cv$BaseModel <- mgsub::mgsub(all_models_cv$BaseModel, c("0906","1106","1306","1506","1706","1806"),c("YR06","YR06","YR06","YR06","YR06","YR06"))
all_models_cv$BaseModel <- mgsub::mgsub(all_models_cv$BaseModel, c("0909","1109","1309","1509","1709","1809"),c("YR09","YR09","YR09","YR09","YR09","YR09"))
unique(all_models_cv$BaseModel)

# cbind is used to populate IndVars & OrigModelYr to above allrfcv* file
# 20201026 join Indvars & OrigmodelYr from previous table (N:\project\monitoring_vol\SoCalShrubs\EE_Extractions\cvtables5_100_3var\allrfcv_20200528.csv) 
# 20201102 - forget the xls work & just cbind here
# 20201117 - added field for meanyrmse
oldtab <- read.csv("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables5_100_3var\\allrfcv_20200528_needinscript.csv")
all_models_cv <- cbind(all_models_cv,oldtab[,3:4])
head(all_models_cv)

all_models_cv$OrigModelYr_ch <- as.character(all_models_cv$OrigModelYr)
all_models_cv$DataYear_ch <- as.character(all_models_cv$year)
colnames(all_models_cv) <- c("file_name","p_withheld", "fitvarexp","fitmse","medyrmse","varyrmse","medymbe","medymae","medmodelmse","meanyrmse","DataYear","var1","var2","var3","iter","BaseModel","IndVars" ,"OrigModelYr","OrigModelYr_ch","DataYear_ch")
all_models_cv$p_withheld <- all_models_cv$p_withheld * 10

head(all_models_cv)
write.csv(all_models_cv, file = "allrfcv_20201117.csv")

uniquebm <- unique(all_models_cv$BaseModel)
uniquedatayr <- unique(all_models_cv$DataYear)
#all_models_melt <- melt(all_models_2more, id.var = c("IndVars","OrigModelYr","p_withheld"))
#head(all_models_melt)


# OLDFIRES INCLUDE ALL ##### Calculate med & mean y rmse at percent withheld 10%, 30%, 50% for each ModelYr-DataYear combo AND ---------
# Calculate ratio of medyrmse at percent withheld 10%:50% for each ModelYr-DataYear combo
# 20201117 Add fields to hold 30% medyrmse & 10% & 30% mean y rmse
setwd("N:/project/monitoring_vol/SoCalShrubs/EE_Extractions/cvtables_20201104")
all_models_cv$medyrmseratio <- 0; all_models_cv$medyrmsepw10 <- 0; all_models_cv$medyrmsepw50 <- 0
for (eachbm in uniquebm){
  print(eachbm)
  for (eachdy in uniquedatayr) {
    print(uniquedatayr)
    ## oiginal code to just calculate the ratio
    # print(eachdy)
    # fitvar_pw10 <- all_models_cv[all_models_cv$p_withheld == 10 & all_models_cv$BaseModel == eachbm & all_models_cv$DataYear == eachdy,"medyrmse"]
    # print(fitvar_pw10)
    # fitvar_pw50 <- all_models_cv[all_models_cv$p_withheld == 50 & all_models_cv$BaseModel == eachbm & all_models_cv$DataYear == eachdy,"medyrmse"]
    # print(fitvar_pw50)
    # medyrmseratio <- fitvar_pw10/fitvar_pw50
    ### New section where fitvar_pw10 & fitvar_pw50 are saved in table
    medyrmsepw10 <- all_models_cv[all_models_cv$p_withheld == 10 & all_models_cv$BaseModel == eachbm & all_models_cv$DataYear == eachdy,"medyrmse"]
    print(medyrmsepw10)
    medyrmsepw50 <- all_models_cv[all_models_cv$p_withheld == 50 & all_models_cv$BaseModel == eachbm & all_models_cv$DataYear == eachdy,"medyrmse"]
    print(medyrmsepw50)
    medyrmsepw30 <- all_models_cv[all_models_cv$p_withheld == 30 & all_models_cv$BaseModel == eachbm & all_models_cv$DataYear == eachdy,"medyrmse"]
    print(medyrmsepw30)
    meanyrmsepw10 <- all_models_cv[all_models_cv$p_withheld == 10 & all_models_cv$BaseModel == eachbm & all_models_cv$DataYear == eachdy,"meanyrmse"]
    print(meanyrmsepw10)
    meanyrmsepw30 <- all_models_cv[all_models_cv$p_withheld == 30 & all_models_cv$BaseModel == eachbm & all_models_cv$DataYear == eachdy,"meanyrmse"]
    print(meanyrmsepw30)
    
    ###calc medyrmse ratio 10/50
    medyrmseratio <- medyrmsepw10/medyrmsepw50
    medyrmseratio_ch <- paste(medyrmseratio, ":1",sep='')
    print(medyrmseratio_ch)
    # Add to table
    all_models_cv$medyrmseratio[all_models_cv$BaseModel == eachbm & all_models_cv$DataYear == eachdy] <- medyrmseratio
    all_models_cv$medyrmsepw10[all_models_cv$BaseModel == eachbm & all_models_cv$DataYear == eachdy] <- medyrmsepw10
    all_models_cv$medyrmsepw50[all_models_cv$BaseModel == eachbm & all_models_cv$DataYear == eachdy] <- medyrmsepw50
    all_models_cv$medyrmsepw30[all_models_cv$BaseModel == eachbm & all_models_cv$DataYear == eachdy] <- medyrmsepw30
    all_models_cv$meanyrmsepw10[all_models_cv$BaseModel == eachbm & all_models_cv$DataYear == eachdy] <- meanyrmsepw10
    all_models_cv$meanyrmsepw30[all_models_cv$BaseModel == eachbm & all_models_cv$DataYear == eachdy] <- meanyrmsepw30

    #all_models_2more$medyrmseratio <- ifelse(all_models_2more$p_withheld == 10 & all_models_2more$BaseModel == eachbm & all_models_2more$DataYear == eachdy, medyrmseratio, all_models_2more$medyrmseratio)
  }
}
unique(all_models_cv$medyrmseratio)
unique(all_models_cv$medyrmsepw10)
unique(all_models_cv$medyrmsepw50)
unique(all_models_cv$meanyrmse)
head(all_models_cv); tail(all_models_cv)
write.csv(all_models_cv, file = "all_models_cv.csv")
all_models_cv <- read.table("all_models_cv.csv", header=T,sep=',',stringsAsFactors = FALSE)


# OLDFIRES INCLUDE ALL ##### Plot for only 10%  (medyrmsepw10) - median Y RMSE where 10% is withheld -----------
setwd("N:/project/monitoring_vol/SoCalShrubs/EE_Extractions/cvtables_20201104")
rounded <-round(reshape2::acast(all_models_cv,OrigModelYr~DataYear,fun=mean,value.var="medyrmsepw10"),2)
#means <- as.matrix(rowMeans(z_medyrmseratio))
#out <- data.frame(z_medyrmseratio,means)
tiff("OLDFIRES_INCLUDE_ALL_medyrmsepw10.tiff", units="in", width=6, height=6, res=300)
medyrmsepw10_plot <- ggplot(data = all_models_cv, mapping = aes(y = OrigModelYr_ch, x = DataYear_ch,fill = medyrmsepw10)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(medyrmsepw10,2))) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title="Model performance as a function of training data withheld", subtitle=stringr::str_wrap("Cell values =  med y rmse where 10% is withheld",width = 65)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(rounded)[6],2), color="black",size=4) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(rounded)[5],2), color="black",size=4) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(rounded)[4],2), color="black",size=4) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(rounded)[3],2), color="black",size=4) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(rounded)[2],2), color="black",size=4) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(rounded)[1],2), color="black",size=4) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(medyrmsepw10_plot)
dev.off()

# OLDFIRES INCLUDE ALL ##### Plot for only 30%  (medyrmsepw10) - median Y RMSE where 10% is withheld -----------
setwd("N:/project/monitoring_vol/SoCalShrubs/EE_Extractions/cvtables_20201104")
rounded <-round(reshape2::acast(all_models_cv,OrigModelYr~DataYear,fun=mean,value.var="medyrmsepw30"),2)
#means <- as.matrix(rowMeans(z_medyrmseratio))
#out <- data.frame(z_medyrmseratio,means)
tiff("OLDFIRES_INCLUDE_ALL_medyrmsepw30.tiff", units="in", width=6, height=6, res=300)
plot1 <- ggplot(data = all_models_cv, mapping = aes(y = OrigModelYr_ch, x = DataYear_ch,fill = medyrmsepw10)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(medyrmsepw10,2))) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title="Model performance as a function of training data withheld", subtitle=stringr::str_wrap("Cell values =  med y rmse where 30% is withheld",width = 65)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(rounded)[6],2), color="black",size=4) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(rounded)[5],2), color="black",size=4) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(rounded)[4],2), color="black",size=4) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(rounded)[3],2), color="black",size=4) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(rounded)[2],2), color="black",size=4) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(rounded)[1],2), color="black",size=4) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(plot1)
dev.off()

# OLDFIRES INCLUDE ALL ##### Plot for only 10%  (medyrmsepw10) - median Y RMSE where 10% is withheld -----------
setwd("N:/project/monitoring_vol/SoCalShrubs/EE_Extractions/cvtables_20201104")
rounded <-round(reshape2::acast(all_models_cv,OrigModelYr~DataYear,fun=mean,value.var="meanyrmsepw10"),2)
#means <- as.matrix(rowMeans(z_medyrmseratio))
#out <- data.frame(z_medyrmseratio,means)
tiff("OLDFIRES_INCLUDE_ALL_meanyrmsepw10.tiff", units="in", width=6, height=6, res=300)
plot1 <- ggplot(data = all_models_cv, mapping = aes(y = OrigModelYr_ch, x = DataYear_ch,fill = meanyrmsepw10)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(meanyrmsepw10,2))) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title="", subtitle=stringr::str_wrap("Cell values =  MEAN y rmse where 10% is withheld",width = 65)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(rounded)[6],2), color="black",size=4) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(rounded)[5],2), color="black",size=4) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(rounded)[4],2), color="black",size=4) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(rounded)[3],2), color="black",size=4) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(rounded)[2],2), color="black",size=4) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(rounded)[1],2), color="black",size=4) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(plot1)
dev.off()

# OLDFIRES INCLUDE ALL ##### Plot for only 30%  (medyrmsepw10) - median Y RMSE where 10% is withheld -----------
setwd("N:/project/monitoring_vol/SoCalShrubs/EE_Extractions/cvtables_20201104")
rounded <-round(reshape2::acast(all_models_cv,OrigModelYr~DataYear,fun=mean,value.var="meanyrmsepw30"),2)
#means <- as.matrix(rowMeans(z_medyrmseratio))
#out <- data.frame(z_medyrmseratio,means)
tiff("OLDFIRES_INCLUDE_ALL_meanyrmsepw30.tiff", units="in", width=6, height=6, res=300)
plot1 <- ggplot(data = all_models_cv, mapping = aes(y = OrigModelYr_ch, x = DataYear_ch,fill = meanyrmsepw30)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(meanyrmsepw30,2))) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title="", subtitle=stringr::str_wrap("Cell values =  MEAN y rmse where 30% is withheld",width = 65)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(rounded)[6],2), color="black",size=4) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(rounded)[5],2), color="black",size=4) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(rounded)[4],2), color="black",size=4) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(rounded)[3],2), color="black",size=4) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(rounded)[2],2), color="black",size=4) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(rounded)[1],2), color="black",size=4) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(plot1)
dev.off()



# OLDFIRES INCLUDE ALL ##### Heatmap for 10 - 50 median Y RMSE ratio -------------
z_medyrmseratio <-round(reshape2::acast(all_models_cv,OrigModelYr~DataYear,fun=mean,value.var="medyrmseratio"),2)
means <- as.matrix(rowMeans(z_medyrmseratio))
out <- data.frame(z_medyrmseratio,means)
tiff("OLDFIRES_INCLUDE_ALL_rmse_ratio.tiff", units="in", width=6, height=6, res=300)
rmse_ratio <- ggplot(data = all_models_cv, mapping = aes(y = OrigModelYr_ch, x = DataYear_ch,fill = medyrmseratio)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(medyrmseratio,2))) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title="Model performance as a function of training data withheld", subtitle=stringr::str_wrap("Cell values = Ratio of median rmse for 100 iterations of 10% training data withheld / 50% withheld",width = 65)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.85, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(z_medyrmseratio)[6],2), color="black",size=4) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(z_medyrmseratio)[5],2), color="black",size=4) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(z_medyrmseratio)[4],2), color="black",size=4) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(z_medyrmseratio)[3],2), color="black",size=4) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(z_medyrmseratio)[2],2), color="black",size=4) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(z_medyrmseratio)[1],2), color="black",size=4) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1)) +
  theme(plot.subtitle=element_text(size=13, hjust=0.5, vjust=0.5, face="italic", colour="black"))
print(rmse_ratio)
dev.off()




# OLDFIRES INCLUDE ALL ##### Heatmap  for fitvarexp - fir variance explained -----------

z_fitvarexp <-round(reshape2::acast(all_models_cv,OrigModelYr~DataYear,fun=mean,value.var="fitvarexp"),2)
fitmeans <- as.matrix(rowMeans(z_fitvarexp))
fitout <- data.frame(z_fitvarexp,fitmeans)
#pdf("OLDFIRES_INCLUDE_ALL_fitvarexp.pdf")
tiff("OLDFIRES_INCLUDE_ALL_fitvarexp.tiff", units="in", width=6, height=6, res=300)
plot1 <- ggplot(data = all_models_cv, mapping = aes(y = OrigModelYr_ch, x = DataYear_ch,fill = fitvarexp)) +
  geom_tile(aes(width=1,height=1),show.legend = TRUE) + geom_text(aes(label=round(fitvarexp,0))) + 
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x="Data Year", y="Base Model Year",title="Model % Variance Explained") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(legend.title = element_text(size=10, hjust=0, face=1, colour="black"),legend.justification = "center") + labs(fill = "% Variance Explained") +
  theme(legend.position = "bottom") + 
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  annotate("text", x=6.8, y=6.9, label="Base Model", color="black",size=3.5) +
  annotate("text", x=6.86, y=6.62, label="Means", color="black",size=3.5) +
  annotate("text", x=6.8, y=6, label=round(rowMeans(z_fitvarexp)[6],0), color="black",size=4) + 
  annotate("text", x=6.8, y=5, label=round(rowMeans(z_fitvarexp)[5],0), color="black",size=4) +
  annotate("text", x=6.8, y=4, label=round(rowMeans(z_fitvarexp)[4],0), color="black",size=4) +
  annotate("text", x=6.8, y=3, label=round(rowMeans(z_fitvarexp)[3],0), color="black",size=4) +
  annotate("text", x=6.8, y=2, label=round(rowMeans(z_fitvarexp)[2],0), color="black",size=4) +
  annotate("text", x=6.8, y=1, label=round(rowMeans(z_fitvarexp)[1],00), color="black",size=4) +
  coord_cartesian(clip = "off") + 
  theme(plot.title=element_text(size=15, hjust=0.5, face=1, colour="black", vjust=-1))
print(plot1)
dev.off()




# OLDFIRES INCLUDE ALL ##### Plot panel of Observed vs predicted  ------
setwd("N:/project/monitoring_vol/SoCalShrubs/EE_Extractions/cvtables_20201104")
rfts <- randomForest(heb_2009 ~ ARGmax0903 + NBRstd0903 + ARGmed0906 + ts2009, 
                     data = y09, importance = TRUE, type = "regression", na.action = na.exclude)
y09$pred_m09d09 <- predict(rfts)

rfts <- randomForest(heb_2011 ~ ARGmax1103 + NBRstd1103 + ARGmed1106 + ts2011, 
                     data = y11, importance = TRUE, type = "regression", na.action = na.exclude)
y11$pred_m09d11 <- predict(rfts)

rfts <- randomForest(heb_2013 ~ ARGmax1303 + NBRstd1303 + ARGmed1306 + ts2013, 
                     data = y13, importance = TRUE, type = "regression", na.action = na.exclude)
y13$pred_m09d13 <- predict(rfts)

rfts <- randomForest(heb_2015 ~ ARGmax1503 + NBRstd1503 + ARGmed1506 + ts2015, 
                     data = y15, importance = TRUE, type = "regression", na.action = na.exclude)
y15$pred_m09d15 <- predict(rfts)

rfts <- randomForest(heb_2017 ~ ARGmax1703 + NBRstd1703 + ARGmed1706 + ts2017, 
                     data = y17, importance = TRUE, type = "regression", na.action = na.exclude)
y17$pred_m09d17 <- predict(rfts)

rfts <- randomForest(Herb_2018 ~ ARGmax1803 + NBRstd1803 + ARGmed1806 + ts2018, 
                     data = y18, importance = TRUE, type = "regression", na.action = na.exclude)
y18$pred_m09d18 <- predict(rfts)
library(gridExtra)
library(cowplot)

d09 <- ggplot(y09, aes(x = heb_2009, y = pred_m09d09)) +
  #geom_point() + 
  geom_hex() +
  scale_fill_gradient(low = "yellow", high = "red") + #, trans = "log10") +
  geom_point(shape = 16) +
  stat_smooth(method=lm) +
  geom_segment(aes(x=0,y=0,xend=100,yend=100,colour='red'),size=1) + #1:1 line
  theme(legend.position = "none") +
  xlab(label= "") + 
  ylab(label= "") +
  annotate("text", label = "2009", size=6, x = 45, y = 94) +
  #axis.text = element_blank() +
  guides(x = "none", y = "none") + # remove the x & y axis numbers & ticks
  #ggtitle("2009") +
  #theme(plot.margin = unit(c(0, -10, -10, -10), "pt"))
  #theme(axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black")) +
  theme(plot.title = element_text(margin = ggplot2::margin(t = 0, r=-10, b = -10, l=-10))) # top, right, bottom, left
    #plot.title = element_text(color="blue", size=14, face=1,hjust = 0.5),
    #axis.text = element_blank()
    #axis.title.x = element_text(color="blue", size=14, face="bold"),
    #axis.title.y = element_text(color="blue", size=14, face="bold")
plot(d09)

d11 <- ggplot(y11, aes(x = heb_2011, y = pred_m09d11)) +
  #geom_point() + 
  geom_hex() +
  scale_fill_gradient(low = "yellow", high = "red") + #, trans = "log10") +
  geom_point(shape = 16) +
  stat_smooth(method=lm) +
  geom_segment(aes(x=0,y=0,xend=100,yend=100,colour='red'),size=1) + #1:1 line
  theme(legend.position = "none") +
  xlab(label= "") + 
  ylab(label= "") +
  annotate("text", label = "2011", size=6, x = 45, y = 94) +
  #axis.text = element_blank() +
  guides(x = "none", y = "none") + # remove the x & y axis numbers & ticks
  #ggtitle("2009") +
  #theme(plot.margin = unit(c(0, -10, -10, -10), "pt"))
  #theme(axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black")) +
  theme(plot.title = element_text(margin = ggplot2::margin(t = 0, r=-10, b = -10, l=-10))) # top, right, bottom, left
#plot.title = element_text(color="blue", size=14, face=1,hjust = 0.5),
#axis.text = element_blank()
#axis.title.x = element_text(color="blue", size=14, face="bold"),
#axis.title.y = element_text(color="blue", size=14, face="bold")

d13 <- ggplot(y13, aes(x = heb_2013, y = pred_m09d13)) +
  #geom_point() + 
  geom_hex() +
  scale_fill_gradient(low = "yellow", high = "red") + #, trans = "log10") +
  geom_point(shape = 16) +
  stat_smooth(method=lm) +
  geom_segment(aes(x=0,y=0,xend=100,yend=100,colour='red'),size=1) + #1:1 line
  theme(legend.position = "none") +
  xlab(label= "") + 
  ylab(label= "") +
  annotate("text", label = "2013", size=6, x = 45, y = 94) +
  #axis.text = element_blank() +
  guides(x = "none", y = "none") + # remove the x & y axis numbers & ticks
  #ggtitle("2009") +
  #theme(plot.margin = unit(c(0, -10, -10, -10), "pt"))
  #theme(axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black")) +
  theme(plot.title = element_text(margin = ggplot2::margin(t = 0, r=-10, b = -10, l=-10))) # top, right, bottom, left
#plot.title = element_text(color="blue", size=14, face=1,hjust = 0.5),
#axis.text = element_blank()
#axis.title.x = element_text(color="blue", size=14, face="bold"),
#axis.title.y = element_text(color="blue", size=14, face="bold")

d15 <- ggplot(y15, aes(x = heb_2015, y = pred_m09d15)) +
  #geom_point() + 
  geom_hex() +
  scale_fill_gradient(low = "yellow", high = "red") + #, trans = "log10") +
  geom_point(shape = 16) +
  stat_smooth(method=lm) +
  geom_segment(aes(x=0,y=0,xend=100,yend=100,colour='red'),size=1) + #1:1 line
  theme(legend.position = "none") +
  xlab(label= "") + 
  ylab(label= "") +
  annotate("text", label = "2015", size=6, x = 45, y = 94) +
  #axis.text = element_blank() +
  guides(x = "none", y = "none") + # remove the x & y axis numbers & ticks
  #ggtitle("2009") +
  #theme(plot.margin = unit(c(0, -10, -10, -10), "pt"))
  #theme(axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black")) +
  theme(plot.title = element_text(margin = ggplot2::margin(t = 0, r=-10, b = -10, l=-10))) # top, right, bottom, left
#plot.title = element_text(color="blue", size=14, face=1,hjust = 0.5),
#axis.text = element_blank()
#axis.title.x = element_text(color="blue", size=14, face="bold"),
#axis.title.y = element_text(color="blue", size=14, face="bold")

d17 <- ggplot(y17, aes(x = heb_2017, y = pred_m09d17)) +
  #geom_point() + 
  geom_hex() +
  scale_fill_gradient(low = "yellow", high = "red") + #, trans = "log10") +
  geom_point(shape = 16) +
  stat_smooth(method=lm) +
  geom_segment(aes(x=0,y=0,xend=100,yend=100,colour='red'),size=1) + #1:1 line
  theme(legend.position = "none") +
  xlab(label= "") + 
  ylab(label= "") +
  annotate("text", label = "2017", size=6, x = 45, y = 94) +
  #axis.text = element_blank() +
  guides(x = "none", y = "none") + # remove the x & y axis numbers & ticks
  #ggtitle("2009") +
  #theme(plot.margin = unit(c(0, -10, -10, -10), "pt"))
  #theme(axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black")) +
  theme(plot.title = element_text(margin = ggplot2::margin(t = 0, r=-10, b = -10, l=-10))) # top, right, bottom, left
#plot.title = element_text(color="blue", size=14, face=1,hjust = 0.5),
#axis.text = element_blank()
#axis.title.x = element_text(color="blue", size=14, face="bold"),
#axis.title.y = element_text(color="blue", size=14, face="bold")

d18 <- ggplot(y18, aes(x = Herb_2018, y = pred_m09d18)) +
  #geom_point() + 
  geom_hex() +
  scale_fill_gradient(low = "yellow", high = "red") + #, trans = "log10") +
  geom_point(shape = 16,size=1) +
  stat_smooth(method=lm) +
  geom_segment(aes(x=0,y=0,xend=100,yend=100,colour='red'),size=1) + #1:1 line
  theme(legend.position = "none") +
  xlab(label= "") + 
  ylab(label= "") +
  annotate("text", label = "2018", size=4, x = 45, y = 94) +
  #axis.text = element_blank() +
  guides(x = "none", y = "none") + # remove the x & y axis numbers & ticks
  #ggtitle("2009") +
  #theme(plot.margin = unit(c(0, -10, -10, -10), "pt"))
  #theme(axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black")) +
  theme(plot.title = element_text(margin = ggplot2::margin(t = 0, r=-10, b = -10, l=-10))) # top, right, bottom, left
#plot.title = element_text(color="blue", size=14, face=1,hjust = 0.5),
#axis.text = element_blank()
#axis.title.x = element_text(color="blue", size=14, face="bold"),
#axis.title.y = element_text(color="blue", size=14, face="bold")
# topt = textGrob(
#   "2009 Model Applied to Each Year",
#   gp = gpar(fontface = 1, fontsize = 17,col='blue')
#   #hjust = .9,
#   #x = 1
# )
bottomt = textGrob(
  "Observed Herbaceous Cover (%)",
  gp = gpar(fontface = 1, fontsize = 10,col='blue'),
  vjust = -.5
  #x = 1
)

leftt = textGrob(
  "Predicted Herbaceous Cover (%)",
  gp = gpar(fontface = 1, fontsize = 10,col='blue'),rot = 90,
  vjust = 1.4
  #x = 1
)
tiff("OLDFIRES_INCLUDE_ALL_2009PlotPanel.tiff", units="in", width=5, height=3, res=300)
plot1 <- grid.arrange(d09,d11,d13,d15,d17,d18,layout_matrix = rbind(c(1, 2, 3),c(4, 5, 6)), bottom=bottomt,left=leftt)
print(plot1)
dev.off()



# 20201117 OLDFIRES & ASPECT ##### Extract aspect and add to incsv--------
# use points themselves & NOT buffered points
library(rgdal); library(maptools)
setwd("N:\\project\\monitoring_vol\\SoCalShrubs\\points")
aspect <- raster("N://rasterlib//Climate//UNR//modern//elev//ne_sw_aspect.dat")
# u11
u11points <- readOGR(dsn="N:\\project\\monitoring_vol\\SoCalShrubs\\points",layer="u11_allpts_j")
u11df <- as.data.frame(u11points)
#u11samp <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\u11\\ARGmax0903_0_ARG_max_u11.tif")
#asp11 <- projectRaster(aspect, u11samp, method="bilinear", alignOnly=FALSE, over=FALSE)
#writeRaster(asp11, "aspectu11.tif",overwrite=T) # for QC in arcmap 

asp11 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\points\\aspectu11.tif"); names(asp11) <- 'aspect'

# output from extract will not have correct ID numbers so need to join it to u11df with correct FIDp1 values
vmax <- extract(asp11,u11points, nl = 1, df=TRUE,na.rm=TRUE, buffer=15,fun=mean)
u11aspect <- cbind(u11df,vmax)
u11aspect$pltnum <- as.integer(u11aspect$FIDp1) + 11000

# u10
u10points <- readOGR(dsn="N:\\project\\monitoring_vol\\SoCalShrubs\\points",layer="u10_allpts_j")
u10df <- as.data.frame(u10points)

#u10samp <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\u10\\ARGmax0903_0_ARG_max_u10.tif")
#asp10 <- projectRaster(aspect, u10samp, method="bilinear", alignOnly=FALSE, over=FALSE)
#writeRaster(asp10, "aspectu10.tif",overwrite=T) # for QC in arcmap 

asp10 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\points\\aspectu10.tif"); names(asp10) <- 'aspect'


# output from extract will not have correct ID numbers so need to join it to u10df with correct FIDp1 values
vmax <- extract(asp10,u10points, nl = 1, df=TRUE,na.rm=TRUE, buffer=15,fun=mean)
u10aspect <- cbind(u10df,vmax)
u10aspect$pltnum <- as.integer(u10aspect$FIDp1) + 10000

u1011aspectdf <- rbind(u11aspect,u10aspect)
formerge <- u1011aspectdf[,267:268]
incsv <- read.table("N:/project/monitoring_vol/SoCalShrubs/points/u1011_2009_2018_mar2novby3_20201104.csv", header = TRUE, sep=",")
csvasp <- merge(incsv,formerge, by.x='pltnum',by.y='pltnum',sort=TRUE)
write.csv(csvasp,"N:/project/monitoring_vol/SoCalShrubs/points/include_aspect_20201124.csv")

# 20201117 OLDFIRES & ASPECT ##### Continue prepping data after aspect has been added ---------
incsv <- csvasp
y09 <- incsv[c(2:77,453,459,463,469,1)] # aspect = col 469
y11 <- incsv[c(2,78:152,453,458,464,469,1)]
y13 <- incsv[c(2,153:227,453,457,465,469,1)]
y15 <- incsv[c(2,228:302,453,456,466,469,1)]
y17 <- incsv[c(2,303:377,453,455,467,469,1)]
y18 <- incsv[c(2,378:454,468,469,1)]
y09 <- y09[complete.cases(y09),]
y11 <- y11[complete.cases(y11),]
y13 <- y13[complete.cases(y13),]
y15 <- y15[complete.cases(y15),]
y17 <- y17[complete.cases(y17),]
y18 <- y18[complete.cases(y18),]



# 20201117 OLDFIRES & ASPECT ##### Generate residuals for model yr 2009 for all years & save tables in cvtables_aspect_20201117 #####-------
# These input data include all records, regardless of whether there was a fire or not, and aspect
require(broom)
# OLDFIRES & ASPECT - Generate residuals for model yr 2009 with 2009 data
rfts <- randomForest(heb_2009 ~ ARGmax0903 + NBRstd0903 + ARGmed0906 + ts2009 + aspect, 
                     data = y09, importance = TRUE, type = "regression", na.action = na.exclude)
importance(rfts)

y09$pred_m09d09 <- predict(rfts)
y09$residm09d09 <- y09$pred_m09d09 - y09$heb_2009
# Separate by presence or absence of past fires. 'timesince' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2009 # data year
cvtabof[1,2] <- 2009 # model year
cvtabof[1,3] <- mean(y09$residm09d09[y09$ts2009 > 199])
cvtabof[1,4] <- mean(y09$residm09d09[y09$ts2009 < 200])
d <- lm(y09$residm09d09 ~ y09$ts2009) # lm for ALL recs
# add ts2009 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y09$residm09d09[y09$ts2009 < 200] ~ y09$ts2009[y09$ts2009 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y09$residm09d09)
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_aspect_20201117\\rfcv_m09d09_mresid_all.csv"
write.csv(cvtabof,tname)

# OLDFIRES & ASPECT - Generate residuals for model year 2009 with 2011 data
rfts <- randomForest(heb_2011 ~ ARGmax1103 + NBRstd1103 + ARGmed1106 + ts2011 + aspect, 
                     data = y11, importance = TRUE, type = "regression", na.action = na.exclude)
# Generate residuals
y11$pred_m09d11 <- predict(rfts)
y11$residm09d11 <- y11$pred_m09d11 - y11$heb_2011
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2011 # data year
cvtabof[1,2] <- 2009 # model year
cvtabof[1,3] <- mean(y11$residm09d11[y11$ts2011 > 199]) # no past fires
cvtabof[1,4] <- mean(y11$residm09d11[y11$ts2011 < 200]) # yes past fires
d <- lm(y11$residm09d11 ~ y11$ts2011)
# add ts2011 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y11$residm09d11[y11$ts2011 < 200] ~ y11$ts2011[y11$ts2011 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y11$residm09d11)
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_aspect_20201117\\rfcv_m09d11_mresid_all.csv"
write.csv(cvtabof,tname)

# OLDFIRES & ASPECT - Generate residuals for model year 2009 with 2013 data
rfts <- randomForest(heb_2013 ~ ARGmax1303 + NBRstd1303 + ARGmed1306 + ts2013 + ne_sw_aspect, 
                     data = y13, importance = TRUE, type = "regression", na.action = na.exclude)
y13$pred_m09d13 <- predict(rfts)
y13$residm09d13 <- y13$pred_m09d13 - y13$heb_2013
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2013 # data year
cvtabof[1,2] <- 2009 # model year
cvtabof[1,3] <- mean(y13$residm09d13[y13$ts2013 > 199])
cvtabof[1,4] <- mean(y13$residm09d13[y13$ts2013 < 200])
d <- lm(y13$residm09d13 ~ y13$ts2013)
# add ts2013 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y13$residm09d13[y13$ts2013 < 200] ~ y13$ts2013[y13$ts2013 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y13$residm09d13)
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_aspect_20201117\\rfcv_m09d13_mresid_all.csv"
write.csv(cvtabof,tname)

# OLDFIRES & ASPECT - Generate residuals for model year 2009 with 2015 data
rfts <- randomForest(heb_2015 ~ ARGmax1503 + NBRstd1503 + ARGmed1506 + ts2015 + ne_sw_aspect, 
                     data = y15, importance = TRUE, type = "regression", na.action = na.exclude)
y15$pred_m09d15 <- predict(rfts)
y15$residm09d15 <- y15$pred_m09d15 - y15$heb_2015
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2015 # data year
cvtabof[1,2] <- 2009 # model year
cvtabof[1,3] <- mean(y15$residm09d15[y15$ts2015 > 199])
cvtabof[1,4] <- mean(y15$residm09d15[y15$ts2015 < 200])
d <- lm(y15$residm09d15 ~ y15$ts2015)
# add ts2015 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y15$residm09d15[y15$ts2015 < 200] ~ y15$ts2015[y15$ts2015 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y15$residm09d15)
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_aspect_20201117\\rfcv_m09d15_mresid_all.csv"
write.csv(cvtabof,tname)

#  OLDFIRES & ASPECT - Generate residuals for model year 2009 with 2017 data
rfts <- randomForest(heb_2017 ~ ARGmax1703 + NBRstd1703 + ARGmed1706 + ts2017 + ne_sw_aspect, 
                     data = y17, importance = TRUE, type = "regression", na.action = na.exclude)
y17$pred_m09d17 <- predict(rfts)
y17$residm09d17 <- y17$pred_m09d17 - y17$heb_2017
# Separate by presence or absence of past fires. 'ts2017' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2017 # data year
cvtabof[1,2] <- 2009 # model year
cvtabof[1,3] <- mean(y17$residm09d17[y17$ts2017 > 199]) # Select where we have no fire history
cvtabof[1,4] <- mean(y17$residm09d17[y17$ts2017 < 200]) # select where there are previous fires
d <- lm(y17$residm09d17 ~ y17$ts2017)
# add ts2017 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y17$residm09d17[y17$ts2017 < 200] ~ y17$ts2017[y17$ts2017 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y17$residm09d17)
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_aspect_20201117\\rfcv_m09d17_mresid_all.csv"
write.csv(cvtabof,tname)

# OLDFIRES & ASPECT - Generate residuals for model year 2009 with 2018 data
rfts <- randomForest(Herb_2018 ~ ARGmax1803 + NBRstd1803 + ARGmed1806 + ts2018 + ne_sw_aspect, 
                     data = y18, importance = TRUE, type = "regression", na.action = na.exclude)
y18$pred_m09d18 <- predict(rfts)
y18$residm09d18 <- y18$pred_m09d18 - y18$Herb_2018
# Separate by presence or absence of past fires. 'ts2011' for NO past fires = 200
cvtabof <- data.frame()
cvtabof[1,1] <- 2018 # data year
cvtabof[1,2] <- 2009 # model year
cvtabof[1,3] <- mean(y18$residm09d18[y18$ts2018 > 199])
cvtabof[1,4] <- mean(y18$residm09d18[y18$ts2018 < 200])
d <- lm(y18$residm09d18 ~ y18$ts2018)
# add ts2018 Estimate & p-value from lm summary to the table
cvtabof[1,5] <- d$coefficients[2]
cvtabof[1,6] <- glance(d)$p.value
dYof <- lm(y18$residm09d18[y18$ts2018 < 200] ~ y18$ts2018[y18$ts2018 < 200]) #lm for ONLY past fires
cvtabof[1,7] <- dYof$coefficients[2] #YES past fires
cvtabof[1,8] <- glance(dYof)$p.value # YES past fires
cvtabof[1,9] <- mean(y18$residm09d18)
colnames(cvtabof) <- c("datayr","modelyr","mresid_NoPF","mresid_YPF","Alllm_estimate","Allpvalue","YPFlm_estimate","YPFpvalue","mresid_All")
tname <- "N:\\project\\monitoring_vol\\SoCalShrubs\\EE_Extractions\\cvtables_aspect_20201117\\rfcv_m09d18_mresid_all.csv"
write.csv(cvtabof,tname)






#################################################################################################

# -------- Subtract predicted rasters u11 ----------
setwd("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\predictions20200616")
rf2009_2013bm_u11 <- raster("rf2009_2013bm_u11.tif")
rf2011_2011bm_u11 <- raster("rf2011_2011bm_u11.tif")
rf2013_2013bm_u11 <- raster("rf2013_2013bm_u11.tif")
rf2015_2015bm_u11 <- raster("rf2015_2015bm_u11.tif")
rf2017_2017bm_u11 <- raster("rf2017_2017bm_u11.tif")
rf2018_2017bm_u11 <- raster("rf2018_2017bm_u11.tif")
bestyrs_u11 <- stack(rf2009_2013bm_u11,rf2011_2011bm_u11,rf2013_2013bm_u11,rf2015_2015bm_u11,rf2017_2017bm_u11,rf2018_2017bm_u11)

dy09m11 <- bestyrs_u11$rf2009_2013bm_u11 - bestyrs_u11$rf2011_2011bm_u11
dy09m13 <- bestyrs_u11$rf2009_2013bm_u11 - bestyrs_u11$rf2013_2013bm_u11
dy09m15 <- bestyrs_u11$rf2009_2013bm_u11 - bestyrs_u11$rf2015_2015bm_u11
dy09m17 <- bestyrs_u11$rf2009_2013bm_u11 - bestyrs_u11$rf2017_2017bm_u11
dy09m18 <- bestyrs_u11$rf2009_2013bm_u11 - bestyrs_u11$rf2018_2017bm_u11

dy11m09 <- bestyrs_u11$rf2011_2011bm_u11 - bestyrs_u11$rf2009_2013bm_u11
dy11m13 <- bestyrs_u11$rf2011_2011bm_u11 - bestyrs_u11$rf2013_2013bm_u11
dy11m15 <- bestyrs_u11$rf2011_2011bm_u11 - bestyrs_u11$rf2015_2015bm_u11
dy11m17 <- bestyrs_u11$rf2011_2011bm_u11 - bestyrs_u11$rf2017_2017bm_u11
dy11m18 <- bestyrs_u11$rf2011_2011bm_u11 - bestyrs_u11$rf2018_2017bm_u11

dy13m09 <- bestyrs_u11$rf2013_2013bm_u11 - bestyrs_u11$rf2009_2013bm_u11
dy13m11 <- bestyrs_u11$rf2013_2013bm_u11 - bestyrs_u11$rf2011_2011bm_u11
dy13m15 <- bestyrs_u11$rf2013_2013bm_u11 - bestyrs_u11$rf2015_2015bm_u11
dy13m17 <- bestyrs_u11$rf2013_2013bm_u11 - bestyrs_u11$rf2017_2017bm_u11
dy13m18 <- bestyrs_u11$rf2013_2013bm_u11 - bestyrs_u11$rf2018_2017bm_u11

dy15m09 <- bestyrs_u11$rf2015_2015bm_u11 - bestyrs_u11$rf2009_2013bm_u11
dy15m11 <- bestyrs_u11$rf2015_2015bm_u11 - bestyrs_u11$rf2011_2011bm_u11
dy15m13 <- bestyrs_u11$rf2015_2015bm_u11 - bestyrs_u11$rf2013_2013bm_u11
dy15m17 <- bestyrs_u11$rf2015_2015bm_u11 - bestyrs_u11$rf2017_2017bm_u11
dy15m18 <- bestyrs_u11$rf2015_2015bm_u11 - bestyrs_u11$rf2018_2017bm_u11

dy17m09 <- bestyrs_u11$rf2017_2017bm_u11 - bestyrs_u11$rf2009_2013bm_u11
dy17m11 <- bestyrs_u11$rf2017_2017bm_u11 - bestyrs_u11$rf2011_2011bm_u11
dy17m13 <- bestyrs_u11$rf2017_2017bm_u11 - bestyrs_u11$rf2013_2013bm_u11
dy17m15 <- bestyrs_u11$rf2017_2017bm_u11 - bestyrs_u11$rf2015_2015bm_u11
dy17m18 <- bestyrs_u11$rf2017_2017bm_u11 - bestyrs_u11$rf2018_2017bm_u11

dy18m09 <- bestyrs_u11$rf2018_2017bm_u11 - bestyrs_u11$rf2009_2013bm_u11
dy18m11 <- bestyrs_u11$rf2018_2017bm_u11 - bestyrs_u11$rf2011_2011bm_u11
dy18m13 <- bestyrs_u11$rf2018_2017bm_u11 - bestyrs_u11$rf2013_2013bm_u11
dy18m15 <- bestyrs_u11$rf2018_2017bm_u11 - bestyrs_u11$rf2015_2015bm_u11
dy18m17 <- bestyrs_u11$rf2018_2017bm_u11 - bestyrs_u11$rf2017_2017bm_u11


# -------- Compare predicted rasters u10 ----------
setwd("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\predictions20200616")
rf2009_2013bm_u10 <- raster("rf2009_2013bm_u10.tif")
rf2011_2011bm_u10 <- raster("rf2011_2011bm_u10.tif")
rf2013_2013bm_u10 <- raster("rf2013_2013bm_u10.tif")
rf2015_2015bm_u10 <- raster("rf2015_2015bm_u10.tif")
rf2017_2017bm_u10 <- raster("rf2017_2017bm_u10.tif")
rf2018_2017bm_u10 <- raster("rf2018_2017bm_u10.tif")
bestyrs_u10 <- stack(rf2009_2013bm_u10,rf2011_2011bm_u10,rf2013_2013bm_u10,rf2015_2015bm_u10,rf2017_2017bm_u10,rf2018_2017bm_u10)

dy09m11 <- bestyrs_u10$rf2009_2013bm_u10 - bestyrs_u10$rf2011_2011bm_u10
dy09m13 <- bestyrs_u10$rf2009_2013bm_u10 - bestyrs_u10$rf2013_2013bm_u10
dy09m15 <- bestyrs_u10$rf2009_2013bm_u10 - bestyrs_u10$rf2015_2015bm_u10
dy09m17 <- bestyrs_u10$rf2009_2013bm_u10 - bestyrs_u10$rf2017_2017bm_u10
dy09m18 <- bestyrs_u10$rf2009_2013bm_u10 - bestyrs_u10$rf2018_2017bm_u10

dy11m09 <- bestyrs_u10$rf2011_2011bm_u10 - bestyrs_u10$rf2009_2013bm_u10
dy11m13 <- bestyrs_u10$rf2011_2011bm_u10 - bestyrs_u10$rf2013_2013bm_u10
dy11m15 <- bestyrs_u10$rf2011_2011bm_u10 - bestyrs_u10$rf2015_2015bm_u10
dy11m17 <- bestyrs_u10$rf2011_2011bm_u10 - bestyrs_u10$rf2017_2017bm_u10
dy11m18 <- bestyrs_u10$rf2011_2011bm_u10 - bestyrs_u10$rf2018_2017bm_u10

dy13m09 <- bestyrs_u10$rf2013_2013bm_u10 - bestyrs_u10$rf2009_2013bm_u10
dy13m11 <- bestyrs_u10$rf2013_2013bm_u10 - bestyrs_u10$rf2011_2011bm_u10
dy13m15 <- bestyrs_u10$rf2013_2013bm_u10 - bestyrs_u10$rf2015_2015bm_u10
dy13m17 <- bestyrs_u10$rf2013_2013bm_u10 - bestyrs_u10$rf2017_2017bm_u10
dy13m18 <- bestyrs_u10$rf2013_2013bm_u10 - bestyrs_u10$rf2018_2017bm_u10

dy15m09 <- bestyrs_u10$rf2015_2015bm_u10 - bestyrs_u10$rf2009_2013bm_u10
dy15m11 <- bestyrs_u10$rf2015_2015bm_u10 - bestyrs_u10$rf2011_2011bm_u10
dy15m13 <- bestyrs_u10$rf2015_2015bm_u10 - bestyrs_u10$rf2013_2013bm_u10
dy15m17 <- bestyrs_u10$rf2015_2015bm_u10 - bestyrs_u10$rf2017_2017bm_u10
dy15m18 <- bestyrs_u10$rf2015_2015bm_u10 - bestyrs_u10$rf2018_2017bm_u10

dy17m09 <- bestyrs_u10$rf2017_2017bm_u10 - bestyrs_u10$rf2009_2013bm_u10
dy17m11 <- bestyrs_u10$rf2017_2017bm_u10 - bestyrs_u10$rf2011_2011bm_u10
dy17m13 <- bestyrs_u10$rf2017_2017bm_u10 - bestyrs_u10$rf2013_2013bm_u10
dy17m15 <- bestyrs_u10$rf2017_2017bm_u10 - bestyrs_u10$rf2015_2015bm_u10
dy17m18 <- bestyrs_u10$rf2017_2017bm_u10 - bestyrs_u10$rf2018_2017bm_u10

dy18m09 <- bestyrs_u10$rf2018_2017bm_u10 - bestyrs_u10$rf2009_2013bm_u10
dy18m11 <- bestyrs_u10$rf2018_2017bm_u10 - bestyrs_u10$rf2011_2011bm_u10
dy18m13 <- bestyrs_u10$rf2018_2017bm_u10 - bestyrs_u10$rf2013_2013bm_u10
dy18m15 <- bestyrs_u10$rf2018_2017bm_u10 - bestyrs_u10$rf2015_2015bm_u10
dy18m17 <- bestyrs_u10$rf2018_2017bm_u10 - bestyrs_u10$rf2017_2017bm_u10
dy18m17_df <- as.data.frame(dy18m17)
histogram(dy18m17_df)



# -------- Compare extra rasters u10 ---------added 20200714-----
setwd("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\predictions20200616")
# data year 2017, with 2017 base minus 2013 base model
rf2017_2017bm_u10 <- raster("rf2017_2017bm_u10.tif")
rf2017_2013bm_u10 <- raster("rf2017_2013bm_u10.tif")

# data year 2013, with 2017 base minus 2013 base model 
rf2013_2017bm_u10 <- raster("rf2013_2017bm_u10.tif")
rf2013_2013bm_u10 <- raster("rf2013_2013bm_u10.tif")

# data year 2013, with 2013 base minus 2009 base model
rf2013_2009bm_u10 <- raster("rf2013_2009bm_u10.tif")
#rf2013_2013bm_u10 <- raster("rf2013_2013bm_u10.tif")

# data year 2009, with 2013 base minus 2009 base model
rf2009_2009bm_u10 <- raster("rf2009_2009bm_u10.tif")
rf2009_2013bm_u10 <- raster("rf2009_2013bm_u10.tif")

# data year 2017, with 2017 base minus 2009 base model
rf2017_2009bm_u10 <- raster("rf2017_2009bm_u10.tif")
#rf2017_2017bm_u10 <- raster("rf2017_2017bm_u10.tif")

# data year 2009, with 2017 base minus 2009 base model
#rf2009_2009bm_u10 <- raster("rf2009_2009bm_u10.tif")
rf2009_2017bm_u10 <- raster("rf2009_2017bm_u10.tif")

bestyrs_u10 <- stack(rf2017_2017bm_u10,rf2017_2013bm_u10,rf2013_2017bm_u10,rf2013_2013bm_u10,rf2013_2009bm_u10,rf2009_2009bm_u10,rf2009_2013bm_u10,rf2017_2009bm_u10,rf2009_2017bm_u10)

# data year 2017, with 2017 base minus 2013 base model
dy17_bm17mbm13_u10 <- bestyrs_u10$rf2017_2017bm_u10 - bestyrs_u10$rf2017_2013bm_u10
writeRaster(dy17_bm17mbm13_u10,"dy17_bm17mbm13_u10.tif",overwrite=TRUE)

# data year 2013, with 2017 base minus 2013 base model
dy13_bm17mbm13_u10 <- bestyrs_u10$rf2013_2017bm_u10 - bestyrs_u10$rf2013_2013bm_u10
writeRaster(dy13_bm17mbm13_u10,"dy13_bm17mbm13_u10.tif",overwrite=TRUE)

# data year 2013, with 2013 base minus 2009 base model
dy13_bm13mbm09_u10 <- bestyrs_u10$rf2013_2013bm_u10 - bestyrs_u10$rf2013_2009bm_u10
writeRaster(dy13_bm13mbm09_u10,"dy13_bm13mbm09_u10.tif",overwrite=TRUE)

# data year 2009, with 2013 base minus 2009 base model
dy09_bm13mbm09_u10 <- bestyrs_u10$rf2009_2013bm_u10 - bestyrs_u10$rf2009_2009bm_u10
writeRaster(dy09_bm13mbm09_u10,"dy09_bm13mbm09_u10.tif",overwrite=TRUE)

# data year 2017, with 2017 base minus 2009 base model
dy17_bm17mbm09_u10 <- bestyrs_u10$rf2017_2017bm_u10 - bestyrs_u10$rf2017_2009bm_u10
writeRaster(dy17_bm17mbm09_u10,"dy17_bm17mbm09_u10.tif",overwrite=TRUE)

# data year 2009, with 2017 base minus 2009 base model
dy09_bm17mbm09_u10 <- bestyrs_u10$rf2009_2017bm_u10 - bestyrs_u10$rf2009_2009bm_u10
writeRaster(dy09_bm17mbm09_u10,"dy09_bm17mbm09_u10.tif",overwrite=TRUE)


# -------- Compare extra rasters u11 ---------added 20200715-----
setwd("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\predictions20200616")
# data year 2017, with 2017 base minus 2013 base model
rf2017_2017bm_u11 <- raster("rf2017_2017bm_u11.tif")
rf2017_2013bm_u11 <- raster("rf2017_2013bm_u11.tif")

# data year 2013, with 2017 base minus 2013 base model 
rf2013_2017bm_u11 <- raster("rf2013_2017bm_u11.tif")
rf2013_2013bm_u11 <- raster("rf2013_2013bm_u11.tif")

# data year 2013, with 2013 base minus 2009 base model
rf2013_2009bm_u11 <- raster("rf2013_2009bm_u11.tif")
#rf2013_2013bm_u11 <- raster("rf2013_2013bm_u11.tif")

# data year 2009, with 2013 base minus 2009 base model
rf2009_2009bm_u11 <- raster("rf2009_2009bm_u11.tif")
rf2009_2013bm_u11 <- raster("rf2009_2013bm_u11.tif")

# data year 2017, with 2017 base minus 2009 base model
rf2017_2009bm_u11 <- raster("rf2017_2009bm_u11.tif")
#rf2017_2017bm_u11 <- raster("rf2017_2017bm_u11.tif")

# data year 2009, with 2017 base minus 2009 base model
#rf2009_2009bm_u11 <- raster("rf2009_2009bm_u11.tif")
rf2009_2017bm_u11 <- raster("rf2009_2017bm_u11.tif")

bestyrs_u11 <- stack(rf2017_2017bm_u11,rf2017_2013bm_u11,rf2013_2017bm_u11,rf2013_2013bm_u11,rf2013_2009bm_u11,rf2009_2009bm_u11,rf2009_2013bm_u11,rf2017_2009bm_u11,rf2009_2017bm_u11)

# data year 2017, with 2017 base minus 2013 base model
dy17_bm17mbm13_u11 <- bestyrs_u11$rf2017_2017bm_u11 - bestyrs_u11$rf2017_2013bm_u11
writeRaster(dy17_bm17mbm13_u11,"dy17_bm17mbm13_u11.tif",overwrite=TRUE)

# data year 2013, with 2017 base minus 2013 base model
dy13_bm17mbm13_u11 <- bestyrs_u11$rf2013_2017bm_u11 - bestyrs_u11$rf2013_2013bm_u11
writeRaster(dy13_bm17mbm13_u11,"dy13_bm17mbm13_u11.tif",overwrite=TRUE)

# data year 2013, with 2013 base minus 2009 base model
dy13_bm13mbm09_u11 <- bestyrs_u11$rf2013_2013bm_u11 - bestyrs_u11$rf2013_2009bm_u11
writeRaster(dy13_bm13mbm09_u11,"dy13_bm13mbm09_u11.tif",overwrite=TRUE)

# data year 2009, with 2013 base minus 2009 base model
dy09_bm13mbm09_u11 <- bestyrs_u11$rf2009_2013bm_u11 - bestyrs_u11$rf2009_2009bm_u11
writeRaster(dy09_bm13mbm09_u11,"dy09_bm13mbm09_u11.tif",overwrite=TRUE)

# data year 2017, with 2017 base minus 2009 base model
dy17_bm17mbm09_u11 <- bestyrs_u11$rf2017_2017bm_u11 - bestyrs_u11$rf2017_2009bm_u11
writeRaster(dy17_bm17mbm09_u11,"dy17_bm17mbm09_u11.tif",overwrite=TRUE)

# data year 2009, with 2017 base minus 2009 base model
dy09_bm17mbm09_u11 <- bestyrs_u11$rf2009_2017bm_u11 - bestyrs_u11$rf2009_2009bm_u11
writeRaster(dy09_bm17mbm09_u11,"dy09_bm17mbm09_u11.tif",overwrite=TRUE)

# -------------- Histograms from subtracted models - 20200925 ---------------------
# Need to combine then convert to data frame

# -but use 5 percentile bins
# -Divide by fire vs non-fire using an 'oldfire' field
# Need 12 histograms in total
#-histograms of differences
library(raster)
library(ggplot2)

#  Separate difference rasters by 'old fires' & 'not old fires' & create histograms ----------

# Open oldfire rasters & fix extents
oldfires_u10 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\SeparateByOldfire20200928\\oldfires_u10.tif")
setwd("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\predictions20200616")
dy17_bm17mbm13_u10 <- raster("dy17_bm17mbm13_u10.tif")
oldfires_u10_clip <- crop(oldfires_u10,dy17_bm17mbm13_u10)
oldfires_u10e <- raster(vals=values(oldfires_u10_clip),ext=extent(dy17_bm17mbm13_u10),crs=crs(dy17_bm17mbm13_u10),
              nrows=dim(dy17_bm17mbm13_u10)[1],ncols=dim(dy17_bm17mbm13_u10)[2])


oldfires_u11 <- raster("N:\\project\\monitoring_vol\\SoCalShrubs\\EE_rasters\\SeparateByOldfire20200928\\oldfires_u11.tif")
#dy17_bm17mbm13_u11 <- raster("dy17_bm17mbm13_u11.tif")
#oldfires_u11_clip <- crop(oldfires_u11,dy17_bm17mbm13_u11)
#oldfires_u11e <- raster(vals=values(oldfires_u11_clip),ext=extent(dy17_bm17mbm13_u11),crs=crs(dy17_bm17mbm13_u11),
#                        nrows=dim(dy17_bm17mbm13_u11)[1],ncols=dim(dy17_bm17mbm13_u11)[2])

# separate by 0 & 1
oldfires_u11_1s <- (oldfires_u11 ==1); oldfires_u11_1s[oldfires_u11_1s == 0] <- NA
oldfires_u11_0s <- oldfires_u11e == 0; oldfires_u11_0s[oldfires_u11_0s == 0] <- NA; oldfires_u11_0s[oldfires_u11_0s == 1] <- 0

oldfires_u10_1s <- (oldfires_u10 ==1); oldfires_u10_1s[oldfires_u10_1s == 0] <- NA
oldfires_u10_0s <- oldfires_u10e == 0; oldfires_u10_0s[oldfires_u10_0s == 0] <- NA; oldfires_u10_0s[oldfires_u10_0s == 1] <- 0

st_u11 <- stack("dy17_bm17mbm13_u11.tif","dy13_bm17mbm13_u11.tif","dy13_bm13mbm09_u11.tif","dy09_bm13mbm09_u11.tif","dy17_bm17mbm09_u11.tif","dy09_bm17mbm09_u11.tif")
names(st_u11) <- c("dy17_bm17mbm13","dy13_bm17mbm13","dy13_bm13mbm09","dy09_bm13mbm09","dy17_bm17mbm09","dy09_bm17mbm09")

st_u11_1s <- mask(st_u11,oldfires_u11_1s)
st_u11_0s <- mask(st_u11,oldfires_u11_0s)

st_u10 <- stack("dy17_bm17mbm13_u10.tif","dy13_bm17mbm13_u10.tif","dy13_bm13mbm09_u10.tif","dy09_bm13mbm09_u10.tif","dy17_bm17mbm09_u10.tif","dy09_bm17mbm09_u10.tif")
names(st_u10) <- c("dy17_bm17mbm13","dy13_bm17mbm13","dy13_bm13mbm09","dy09_bm13mbm09","dy17_bm17mbm09","dy09_bm17mbm09")

st_u10_1s <- mask(st_u10,oldfires_u10_1s)
st_u10_0s <- mask(st_u10,oldfires_u10_0s)





st_u11_1s <- mask(st_u10)

# Generate matrics - too large in some cases
# stu10_oldfires <- st_u10[oldfires_u10 == 1,]  # old fire areas
# stu10_nooldfires <- st_u10[oldfires_u10 == 0,]  # no old fire area
# stu11_oldfires <- st_u11[oldfires_u11 == 1,]  # old fire areas
# stu11_nooldfires <- st_u11[oldfires_u11 == 0,]  # no old fire areas

u10_nf <- stackSelect(st_u10,oldfires_u10e ==1)


dy17_bm17mbm13_df <- rbind(as.data.frame(st_u10$dy17_bm17mbm13, na.rm=TRUE),as.data.frame(st_u11$dy17_bm17mbm13, na.rm=TRUE))
dy13_bm17mbm13_df <- rbind(as.data.frame(st_u10$dy13_bm17mbm13, na.rm=TRUE),as.data.frame(st_u11$dy13_bm17mbm13, na.rm=TRUE))
dy13_bm13mbm09_df <- rbind(as.data.frame(st_u10$dy13_bm13mbm09, na.rm=TRUE),as.data.frame(st_u11$dy13_bm13mbm09, na.rm=TRUE))
dy09_bm13mbm09_df <- rbind(as.data.frame(st_u10$dy09_bm13mbm09, na.rm=TRUE),as.data.frame(st_u11$dy09_bm13mbm09, na.rm=TRUE))
dy17_bm17mbm09_df <- rbind(as.data.frame(st_u10$dy17_bm17mbm09, na.rm=TRUE),as.data.frame(st_u11$dy17_bm17mbm09, na.rm=TRUE))
dy09_bm17mbm09_df <- rbind(as.data.frame(st_u10$dy09_bm17mbm09, na.rm=TRUE),as.data.frame(st_u11$dy09_bm17mbm09, na.rm=TRUE))
head(dy09_bm17mbm09_df)
names(st_u10) <- c("dy17_bm17mbm13","dy13_bm17mbm13","dy13_bm13mbm09","dy09_bm13mbm09","dy17_bm17mbm09","dy09_bm17mbm09")
rs <- c("dy17_bm17mbm13","dy13_bm17mbm13","dy13_bm13mbm09","dy09_bm13mbm09","dy17_bm17mbm09","dy09_bm17mbm09")

# st_u10df <- as.data.frame(st_u10,na.rm=TRUE)
# st_u11df <- as.data.frame(st_u11,na.rm=TRUE)
# st_df <- rbind(st_u11df,st_u10df)

# ****** Access column name with variable. Use df[,"nominalprice"] to obtain that column******


##  dy17_bm17mbm13_df histogram -----
dy17_bm17mbm13_df_s500000 <- as.data.frame(dy17_bm17mbm13_df[sample(nrow(dy17_bm17mbm13_df), 500000), ]) # take random sample to reduce size
names(dy17_bm17mbm13_df_s500000) <- c("dy17_bm17mbm13")
plot.new()
pdf("NoFires_dy17_bm17mbm13_500000_randomsample.pdf")
ggplot(dy17_bm17mbm13_df_s500000, aes(x=dy17_bm17mbm13)) + 
  geom_histogram(binwidth=2,color="black", fill="light blue") + 
  xlab("dy17_bm17mbm13") + 
  theme(axis.title.x = element_text(face="plain",color="black",size=16)) +
  ylab("Count") +
  xlim(c(-50,50)) +
  theme(axis.title.y = element_text(face="plain",color="black",size=16))  +
  #annotate("text", 50, -50, label="No Fires dy17_bm17mbm13 50000 random sample", color = "black",hjust = 1,size = 6)
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("No Fires dy17_bm17mbm13 500,000 random sample")
dev.off()

##  dy13_bm17mbm13_df histogram -----
dy13_bm17mbm13_df_s500000 <- as.data.frame(dy13_bm17mbm13_df[sample(nrow(dy13_bm17mbm13_df), 500000), ])
names(dy13_bm17mbm13_df_s500000) <- c("dy13_bm17mbm13")
plot.new()
pdf("NoFires_dy13_bm17mbm13_500000_randomsample.pdf")
ggplot(dy13_bm17mbm13_df_s500000, aes(x=dy13_bm17mbm13)) + 
  geom_histogram(binwidth=2,color="black", fill="light blue") + 
  xlab("dy13_bm17mbm13") + 
  theme(axis.title.x = element_text(face="plain",color="black",size=16)) +
  ylab("Count") +
  xlim(c(-50,50)) +
  theme(axis.title.y = element_text(face="plain",color="black",size=16))  +
  #annotate("text", 50, -50, label="No Fires dy17_bm17mbm13 50000 random sample", color = "black",hjust = 1,size = 6)
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("No Fires dy13_bm17mbm13 500,000 random sample")
dev.off()

##  dy13_bm13mbm09_df histogram -----
dy13_bm13mbm09_df_s500000 <- as.data.frame(dy13_bm13mbm09_df[sample(nrow(dy13_bm13mbm09_df), 500000), ])
names(dy13_bm13mbm09_df_s500000) <- c("dy13_bm13mbm09")
plot.new()
pdf("NoFires_dy13_bm13mbm09_500000_randomsample.pdf")
ggplot(dy13_bm13mbm09_df_s500000, aes(x=dy13_bm13mbm09)) + 
  geom_histogram(binwidth=2,color="black", fill="light blue") + 
  xlab("dy13_bm13mbm09") + 
  theme(axis.title.x = element_text(face="plain",color="black",size=16)) +
  ylab("Count") +
  xlim(c(-45,45)) +
  theme(axis.title.y = element_text(face="plain",color="black",size=16))  +
  #annotate("text", 50, -50, label="No Fires dy17_bm17mbm13 50000 random sample", color = "black",hjust = 1,size = 6)
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("No Fires dy13_bm13mbm09 500,000 random sample")
dev.off()

##  dy09_bm13mbm09_df histogram -----
dy09_bm13mbm09_df_s500000 <- as.data.frame(dy09_bm13mbm09_df[sample(nrow(dy09_bm13mbm09_df), 500000), ])
names(dy09_bm13mbm09_df_s500000) <- c("dy09_bm13mbm09")
plot.new()
pdf("NoFires_dy09_bm13mbm09_500000_randomsample.pdf")
ggplot(dy09_bm13mbm09_df_s500000, aes(x=dy09_bm13mbm09)) + 
  geom_histogram(binwidth=2,color="black", fill="light blue") + 
  xlab("dy13_bm13mbm09") + 
  theme(axis.title.x = element_text(face="plain",color="black",size=16)) +
  ylab("Count") +
  xlim(c(-45,45)) +
  theme(axis.title.y = element_text(face="plain",color="black",size=16))  +
  #annotate("text", 50, -50, label="No Fires dy17_bm17mbm13 50000 random sample", color = "black",hjust = 1,size = 6)
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("No Fires dy09_bm13mbm09 500,000 random sample")
dev.off()

##  dy17_bm17mbm09_df histogram -----
dy17_bm17mbm09_df_s500000 <- as.data.frame(dy17_bm17mbm09_df[sample(nrow(dy17_bm17mbm09_df), 500000), ])
names(dy17_bm17mbm09_df_s500000) <- c("dy17_bm17mbm09")
plot.new()
pdf("NoFires_dy17_bm17mbm09_500000_randomsample.pdf")
ggplot(dy17_bm17mbm09_df_s500000, aes(x=dy17_bm17mbm09)) + 
  geom_histogram(binwidth=2,color="black", fill="light blue") + 
  xlab("dy13_bm13mbm09") + 
  theme(axis.title.x = element_text(face="plain",color="black",size=16)) +
  ylab("Count") +
  xlim(c(-45,45)) +
  theme(axis.title.y = element_text(face="plain",color="black",size=16))  +
  #annotate("text", 50, -50, label="No Fires dy17_bm17mbm13 50000 random sample", color = "black",hjust = 1,size = 6)
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("No Fires dy17_bm17mbm09 500,000 random sample")
dev.off()

##  dy09_bm17mbm09_df histogram -----
dy09_bm17mbm09_df_s500000 <- as.data.frame(dy09_bm17mbm09_df[sample(nrow(dy09_bm17mbm09_df), 500000), ])
names(dy09_bm17mbm09_df_s500000) <- c("dy09_bm17mbm09")
plot.new()
pdf("NoFires_dy09_bm17mbm09_500000_randomsample.pdf")
ggplot(dy09_bm17mbm09_df_s500000, aes(x=dy09_bm17mbm09)) + 
  geom_histogram(binwidth=2,color="black", fill="light blue") + 
  xlab("dy13_bm13mbm09") + 
  theme(axis.title.x = element_text(face="plain",color="black",size=16)) +
  ylab("Count") +
  xlim(c(-45,45)) +
  theme(axis.title.y = element_text(face="plain",color="black",size=16))  +
  #annotate("text", 50, -50, label="No Fires dy17_bm17mbm13 50000 random sample", color = "black",hjust = 1,size = 6)
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("No Fires dy09_bm17mbm09 500,000 random sample")
dev.off()




# Create histograms in loop - not used
for (i in rs){
  print(i)
  dfname <- paste0(i,"_df_s50000")
  sampledf <- as.data.frame(paste0(i,"_df")[sample(nrow(paste0(i,"_df")), 100), ])
  #plot.new()
  #pdf(paste0(getwd(),"\\hist_",i,".pdf"))
  ggplot(dy17_bm17mbm13_df_s50000,aes(x=dy17_bm17mbm13_df_s50000[,"dy17_bm17mbm13"])) + 
    geom_histogram(binwidth=2,color="black", fill="white")

  
}





























# -- http://www.countbio.com/web_pages/left_object/R_for_biology/R_fundamentals/3D_surface_plot_R.html ------

##  2D Gaussian Kernal plot

## Generate x and y coordinates as sequences
x = seq(-4,4,0.2)
y = seq(-4,4,0.2)

# An empty matrix z
z = matrix(data=NA, nrow=length(x), ncol=length(x))

### Gaussian kernal generation to fill the z matrix.
sigma = 1.0
mux = 0.0
muy = 0.0
A = 1.0

for(i in 1:length(x)){
  for(j in 1:length(y))
  {
    
    z[i,j] = A * (1/(2*pi*sigma^2)) * exp( -((x[i]-mux)^2 + (y[j]-muy)^2)/(2*sigma^2)) 
  }
}

### Now z is a matrix of dimension length(x) * length(y)

# Plotting surface with persp() function of default "Graphics" library in R.
## (Note: At this point, just give any x,y vector and z matrix of terrain to 
##       get your plot)

persp(x,y,z, theta=30, phi=50, r=2, shade=0.4, axes=TRUE,scale=TRUE, box=TRUE, 
      nticks=5, ticktype="detailed", zlim = c(0,0.2) , col="cyan", xlab="X-value", 
      ylab="Y-value", zlab="Z-value", main="Gaussian Kernal with persp()")
typeof(z) #double


# --- https://stats.stackexchange.com/questions/197455/how-to-plot-3d-partial-dependence-in-gbm -----
library(gbm)
library(reshape2)
data(diabetes, package = 'lars')

y        <- diabetes$y
x        <- diabetes$x
head(x)
class(x) <- 'matrix'
data     <- data.frame(y, as.data.frame(x))

gbm.model <- gbm::gbm(formula = y ~ . , data = data, distribution = 'gaussian', 
                      shrinkage = 1, bag.fraction = 1, n.trees = 100,
                      interaction.depth = 3, verbose = T, keep.data = F)


partial <- plot(gbm.model, i.var = c(1,2), return.grid = T)

colnames(partial)

mat <- reshape2::acast(data = partial, formula = age ~ sex, value.var = 'y')

persp(x = as.numeric(colnames(mat)), y = as.numeric(rownames(mat)), z=mat,
      zlab = 'partial dependence', xlab = 'sex', ylab = 'age', theta = 30)




# -----Other attempts at plotting random forest info -----------------------------

rfsrc_1 <- rfsrc(heb_2011 ~ NDIstd1103 + NBRstd1103 + ARGmax1106 + B3max1106, data=y11_nofires)
varsel_1 <- var.select(rfsrc_1)
# Plot the OOB errors against the growth of the forest.
#gg_e <- gg_error(rfsrc_1)
#plot(gg_e)
# Plot predicted median home values.
plot(gg_rfsrc(rfsrc_1), alpha=.5)+ coord_cartesian(ylim=c(5,49))
plot(gg_vimp(rfsrc_1), lbls=pltnum)
varsel_1 <- var.select(rfsrc_1)
gg_md <- gg_minimal_depth(varsel_1)
plot(gg_md, lbls=pltnum)
plot(gg_minimal_vimp(gg_md))

gg_v <- gg_variable(rfsrc_1)
xvar <- gg_md$topvars

plot(gg_v, xvar=xvar, panel=TRUE,se=.95, span=1.2, alpha=.4) + labs(y=pltnum["medv"], x="")

partial_1 <- plot.variable(rfsrc_1, xvar=gg_md$topvars,partial=TRUE, sorted=FALSE,show.plots = FALSE )
gg_p <- gg_partial(partial_1)

plot(gg_p, xvar=xvar, panel=TRUE, se=FALSE) + labs(y=pltnum["medv"], x="")
# ------------ from randomForestExplainer ------------
min_depth_frame <- min_depth_distribution(rf)
head(min_depth_frame)

# plot_min_depth_distribution(forest) # gives the same result as below but takes longer
plot_min_depth_distribution(min_depth_frame)

importance_frame <- measure_importance(rf)
importance_frame

# plot_multi_way_importance(forest, size_measure = "no_of_nodes") # gives the same result as below but takes longer
plot_multi_way_importance(importance_frame, size_measure = "no_of_nodes")

plot_multi_way_importance(importance_frame, x_measure = "mse_increase", y_measure = "node_purity_increase", size_measure = "p_value", no_of_labels = 5)
#plot_importance_ggpairs(rf) # gives the same result as below but takes longer
plot_importance_ggpairs(importance_frame)

# plot_importance_rankings(forest) # gives the same result as below but takes longer
plot_importance_rankings(importance_frame)

(vars <- important_variables(importance_frame, k = 5, measures = c("mean_min_depth", "no_of_trees")))
interactions_frame <- min_depth_interactions(rf, vars)
head(interactions_frame[order(interactions_frame$occurrences, decreasing = TRUE), ])
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Determine the marginal effect of a variable on the response
partialPlot(rf, y15_nofires, ARGmax1506, "versicolor")
plot(rf)
detach()
attach(y17_nofires)
c <- lm(heb_2017 ~ ARGmea1706 + ARGstd1706 + B3min1706)
summary(c); vif(c)
plot(c, which = 2)
model.diag.metrics <- augment(c)
#head(model.diag.metrics)
coeff=coefficients(c)
y15_nofires$fits <- fitted.values(c)
eq = "ARGmax1503 + ARGmax1506"
#plot.new()
#pdf("NBRstd0903ARGmed0906B3max0909_y09_nofires.pdf") 
ggplot(y15_nofires, aes(x=heb_2015, y=fits)) + # + mam_NDVI_max + mam_RGA_sd,colour=Herb_PCT)) +
  geom_smooth(method = lm, se = FALSE) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  ggtitle(eq)
#geom_text(aes(label=NewID)) 
#dev.off() 

library(party)
crf <- cforest(heb_2011 ~ NDIstd1103 + NBRstd1103 + ARGmax1106 +      B3max1106, data=y11_nofires, controls=cforest_control(mtry=2, mincriterion=0))
summary(crf)

# ---Extras in betareg -----------------------------------------------------------------------

detach()
detach(y18_nofires)
attach(y18_nofires)
lmodel <- lm( hpsc ~ NBRmea1103 + NDImin1103 + B3med1106 + NDImed1103)
summary(lmodel); vif(lmodel)
ggplot(y11_nofires, aes(x=heb_2011, y= NDVstd1103+ARGmea1106 ))+
  geom_smooth(method = lm, se = FALSE) +
  geom_point() + geom_text(aes(label=pltnum))


plot(beta_model)


pdf("betareg_jja_RGA_mean_mam_NDVI_max_mam_RGA_sd.pdf")
plot.new()
#par(mfrow = c(3, 2))
plot(beta_model2, which = 1:4)
dev.off() 

betareg_predict <- cbind(
  predict(beta_model2, type = "response"),
  predict(beta_model2, type = "link"),
  predict(beta_model2, type = "precision"),
  predict(beta_model2, type = "variance"),
  predict(beta_model2, type = "quantile", at = c(0.25, 0.5, 0.75))
)
colnames(betareg_predict) <- c("response","link","precision","variance","quantilep25","quantilep5","quantilep75")
head(betareg_predict)

beta_model2 <- betareg(hpsc ~ ARGmea1506 + NDVmax1503 + ARGstd1503, na.action=na.exclude, data=y15_nofires)
summary(beta_1)
par(mfrow = c(1,1))
plot(beta_1, which = 1:4, id.n = 10, labels.id = names(residuals(beta_1)))

aa <- cbind(betareg_predict,y15_nofires)

#pdf("betareg_0_1_jja_RGA_mean_mam_NDVI_max_mam_RGA_sd.pdf")
#plot.new()
#  plot betareg output w ggplot
ggplot(aa, aes(x=hpsc, y=response)) +
  geom_point() + geom_text(aes(label=pltnum)) +
  scale_fill_grey() +
  abline(a=0,b=1) +
  #geom_line(aes(y = predict(beta_loglog, aa),colour = "log-log", linetype = "log-log")) +
  #geom_line(aes(y = predict(beta_model2, y15_nofires), colour = "logit", linetype = "logit")) +
  scale_colour_manual("", values = c("red", "blue")) +
  scale_linetype_manual("", values = c("solid", "dashed")) +
  geom_smooth(method = "lm", se = TRUE,colour = "red") #+
  abline(a=0,b=1)
theme_bw()


plot(y15_nofires$hpsc, fitted(beta_model2), ylim = c(0,1), xlim = c(0,1),main = "probit")
plot(MMI_4$dcc_prop, fitted(beta_model1), ylim = c(0,0.7), xlim = c(0,0.7), main = "beta model")
abline(a=0,b=1)

mat <- incsv[,c(33,56,67)]
#bfit <- betareg.fit(mat,incsv$new1,na.action=na.exclude)


nofires2$fits <- fitted.values(beta_model2)
plot.new()
pdf("Herb_PCT_and_jja_RGA_mean_mam_NDVI_max_mam_RGA_sd__betareg.pdf") 
eq = "jja_RGA_mean_mam_NDVI_max_mam_RGA_sd"
ggplot(nofires2, aes(x=Herb_PCT, y=fits)) + # + mam_NDVI_max + mam_RGA_sd,colour=Herb_PCT)) +
  geom_smooth(method = lm, se = FALSE) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  ggtitle(eq) +
  geom_text(aes(label=NewID)) 
dev.off() 


logistic <- function(p) log(p / (1-p) +0.01)
lm(logistic(new1)~jja_RGA_mean + mam_NDVI_max + mam_RGA_sd, na.action=na.exclude,data=incsv)




# ---------------  Mosaic rasters
# Run in ShrubGrass_extractededartwithopints.R


# ----------- Predict out new rasters
# -------------- cortest based on pared down table - try later  --------------
newtab <- y13_nofires[c(1,17,46,62)]
colnames(newtab)

cormatrix <- corr.test(newtab[1:ncol(newtab)], use="pairwise",method="pearson",adjust="none")
round(cormatrix$r,2)
#cor.plot(cormatrix)
