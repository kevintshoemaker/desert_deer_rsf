
#################
# Question: is there evidence for a tradeoff whereby females that hang out in safer (but less resource-favorable) places are more likely 
#     to have higher current reproductive success, and females that favor resource acquisition are 
#     prioritizing future reproduction (giving up on the present year)?

rm(list=ls())

###############
# TODO
###############

### rerun for different background ratios and see if coefficients and likelihood stabilize

### run with habitat type as a covar?
### make sure we account for the fact that NDVI changes with time (ndvi differs for the same location by month)
### use a 

#################
# Load packages
#################

library(lme4)
library(glmmTMB)
library(runjags)
library(loo)
library(survival)
library(lubridate)
library(caret)
library(ranger)
#library(buildmer)

#library(ppmlasso)

######################
# LOAD FUNCTIONS
######################

predvar="water"
varname = "dist2water"
status="postmort"
VisualizeRelation <- function(data=df_rsf2,model=mod_final2,predvar="NDVI",varname="NDVI",allvars=predvars,status="atheel"){
  len <- 25
  
  dataclasses <- sapply(data,class)
  
  dim <- data[,predvar]
  range <- seq(min(dim),max(dim),length=len)
  
  realmean <- mean(df_rsf[[predvar]])
  realsd <- sd(df_rsf[[predvar]])
  
  newdata <- data.frame(temp=range)
  names(newdata) <- c(predvar)
  
  othervars <- allvars[!allvars%in%c(predvar,"USED")]
  
  tmp <- table(data$ID_YR,data$STATUS)[,status]
  IDs <- names(tmp)[tmp>100]
  
  
  var = othervars[5]
  for(var in othervars){
    thisvar <- data[[var]]
    if(is.factor(thisvar)){
      tab <- table(thisvar)
      vals <- names(tab)
      levs <- levels(thisvar)
      mostcom <- vals[which.max(tab)]
      newvec <- factor(rep(mostcom,times=nrow(newdata)),levels=levs)
      newdata[,var] <- newvec
    }else{
      newdata[,var] <- 0 #mean(thisvar)
    }
  }
  
  newdata$ID_YR <- NA  #data$ID[1]
  newdata$weights = 1
  newdata$STATUS = factor(stat,levels=levels(data$STATUS))
  # pred <- predict(model,newdata[1,])   # ,type="conditional",se.fit=F
  
  x <- getME(model,"X")
  # str(x)
  # colnames(x)
  # x[1:2,]
  newdata2 <- x[1:len,]
  newdata2[]  <- 0
  newdata2[,"(Intercept)"] <- 1
  #newdata[1:2,]
  newdata2[,predvars] <-  as.matrix(newdata[,predvars]) 
  if(status=="atheel"){
    newdata2[,"STATUSatheel"] <- 1
    for(p in predvars){
      thisvar <- paste0(p,":STATUSatheel")
      newdata2[,thisvar] <- newdata[,p]
    }
    
  }else if(status=="postmort"){
    newdata2[,"STATUSpostmort"] <- 1
    for(p in predvars){
      thisvar <- paste0(p,":STATUSpostmort")
      newdata2[,thisvar] <- newdata[,p]
    }
  }
  
  rr <- ranef(model)
  rr0 <- unlist(rr$cond)
  #str(rr0)
  eta <- c(matrix(newdata2 %*% fixef(model)$cond )) #+ z %*%  rr0))
  mu <- plogis(eta)
  pred = mu
  
  z <- getME(model,"Z")
  zshort <- z[1:len,]
  # str(z)
  # colnames(z)
  pred2 <- matrix(0,nrow=len,ncol=40)
  id=1
  for(id in 1:min(30,length(IDs))){
    #newdata$ID_YR <- IDs[id]  #unique(data$ID_YR)[id]
    
    newran2 <- zshort
    newran2[]  <- 0
    newran2[,IDs[id]] <- 1
    #colnames(newran2)
    
    statid <- paste0(status,":",IDs[id])  
    newran2[,statid] <- 1
    
    #str(rr0)
    eta <- c(matrix(newdata2 %*% fixef(model)$cond + newran2 %*%  rr0))
    mu <- plogis(eta)
    pred2[,id] = mu
    
    #pred2[,id] <- predict(model,newdata,type="conditional",se.fit=F)
    #lines(range,pred2[,id]*1000,col=cols[id],lwd=1)
  }
  
  yrge <- diff(range(pred))*1000*2
  plot(range,pred*1000,xlab=varname,ylab="Sel. Intensity",type="l",lwd=2,xaxt="n",col="black",
       ylim=c(max(min(pred2)*1000,min(pred)-yrge),min(max(pred)*1000+yrge,max(pred2)*1000)))
 # points(range,(pred$fit*1000+pred$se.fit*1000),type="l",lty=2)
 #  points(range,(pred$fit*1000-pred$se.fit*1000),type="l",lty=2)
  ats <- seq(min(range),max(range),length=6)
  axis(1,ats,labels = round(realmean+ats*realsd,ifelse(predvar%in%c("elevation","water"),0,2)))
  rug(jitter(data[seq(1,nrow(data),10000),][[predvar]]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
  id=1
  cols <- gray.colors(10,0.3,0.8)
  for(id in 1:min(30,length(IDs))){
    # newdata$ID_YR <- data$ID_YR[id]
    # pred <- predict(model,newdata,type="conditional",se.fit=F)
    lines(range,pred2[,id]*1000,col="gray",lwd=1)
  }
  
  lines(range,pred*1000,col="black",lwd=2)
}

VisualizeRelation_rf <- function(data,model,predvar,varname,allvars,status="atheel"){
  len <- 25
  
  dataclasses <- sapply(data,class)
  
  dim <- data[,predvar]
  range <- seq(min(dim),max(dim),length=len)
  
  realmean <- mean(df_rsf[[predvar]])
  realsd <- sd(df_rsf[[predvar]])
  
  newdata <- data.frame(temp=range)
  names(newdata) <- c(predvar)
  
  othervars <- allvars[!allvars%in%c(predvar,"USED")]
  
  var = othervars[5]
  for(var in othervars){
    thisvar <- data[[var]]
    if(is.factor(thisvar)){
      habvars <- model$finalModel$xNames[grepl("habitat",model$finalModel$xNames)]
      for(v2 in 1:length(habvars)){
        newdata[[habvars[v2]]] <- 0
      }
    }else{
      newdata[,var] <- 0 #mean(thisvar)
    }
  }
  
  #model$finalModel$xNames[grepl("habitat",model$finalModel$xNames)]
  
  # newdata$USED <- factor("0",levels=c("0","1"))
  # mf <- model.frame(paste0("USED~",paste(othervars,collapse="+")),newdata)
  # model.matrix(as.formula(paste0("USED~",paste(othervars,collapse="+"))),mf)
  # 
  library(ranger)
  preds <- predict(model$finalModel,data=newdata)
  
  pred <- preds$predictions[,2]
  
  yrge <- diff(range(pred))
  plot(range,pred,xlab=varname,ylab="Sel. Intensity",type="l",lwd=2,xaxt="n",col="black",
       ylim=c(max(min(pred),min(pred)-yrge),min(max(pred)+yrge,max(pred))),main=paste0(predvar,", ",status))
  # points(range,(pred$fit*1000+pred$se.fit*1000),type="l",lty=2)
  #  points(range,(pred$fit*1000-pred$se.fit*1000),type="l",lty=2)
  ats <- seq(min(range),max(range),length=6)
  axis(1,ats,labels = round(realmean+ats*realsd,ifelse(predvar%in%c("elevation","water"),0,2)))
  rug(jitter(data[seq(1,nrow(data),1000),][[predvar]]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
  id=1
  cols <- gray.colors(10,0.3,0.8)
  
  
}



######################
# GLMM analysis
#######################

#rm(list=ls())

##############
# data from levi ms thesis

      # takes a while to read with 50:1 background to used ratio
df_rsf <- read.csv("MOJA_Deer_repro_50-1.csv",stringsAsFactors = F)  #     All_repro_use_avail_data_ch2.csv  newer: all_repro_data_with_dates_clean.csv
names(df_rsf)

table(df_rsf$STATUS)
df_rsf$STATUS <- factor(df_rsf$STATUS,levels=c("prepart","atheel","postmort"))
#levels(df_rsf$STATUS) <- c("prepart","atheel","postmort")
df_rsf$habitat <- factor(df_rsf$habitat,levels=c("Joshua Tree Woodland","Big Sagebrush","Burn Recovery","Juniper Woodland","Low Shrub","Pinyon Woodland","Wash","Yucca Shrub"))
#levels(df_rsf$habitat) <- c("Joshua Tree Woodland","Big Sagebrush","Burn Recovery","Juniper Woodland","Low Shrub","Pinyon Woodland","Wash","Yucca Shrub")

table(df_rsf$habitat[df_rsf$USED==1])

#summary(df_rsf)

# head(df_rsf)
names(df_rsf)
# df_rsf$mojaveg_ra   # habitat type
# nrow(df_rsf)

allcovars <- c(
  "habitat",
  "rugged",
  "water",
  "treecvr",
  "slope",
  "sine",
  "shrubcvr",
  "elevation",
  "cosine",
  "NDVI"
)

tostd <- c(
  "rugged",
  "water",
  "treecvr",
  "slope",
  "sine",
  "shrubcvr",
  "elevation",
  "cosine",
  "NDVI"
)

###########
# process dates
###########

df_rsf$DATE2 <- ymd(df_rsf$DATE)
df_rsf$YEAR <- year(df_rsf$DATE)
df_rsf$day <- yday(df_rsf$DATE2) 

## save original df

save(df_rsf,file="originaldf.RData")

#load("originaldf.RData")

# temp <- lapply(tostd,function(t) {df_rsf[[t]] <<- scale(df_rsf[[t]])}   )  # standardize
# 
# save(df_rsf,file="stddf.RData")

#lapply(tostd,function(t) hist(df_rsf[[t]],main=t))

nback <- length(which(df_rsf$USED==0))   # 3.5m background
nused <- length(which(df_rsf$USED==1))   # 68166 used
nback/nused   # 50 background per used


areas = read.csv("HR_area_200m_buff.csv",stringsAsFactors = F)
siteareas = areas$area_km2
names(siteareas) = areas$ID




#df_rsf$ID_YR <- paste(df_rsf$ID,df_rsf$YEAR,sep="_")
#temp <- subset(df_rsf,ID_YR=="112_2016")    # fixes every 1.5 hours?


##############
# Correlation analysis
##############

predvars <- c("rugged","water","NDVI","elevation","slope")

cor(df_rsf[,predvars])

predvars <- c("rugged","water","NDVI","elevation")   # remove ruggedness or slope



############
# Identify how many background points are needed!
############

# do it one by one for each animal?

# ?ppmlasso::findres

# data("BlueMountains")
# BlueMountains$env
# rm(BlueMountains)

allinds <- unique(df_rsf$ID_YR)
thisind <- allinds[1]
thisindd <- strsplit(thisind,"_")[[1]][1]

thisdf <- subset(df_rsf,ID_YR==thisind)
nrow(thisdf)

allmon <- unique(thisdf$month)
thismon <- allmon[2] 
thisdf2 <- subset(thisdf,month==thismon)
nrow(thisdf2)

thisdf2 <- subset(thisdf,STATUS="prepart")

# controls for thinning data
#perobs <- 100
perday <- 6

## thin the dataset for model testing (all individuals)
#thisdf3 <- thisdf
thisdf3 <- thisdf2[rep(FALSE,times=nrow(thisdf2)),]

rs="prepart"
for(rs in levels(thisdf2$STATUS)){
  temp2 <- subset(thisdf2,STATUS==rs)
  if(nrow(temp2)>0){
    #months <- sort(unique(temp2$month))
    days <- sort(unique(temp2$day))
    d=days[2]
    for(d in days){
      temp3 <- subset(temp2,day==d)
      temp4 <- subset(temp3,USED==1)
      usdpt <- temp4[sort(sample(1:nrow(temp4),perday)),]
      #toadd <- rbind(usdpt,rdmpt)
      thisdf3 <- rbind(thisdf3,usdpt)
    }
  }
}
nrow(thisdf3)
thisdf3 <- rbind(thisdf3,thisdf[thisdf$USED==0,])

quad = thisdf3[thisdf3$USED==0,c("UTM_X","UTM_Y",predvars,"habitat","USED")]   # change to thisdf3?
names(quad)[1:2] = c("X","Y")
nrow(quad)
nrow(thisdf3[thisdf3$USED==1,])

# quad$X <- round(quad$X/10)*10
# quad$Y <- round(quad$Y/10)*10
# quad <- quad[order(quad$X,quad$Y),]
# tail(quad)
# summary(quad)
# nrow(quad)

#temp <- sample.quad(quad,10)
#nrow(temp)

# sp.xy = data.frame(X=thisdf3$UTM_X[thisdf$USED==1],
#                    Y=thisdf3$UTM_Y[thisdf$USED==1])
# nrow(sp.xy)
# ppm.form = ~ poly(rugged, water, NDVI, elevation, degree = 2, raw = TRUE) #+ habitat

#ppm.form = ~rugged+water+NDVI+elevation #+habitat

    # wow, this takes a long time, even for smallish datasets...
# system.time(
#   mod <- ppmlasso(formula=ppm.form,sp.xy=sp.xy,env.grid=quad,sp.scale=20,lamb=0)
# )
# summary(mod)
# mod$betas
# mod$lambdas
# mod$likelihoods
# mod

# scales = c(10,20,40,50,60,80,100,150)
# fr <- findres(scales, formula = ppm.form, sp.xy = sp.xy, env.grid = quad)


###############
# use downweighted poisson regression?

n.quad = c(50, 100, 250, 500, 1000, 1500, 2000, 2500)   # per km2
sapply(1:length(n.quad),function(t) n.quad[t]*siteareas)

reps = 25
r=1
for(r in 1:reps){
  
  #n.quad2 = pmin(n.quad,nrow(quad))
  n.quad2 = round(n.quad*siteareas[thisindd])
  
  nrow(quad)
  
  prez <- thisdf3[thisdf3$USED==1,c("UTM_X","UTM_Y",predvars,"habitat","USED")]   # change to thisdf?
  names(prez)[1:2] = c("X","Y")
  nrow(prez)
  
  # quad.inc = sample(1:dim(quad)[1], n.quad2[1])
  # assign(paste("quad.", n.quad[1], sep = ""), quad[quad.inc[1:n.quad[1]],])
  # i=2
  # for (i in 2:length(n.quad2)){
  #   quad.inc = c(quad.inc, sample(setdiff(1:dim(quad)[1], quad.inc),(n.quad2[i] - n.quad2[i - 1])))
  #   assign(paste("quad.", n.quad2[i], sep = ""), quad[quad.inc[1:n.quad2[i]],])
  # }
  
  for (i in 1:length(n.quad2)){
    quad.inc = sample(1:nrow(quad),n.quad2[i],replace = T)
    assign(paste("quad.", n.quad2[i], sep = ""), quad[quad.inc,])
  }
  
  
  loglik = rep(NA, length(n.quad2))
  i=1
  for (i in 1:length(n.quad2)){
    thisquad = get(paste("quad.", n.quad2[i], sep = ""))
    all.dat = data.frame(rbind(prez, thisquad))
    X.des = as.matrix(cbind(poly(all.dat$rugged, all.dat$water, all.dat$NDVI,
                                 all.dat$elevation, degree = 2, raw = TRUE), all.dat$habitat))
    p.wt = rep(1.e-8, dim(all.dat)[1])
    p.wt[all.dat$USED == 0] = siteareas[thisindd]/n.quad2[i]
    z = all.dat$USED/p.wt
    dwpr = glm(z ~ X.des, family = poisson(), weights = p.wt)
    mu = dwpr$fitted
    loglik[i] = sum(p.wt*(z*log(mu) - mu))
  }
  if(r==1){
    min = loglik[length(n.quad)]-abs(min(loglik)-max(loglik))*3
    max = loglik[length(n.quad)]+abs(min(loglik)-max(loglik))*3
    plot(n.quad, loglik, type = "o",ylim=c(min,max))  #log = "x", 
  } else{
    lines(n.quad, loglik, type = "o")
  }
  
}


### Result: 1000 random points per km2 works well
### 500 points works pretty well too! USE 500 POINTS


rm(list=ls()[grepl("quad\\.",ls())])




##############
# THIN the dataset appropriately
##############



# temp <- table(df_rsf$ID_YR,df_rsf$STATUS)
# full_inds <- names(apply(temp,1,function(t) all(t>200) ))[apply(temp,1,function(t) all(t>200) )]
# # 20 ids with good data for all three phases
# 
# 
# df_rsf3 <- df_rsf[df_rsf$ID_YR%in%full_inds,]   # subset for just individuals with good amt of data for all phases
# nrow(df_rsf3)   #1.1m rows

# ## thin the dataset for model testing (only those individuals with data for all phases)
# 
# df_rsf4 <- df_rsf3
# df_rsf4 <- df_rsf4[rep(FALSE,times=nrow(df_rsf4)),]
# id="414_2014"
# for(id in unique(df_rsf3$ID_YR)){
#   temp <- subset(df_rsf3,ID_YR==id)
#   rs="atheel"
#   for(rs in levels(df_rsf3$STATUS)){
#     temp2 <- subset(temp,STATUS==rs)
#     df_rsf4 <- rbind(df_rsf4,temp2[ceiling(seq(0.1,nrow(temp2),length=min(750,nrow(temp2)))),])
#   }
# }
# nrow(df_rsf4)   # 41318 observations
# 
# unique(df_rsf4$ID_YR)
# 
# 
# ## thin the dataset for model testing (all individuals)
# df_rsf2 <- df_rsf
# df_rsf2 <- df_rsf2[rep(FALSE,times=nrow(df_rsf2)),]
# id="414_2014"
# for(id in unique(df_rsf$ID_YR)){
#   temp <- subset(df_rsf,ID_YR==id)
#   rs="atheel"
#   for(rs in levels(df_rsf$STATUS)){
#     temp2 <- subset(temp,STATUS==rs)
#     if(nrow(temp2)>0){
#       df_rsf2 <- rbind(df_rsf2,temp2[ceiling(seq(0.1,nrow(temp2),length=min(500,nrow(temp2)))),]) 
#     }
#   }
# }
# nrow(df_rsf2)   # 76036 observations


# note: 150 and 4 works
# 500 and 10 does not work
# 250 and 6 works (keep this one)

load(file="originaldf.RData")

   # controls for thinning data
perkm2 <- 500   #250    # should be 1000 or 500
perday <- 18      # should be 6 or 10?

## thin the dataset for model testing (all individuals)
df_rsf2 <- df_rsf
df_rsf2 <- df_rsf2[rep(FALSE,times=nrow(df_rsf2)),]
df_blank <- df_rsf2

id="414_2014"
for(id in unique(df_rsf$ID_YR)){
  thisindd <- strsplit(id,"_")[[1]][1]
  temp <- subset(df_rsf,ID_YR==id)
  nrand <- round(perkm2*siteareas[thisindd])
  
  allmon <- as.character(unique(temp$month))
  
  bgs <- list()
  
  thismon <- allmon[2] 
  for(thismon in allmon){
    bgs[[thismon]] <- subset(temp,month==thismon&USED==0)
    #nrow(bg1)
  }
  
  toadd_o <- df_blank
  toadd_bg <- df_blank
  rs="atheel"
  for(rs in levels(df_rsf$STATUS)){
    temp2 <- subset(temp,STATUS==rs&USED==1)
    if(nrow(temp2)>0){
      #months <- sort(unique(temp2$month))
      days <- sort(unique(temp2$day))
      d=days[2]
      for(d in days){
        temp3 <- subset(temp2,day==d)
        usdpt <- temp3[sort(sample(1:nrow(temp3),min(nrow(temp3),perday))),]
        toadd_o <- rbind(toadd_o,usdpt)
      }
      #toadd_o
      montab <- table(toadd_o$month)/nrow(toadd_o)
      m=names(montab)[1]
      for(m in names(montab)){
        thisrdm <- round(nrand*as.numeric(montab[m]))
        bgpt <- bgs[[m]][sample(1:min(thisrdm,nrow(bgs[[m]]))),]
        bgpt$STATUS = rs 
        toadd_bg<- rbind(toadd_bg,bgpt)   
      }
    }
  }
  toadd <- rbind(toadd_o,toadd_bg)
  df_rsf2 <- rbind(df_rsf2,toadd)
}

nrow(df_rsf2)   # 1.9 million observations with 500 background points

nrow(df_rsf)   # 3.5 million observations for the 50:1 dataset!

save(df_rsf2,file="df_forRF.RData")

temp <- lapply(tostd,function(t) {df_rsf2[[t]] <<- scale(df_rsf2[[t]])}   )  # standardize

save(df_rsf2,file="df_forfit.RData")



# hist(df_rsf2$UTM_X)
# hist(df_rsf2$UTM_Y)
# hist(df_rsf2$cosine)
# hist(df_rsf2$elevation)
# 
# hist(df_rsf2$NDVI)
# hist(df_rsf2$elevation)
# hist(df_rsf2$water)
# hist(df_rsf2$treecvr)
# 
# hist(df_rsf2$rugged)
# hist(df_rsf2$shrubcvr)
# table(df_rsf2$habitat)
# table(df_rsf$YEAR)
# table(df_rsf2$USED)     # 3146472/42376  74 random points per used point
# table(df_rsf2$ID_YR)    



###############
# Select important features using CARET package
###############

load("df_forRF.RData")
df_rsf <- df_rsf2

load("df_forfit.RData")
nrow(df_rsf)
nrow(df_rsf2)

df2 <- df_rsf2

summary(df2)

names(df2)


  ##############
  # Correlation analysis

cor(df_rsf[,tostd])   # remove 'slope'

   ########
   # prepare data

allcovars <- allcovars[!allcovars=="slope"]  # remove slope
allcovars

resp.var <- "USED"

df2$USED <- as.factor(df2$USED)
levels(df2$USED) <- c("UNUSED","USED")

levels(df2$habitat) <- c("josh","bigsage","burn","juniper","lowshrub","pinyon","wash","yucca")

#df2$habitat

  #######
  ## run caret

all.vars <- c(allcovars,resp.var)

formula <- as.formula(sprintf("%s~%s",resp.var,paste(allcovars,collapse="+")))

# prepare training scheme
#control <- caret::trainControl(method="repeatedcv", number=10, repeats=3)

#control <- caret::trainControl(method="repeatedcv", number=10, repeats=3,classProbs=TRUE,summaryFunction = twoClassSummary)

# thisndx <- groupKFold(df2$ID_YR,3)
# #thisndx$Fold1
# 
# control <- caret::trainControl(method="cv", number=3,classProbs=TRUE,summaryFunction = twoClassSummary,index = thisndx)

bestTune <- data.frame(mtry=5,splitrule="gini",min.node.size=1)

# train the model
# model <- caret::train(formula, data=df2, method = 'ranger', 
#                       trControl=control, importance = 'impurity',tuneGrid=bestTune)   # method="ranger"

#save(model,file="RFcaretModel_alldata.RData")

# estimate variable importance
# importance <- caret::varImp(model, scale=FALSE)
# 
# save(model, importance, file="RFcaretModel_alldata.RData")
# 
# summary(model)
# 
# model$bestTune    # mtry is 8, splitrule is gini [better mtry is probably 5]
# model$results
# 
# model$finalModel   # this gives the final model
# 
# modcall <- model$finalModel$call

########
#   loop through statuses


allstat <- levels(df2$STATUS)

modlist <- list()
importlist <- list()
stat <- 1
for(stat in 1:length(allstat)){
  dfnew <- subset(df2,STATUS==allstat[stat])
  thisndx <- groupKFold(dfnew$ID_YR,3)
  #thisndx$Fold1
  
  control <- caret::trainControl(method="cv", number=3,classProbs=TRUE,summaryFunction = twoClassSummary,index = thisndx)
  
  bestTune <- data.frame(mtry=5,splitrule="gini",min.node.size=1)
  modlist[[allstat[stat]]] <- caret::train(formula, data=dfnew, method = 'ranger', 
                                           trControl=control, importance = 'impurity',tuneGrid=bestTune)   # method="ranger"
  importlist[[allstat[stat]]] <- caret::varImp(modlist[[allstat[stat]]], scale=FALSE)
}
save(modlist, importlist, file="RFcaretModel_bystatus_indcv.RData")
load("RFcaretModel_bystatus_indcv.RData")

load("RFcaretModel_bystatus.RData")

stat <- 1
for(stat in 1:length(allstat)){
  importlist[[allstat[stat]]]
  #svg(sprintf("importance_%s.svg",allstat[stat]),3,4)
  print(plot(importlist[[allstat[stat]]],main=allstat[stat]))
  #dev.off()
}

stat <- 1
for(stat in 1:length(allstat)){
  print(modlist[[allstat[stat]]])
}


# summarize importance
print(importance)
# plot variable importance
plot(importance)  


  ########
  # univariate plots

dflist <- list()
for(stat in 1:length(allstat)){
  dflist[[allstat[stat]]] <- subset(df2,STATUS==allstat[stat])
}

statuses <- levels(df2$STATUS)
vars <- c("rugged","water","treecvr","cosine","shrubcvr","elevation","NDVI")
stat <- statuses[1]
var <- vars[2]

for(stat in statuses){
  for(var in vars){
    thismod <- modlist[[stat]]
    thisdf <- dflist[[stat]]
    thisterms <- rownames(attributes(thismod$terms)$factors)
    
    svg(sprintf("univarplot_%s_%s.svg",stat,var),4,3)
    VisualizeRelation_rf(data=thisdf,model=thismod,predvar=var,varname=var,allvars=thisterms,status=stat)
    dev.off()
  }
}


########
# Feature selection

# 
# define the control using a random forest selection function

stat <- statuses[1]

thismod <- modlist[[stat]]
thisdf <- dflist[[stat]]
thisndx <- groupKFold(thisdf$ID_YR,3)
control <- caret::rfeControl(functions=rfFuncs, method="cv", number=3,index = thisndx)


thisterms <- rownames(attributes(thismod$terms)$factors)


tick <- Sys.time()
# run the RFE algorithm
results <- caret::rfe(thisdf[,vars], thisdf$USED, sizes=c(1,2,3,4,5,6,7,9), rfeControl=control)
tock <- Sys.time()

save(results,file="rfe_results_prepart.RData")


# summarize the results
print(results)    
# list the chosen features
predictors(results)
# plot the results
plot(results, type=c("g", "o"))
results$optVariables








##############
# Build GLMM_TMB models 
##############




# 
# summary(df_rsf2)
# head(df_rsf2,20)
# 
# length(unique(df_rsf2$ID_YR))    # 68 individuals
# table(df_rsf2$ID_YR,df_rsf2$STATUS)

#### FULL MODEL (find fullest model that converges without warnings)

# mod_comp2.0 <- glmmTMB(USED ~ STATUS + NDVI + elevation + Rugged + dist2WATER +    # doesn't converge
#                          NDVI:STATUS + elevation:STATUS + Rugged:STATUS + dist2WATER:STATUS + 
#                          (1|ID) + #(0+STATUS|ID) +
#                          (0 + NDVI | ID) + (0 + elevation | ID) + (0 + Rugged | ID) + (0 + dist2WATER | ID) +
#                          (0 + NDVI:STATUS | ID) + (0 + elevation:STATUS | ID) + 
#                          (0 + Rugged:STATUS | ID) + (0 + dist2WATER:STATUS | ID),
#                        weights = weights,
#                        family=binomial,
#                        data=df_rsf4,
#                        REML=FALSE)
# summary(mod_comp2.0)
# fixef(mod_comp2.0)
# ranef(mod_comp2.0)

# mod_comp2.1 <- glmmTMB(USED ~ STATUS + NDVI + elevation + Rugged + dist2WATER + 
#                        NDVI:STATUS + elevation:STATUS + Rugged:STATUS + dist2WATER:STATUS + 
#                        (1|ID) + (0+STATUS|ID) +
#                        #(0 + NDVI | ID) + (0 + elevation | ID) + (0 + Rugged | ID) + (0 + dist2WATER | ID) +
#                        (0 + NDVI:STATUS | ID) + (0 + elevation:STATUS | ID) + 
#                        (0 + Rugged:STATUS | ID) + (0 + dist2WATER:STATUS | ID),
#                      weights = weights,
#                      family=binomial,
#                      data=df_rsf2,
#                      REML=FALSE)

# mod_comp2.2 <- glmmTMB(USED ~ STATUS + NDVI + elevation + Rugged + dist2WATER +    # this one works
#                          NDVI:STATUS + elevation:STATUS + Rugged:STATUS + dist2WATER:STATUS + 
#                          (1|ID) + #(0+STATUS|ID) +
#                          #(0 + NDVI | ID) + (0 + elevation | ID) + (0 + Rugged | ID) + (0 + dist2WATER | ID) +
#                          (0 + NDVI | STATUS:ID) + (0 + elevation | STATUS:ID) + 
#                          (0 + Rugged | STATUS:ID) + (0 + dist2WATER | STATUS:ID),
#                        weights = weights,
#                        family=binomial,
#                        data=df_rsf4,
#                        REML=FALSE)

# summary(mod_comp2.2)
# fixef(mod_comp2.2)
# ranef(mod_comp2.2)


#load("df_forfit.RData")

nrow(df_rsf2)

df_rsf2$ID2 <- as.factor(df_rsf2$ID_YR)

df_rsf2$weights <- ifelse(df_rsf2$USED == 1, 1, 1000)   # add weights variable



names(df_rsf2)

mod_comp2.3 <- glmmTMB(USED ~ NDVI + elevation + rugged + water + STATUS +   # this one works   this is the best model
                         NDVI:STATUS + elevation:STATUS + rugged:STATUS + water:STATUS + 
                         (1|ID_YR) +     #(0+STATUS|ID) +
                         (0 + NDVI | ID_YR/STATUS) + (0 + elevation | ID_YR/STATUS) + 
                         (0 + rugged | ID_YR/STATUS) + (0 + water | ID_YR/STATUS),
                       weights = weights,
                       family=binomial,
                       data=df_rsf2,    # try with only individuals with all levels. 
                       REML=FALSE)

summary(mod_comp2.3)
fixef(mod_comp2.3)
ranef(mod_comp2.3)


# do AIC model selection to select a reasonable reduced model
# bbmle::AICtab(cow.mean.nonz, cow.mean, cow.med, cow.max,weights=TRUE,mnames=c("Mean Nonz", "Mean", "Median", "Max") )

numre <- length(colnames(ranef( mod_comp2.3)$cond$ID_YR))+length(colnames(ranef( mod_comp2.3)$cond$`STATUS:ID_YR`)) -1

# fit the model using the Muff et al. 2019 modifications


# make the table of coefficients etc...  



TMBStruc = glmmTMB(USED ~ NDVI + elevation + rugged + water + STATUS +   # this one works   this is the best model
                     NDVI:STATUS + elevation:STATUS + rugged:STATUS + water:STATUS + 
                     (1|ID_YR) +     #(0+STATUS|ID) +
                     (0 + NDVI | ID_YR/STATUS) + (0 + elevation | ID_YR/STATUS) + 
                     (0 + rugged | ID_YR/STATUS) + (0 + water | ID_YR/STATUS),
                   weights = weights,
                   family=binomial,
                   data=df_rsf2,
                   REML=FALSE,
                   doFit=FALSE)



# Fix the standard deviation of the first random term, which is the (1|id) component
# in the above model equation
TMBStruc$parameters$theta[1] = log(1e3)

#TMBStruc$condList$reTrms$Zt@factors

# Tell glmmTMB not to change the first entry of the vector of variances,
# and give all other variances another indicator to make sure they can be freely estimated
TMBStruc$mapArg = list(theta=factor(c(NA,1:numre)))

mod_final = glmmTMB:::fitTMB(TMBStruc)   # takes a long time to run!
summary(mod_final)
mod_final$obj      # note: can't make predictions from this. Must re-run without fixed intercept effect

mod_final2 = mod_comp2.3    # for making predictions

#########

save( mod_final2, file= "rsf_finalmodels4.RData")  # mod_final,

# load("rsf_finalmodels4.RData")

summary(mod_final2)

#################
# MAKE TABLE OF COEFS

cis <- confint(mod_final2,component="cond",method="wald")   # takes a little while.    #,method="profile"
fes <- fixef(mod_final2)$cond
table1 <- data.frame(
  label1 = names(fes)[c(2,9,8, 3,11,10, 4,13,12, 5,15,14, 6,7)],
  label2 = c("NDVI","NDVI (post mort)","NDVI (provis)",
             "Elevation","Elevation (post mort)","Elevation (provis)",
             "Ruggedness","Ruggedness (post mort)","Ruggedness (provis)",
             "Distance to Water","Distance to Water (post mort)","Distance to Water (provis)",
             "Provisioning","Post mort")
)
table1
table1$coef <- fes[match(table1$label1,names(fes))]
t=table1$label1[1]
ndx <- sapply(table1$label1,function(t) which(sapply(rownames(cis),function(v) grepl(t,v)  ))[1] )
table1$lcb <- cis[ndx,1] 
table1$ucb <- cis[ndx,2]
table1

#rownames(table1) <- table1$label2
table1 <- table1[,-1]

table1[,c(2,3,4)] <- round(table1[,c(2,3,4)],3)
table1$sig <- ifelse(sign(sign(table1[3])*sign(table1[,4]))==1,"*","")[,1]
names(table1) <- c("Name","Coef","CI,lower","CI,upper","Signif")
table1$Signif

write.csv(table1,"table1_v4.csv")


##################
# Basic univariate plots
##################

# recover the original data!
#df_rsf <- read.csv("MOJA_Deer_repro_50-1.csv",stringsAsFactors = F)  #     All_repro_use_avail_data_ch2.csv

load(file="originaldf.RData")

predvars <- c("NDVI","elevation","rugged","water")
predvars2 <- c("NDVI","Elevation","Ruggedness","Dist to H20")
predvars
predvars2
rep_stats <- levels(df_rsf2$STATUS)   #[c(3,1,2)]
rep_stats
rep_stats2 <- c("Pre-parturition","Provisioning","Post-loss")



varndx = 1
var = predvars[varndx]
var2 = predvars2[varndx]

ordr <- data.frame(
  var=rep(predvars,each=3),
  var2=rep(predvars2,each=3),
  stat=rep(rep_stats,times=4),
  stringsAsFactors = F
)
ordr2 <- ordr
ordr2

svg("pdp_byrep5.svg",8.5,8)

laymat <- matrix(1:30,nrow=6,byrow = T)
laymat[6,] = 40
laymat[2:5,2] = 39
laymat[] <- as.numeric(as.factor(laymat))
laymat[6,3:5] <- 24

layout(laymat,widths = c(2,1.2,4,4,4),heights=c(1,4,4,4,4,1))          
#layout.show(24)
onleft <- c(6,10,14,18)
ontop <- c(1,2,3,4,5,23)
axlabs <- c(22,24)
#counter=1    
counter2=1
counter=1
panel=0
panel=panel+1
# panel=22
#sapply(1:21,function(t) plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch=""))
while(panel < 25){
  if(panel %in% onleft){
    par(mai=c(0,0,0,0))
    plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
    text(0,0.5, tools::toTitleCase(ordr2$var2[counter]),adj=0,cex=1.5)
    panel=panel+1
  }else if (panel %in% ontop){
    if(panel%in%c(1,2,23)){
      par(mai=c(0,0,0,0))
      plot(0,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      panel=panel+1
    }else{
      par(mai=c(0,0,0,0))
      plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      text(0.5,0.5, rep_stats[counter2],adj=0.5,cex=1.5)
      counter2=counter2+1
      panel=panel+1
    }
  }else if(!panel%in%axlabs){
    #par(mai=base$mai)
    par(mai=c(0.3,0.4,0.0,0.0))
    var <- ordr2$var[counter]
    stat <- ordr2$stat[counter]
    var2 <- ordr2$var2[counter]
    allmains = predvars
    
    VisualizeRelation(data=df_rsf2,model=mod_final2,predvar=var,varname=var2,allvars = allmains, status = stat)      # run and save partial dependence plots for all top variable
    
    counter <- counter+1
    panel=panel+1
    
  }else{
    if(panel==22){
      plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      text(0.5,0.5,"Selection Intensity",srt=90,adj=0.5,cex=1.5)
      panel=panel+1
    }else{
      plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="",pch="")
      text(0.5,0.5,"Environmental gradient",adj=0.5,cex=1.5)
      panel=panel+1
    }
  }
}

dev.off()


#################
#################


#################
# Figure illustrating lack of tradeoff between forage quality and offspring safety
#################


########
# look for 'need for tradeoff' in values of key covariables...

graphics.off()

svg(filename="tradeoffs_v2.svg",width=7,height=4)

# [load up 50:1 dataset for this!]

df <- subset(df_rsf,USED==0)

rugsel <- ranef(mod_final2)$cond$ID[,"rugged"]
names(rugsel) <-  rownames(ranef(mod_final2)$cond$ID)

# ndvisel <- ranef(mod_final)$cond$ID[,"NDVI"]
# names(rugsel) <-  rownames(ranef(mod_final)$cond$ID)

ndviav <- sapply(names(rugsel),function(t) mean(subset(df,ID_YR==t)$NDVI)  )
rugav <- sapply(names(rugsel),function(t) mean(subset(df,ID_YR==t)$rugged)  )

layout(matrix(c(1:2),nrow=1))

ndx <- ceiling(seq(0.1,nrow(df),length=10000))
plot(ndviav~rugav,xlab="Mean Ruggedness for HR",ylab="Mean NDVI for HR")
abline(lm(ndviav~rugav),lwd=2)
ci <- predict(lm(ndviav~rugav),newdata=data.frame(rugav=seq(min(rugav),max(rugav),length=100)),interval="confidence")
# summary(lm(ndviav~rugav))    positive relationship
lines(seq(min(rugav),max(rugav),length=100),ci[,2],lty=2)
lines(seq(min(rugav),max(rugav),length=100),ci[,3],lty=2)
legend("topleft",bty="n",legend="a)",cex=1.3)

# pairs(df_rsf[ndx,predvars])  # no evidence for tradeoffs just looking at correlations
# temp <- lapply(vars,function(t) data[[t]] <<- pmin(pmax(data[[t]],-4),4) )



plot(ndviav~rugsel,xlab="Selection for Ruggedness",ylab="Mean NDVI for HR")
legend("topleft",bty="n",legend="b)",cex=1.3)
# summary(lm((ndviav~rugsel)))  #positive relationship, almost significant
# abline(lm(data$NDVI~data$Rugged_c))

# plot(data$elevation~data$NDVI_c)
# summary(lm(data$elevation~data$NDVI_c))  #no relationship

dev.off()


########### break down by repro status


graphics.off()

svg(filename="tradeoffs_v4.svg",width=6,height=7)


layout(matrix(c(1:6),nrow=3,byrow = T))
r=rep_stats[1]
par(mai=c(0.3,0.3,0.1,0.1))
for(r in rep_stats){
  df <- subset(df_rsf,USED==0&STATUS==r)
  allinds <- rownames(ranef(mod_final2)$cond$ID)
  ndx2 <- grepl(r,rownames(ranef(mod_final2)$cond$'STATUS:ID'))
  ndx3 <- grepl(sprintf("rugged:STATUS%s",r),names(fixef(mod_final2)$cond))
  thisinds <- gsub(sprintf("%s:",r),"", rownames(ranef(mod_final2)$cond$'STATUS:ID'[ndx2,]))
  ndx4 <- sapply(allinds,function(t) t%in%thisinds  )
  
  rugsel <- ranef(mod_final2)$cond$ID[ndx4,"rugged"] + ranef(mod_final2)$cond$'STATUS:ID'[ndx2,"rugged"] + fixef(mod_final2)$cond["rugged"] + ifelse(r!="prepart", fixef(mod_final2)$cond[ndx3],0)
  names(rugsel) <-  rownames(ranef(mod_final2)$cond$ID)[ndx4]
  
  ndviav <- sapply(names(rugsel),function(t) mean(subset(df,ID_YR==t)$NDVI)  )
  rugav <- sapply(names(rugsel),function(t) mean(subset(df,ID_YR==t)$rugged)  )
 
  
  ndx <- ceiling(seq(0.1,nrow(df),length=10000))
  plot(ndviav~rugav,xlab="",ylab="")   # Mean Ruggedness for HR Mean NDVI for HR
  abline(lm(ndviav~rugav),lwd=2)
  ci <- predict(lm(ndviav~rugav),newdata=data.frame(rugav=seq(min(rugav),max(rugav),length=100)),interval="confidence")
  # summary(lm(ndviav~rugav))    positive relationship
  lines(seq(min(rugav),max(rugav),length=100),ci[,2],lty=2)
  lines(seq(min(rugav),max(rugav),length=100),ci[,3],lty=2)
  #legend("topleft",bty="n",legend="a)",cex=1.3)
  
  # pairs(df_rsf[ndx,predvars])  # no evidence for tradeoffs just looking at correlations
  # temp <- lapply(vars,function(t) data[[t]] <<- pmin(pmax(data[[t]],-4),4) )
  
  
  
  plot(ndviav~rugsel,xlab="",ylab="")   # Selection for Ruggedness   Mean NDVI for HR
  #legend("topleft",bty="n",legend="b)",cex=1.3)
  # summary(lm((ndviav~rugsel)))  #positive relationship, almost significant
  # abline(lm(data$NDVI~data$Rugged_c))
  
  # plot(data$elevation~data$NDVI_c)
  # summary(lm(data$elevation~data$NDVI_c))  #no relationship

}
dev.off()





#############################
# look at functional responses... by var and repro status?

svg("funcresp4.svg",7,7)

layout(matrix(1:12,nrow=4,byrow=F))
par(mai=c(0.6,0.6,0.4,0.1))

repro <- "atheel"
var <- "NDVI"

for(repro in rep_stats){
  for(var in predvars){
    used <- subset(df_rsf,USED==1&STATUS==repro)
    avail <- subset(df_rsf,USED==0&STATUS==repro)
    IDz <- unique(used$ID_YR)
    meanz_u <- sapply(IDz,function(t) mean(used[[var]][used$ID_YR==t]))
    meanz_a <- sapply(IDz,function(t) mean(avail[[var]][used$ID_YR==t]))
    
    plot(meanz_u~meanz_a,xlab=sprintf("Mean %s, Available",var),ylab=sprintf("Mean %s, Used",var),main=repro)
    lm1 <- lm(meanz_u~meanz_a)
    abline(lm1,col="gray")
    ci <- predict(lm(meanz_u~meanz_a),newdata=data.frame(meanz_a=seq(min(meanz_a),max(meanz_a),length=100)),interval="confidence")
    # summary(lm(ndviav~rugav))    positive relationship
    lines(seq(min(meanz_a),max(meanz_a),length=100),ci[,2],lty=2,col="gray")
    lines(seq(min(meanz_a),max(meanz_a),length=100),ci[,3],lty=2,col="gray")
    
    abline(0,1,lwd=2,lty=2)
    
  }
}
dev.off()


# graphics.off()
# layout(matrix(1:3,nrow=1,byrow=T))
# par(mai=c(0.6,0.6,0.4,0.1))
# 
# repro <- "atheel"
# for(repro in rep_stats){
#   used <- subset(df_rsf,USED==1&STATUS==repro)
#   avail <- subset(df_rsf,USED==0&STATUS==repro)
#   IDz <- unique(used$ID)
#   meanz_u <- sapply(IDz,function(t) mean(used[["Rugged"]][used$ID==t]*used[["NDVI"]][used$ID==t]))
#   meanz_a <- sapply(IDz,function(t) mean(avail[["Rugged"]][used$ID==t]*avail[["NDVI"]][used$ID==t]))
#   
#   plot(meanz_u~meanz_a,xlab="Mean NDVI*rug, Available",ylab="Mean Mean NDVI*rug, Used",main=repro)
#   lm1 <- lm(meanz_u~meanz_a)
#   abline(lm1,col="gray")
#   ci <- predict(lm(meanz_u~meanz_a),newdata=data.frame(meanz_a=seq(min(meanz_a),max(meanz_a),length=100)),interval="confidence")
#   # summary(lm(ndviav~rugav))    positive relationship
#   lines(seq(min(meanz_a),max(meanz_a),length=100),ci[,2],lty=2,col="gray")
#   lines(seq(min(meanz_a),max(meanz_a),length=100),ci[,3],lty=2,col="gray")
#   
#   abline(0,1,lwd=2,lty=2)
#   
#   
#   
# }






###########
# TODO:
##########

# use the difference in resource selection for NDVI and Ruggedness from pre to at-heel and compare with survival??? [not necessary]

# figure illustrating survival probability as a function of NDVI


# look at functional responses, and include an interaction between ruggedness and NDVI?  


###################
# FAWN SURVIVAL ANALYSIS ETC
###################


#################
# Clear workspace
#################

# rm(list=ls())

#################
# Global vars
#################

ndays <- 120          # total days in the dataset
ndays2 <- 120  #120     # days to consider
ndays3 <- 30


#################
# Load packages
#################

library(lme4)
library(glmmTMB)
library(runjags)
library(loo)
library(survival)
library(lubridate)
library(buildmer)
library(logistf)

#library(ISwR)
#library(R2jags)
#library(car)
#library(tidyr)
#library(raster)
#library(tidyverse)
#library(amt)
#library(INLA)
# library(ranger)
# library(coxme)
# library(ggfortify)



#################
# Recruitment success analysis
##################

##################
# Read and process data
##################

#######
# telemetry data for females with fawns (at heel status only)

df<-read.csv("120day_atheel.csv", header=T,stringsAsFactors = F)   # "atheel_seperate_indv.csv"

head(df)
summary(df)
   #check # of deer
length(unique(df$ID_YR))       # 67 deer/year combinations   

table(df$mojaveg_ra)
table(df$AREA)


IDs <- unique(df$ID_YR)

nrow(df)    # 123,914 observations
head(df)


df$DATE2 <- ymd(df$DATE)    # convert to date format

allvars <- c(
  "Rugged",
  "dist2WATER",
  "treecvr",
  "slope",
  "sine",
  "shrubcvr",
  "elevation",
  "cosine",
  "NDVI"
)

temp <- lapply(allvars,function(t) {df[[t]] <<- scale(df[[t]])}   )  # standardize

lapply(allvars,function(t) hist(df[[t]],main=t))

#######
# recruitment success for females with fawns

df2<-read.csv("recruit_success.csv", header=T,stringsAsFactors = F)
names(df2)[1] <- "ID"
df2$Recruit_Success <- 1-df2$Recruit_Success

setdiff(df2$ID_YR,IDs)
setdiff(IDs,df2$ID_YR)
IDs <- intersect(IDs,df2$ID_YR)

head(df2)

df2$Capture.date2 <- mdy(df2$Capture.date)
df2$Birth.date2 <- mdy(df2$Birth.Date)
df2$Death.date2 <- mdy(df2$Death.Date)

head(df2,20)

nrow(df2)       # 67 rows


# #########
# # data for RSFs
# 
# avail <- read.csv("120_use_rand.csv",header=T,stringsAsFactors = F)
# 
# avail$DATE2 <- ymd(avail$DATE)
# 
# names(avail)
# 
# temp <- lapply(allvars,function(t) {avail[[t]] <<- scale(avail[[t]])}   )  # standardize
# 
# used = avail[which(avail$USED==1),] # grab the used sample
# avail = avail[which(avail$USED==0),] # grab the available sample
# 
# library(dplyr)
# firstdates <- used %>% 
#               group_by(ID_YR) %>%
#               summarize(min(DATE2))
# 
# t=1
# 
# used$dayssince <- sapply(1:nrow(used),function(t) {
#   abs(difftime(used$DATE2[t], firstdates$`min(DATE2)`[which(firstdates$ID_YR==used$ID_YR[t])],units = 'days'))
# } )
# 
# summary(used)
# 
# avail$dayssince <- NA
# 
# nrow(avail)    # 371,742 available points 
# nrow(used)     # 123,914 used points
# 
# #rm(used)
# 
# alldat <- rbind(avail,used)
# 
# names(alldat)
# 
# 
# ##############
# # Make separate GLM models for each individual
# ##############
# 
# #empty list for model storage
# mods<-list()
# mods2 <- list()
# 
# #empty list for data storage
# RSFdata<-list()
# 
# #list of IDs
# IDs
# 
# predvars <- c("Rugged","dist2WATER","NDVI","elevation","slope")
# 
# #loop and run GLM for each indv
# i = IDs[1]
# for (i in IDs){
#   thisused <- subset(used,((ID_YR==i)&(dayssince<=ndays3)))
#   thisavail0 <- subset(avail,(ID_YR==i))
#   numavail <- min(nrow(thisused)*3,nrow(thisavail0)) 
#   thisavail <- thisavail0[sample(1:nrow(thisavail0),numavail,replace=FALSE),]
#   RSFdata[[i]] <- rbind(thisused,thisavail)
#   mod<-logistf(USED~Rugged+dist2WATER+NDVI, data = RSFdata[[i]])
#   mods2[[i]] <- sapply(predvars,function(t) coef(logistf(as.formula(sprintf("USED~%s",t)),data=RSFdata[[i]])) )
#   mods[[i]]<-mod
# }
# 
# summary(mods[[1]])
# coef(mods[[1]])
# 
#    # histograms of backgrounds
# sapply(IDs,function(t) hist(RSFdata[[t]]$dist2WATER[which(RSFdata[[t]]$USED==0)],xlim=c(-3,3),main=t,xlab="dist to water (std)"))
# sapply(IDs,function(t) hist(RSFdata[[t]]$Rugged[which(RSFdata[[t]]$USED==0)],xlim=c(-3,3),main=t,xlab="ruggedness (std)"))
# sapply(IDs,function(t) hist(RSFdata[[t]]$NDVI[which(RSFdata[[t]]$USED==0)],xlim=c(-3,3),main=t,xlab="ndvi (std)"))
# 
# 
# hist(sapply(mods,function(t) coef(t)[2] ) , main="ruggedness")  # ruggedness
# hist(sapply(mods,function(t) coef(t)[3] ) , main="dist2h20")  
# hist(sapply(mods,function(t) coef(t)[4] ) , main="ndvi")  

###########
# assemble data frame for modeling (each row is a different doe/year, each column is a covariate for explaining recruitment success)
###########

data<-data.frame(ID=IDs,stringsAsFactors = F)
rownames(data) <- data$ID

head(data)
# 
# ###### 
# # Add columns to modeling data frame
# 
vars <- sprintf("%s_c",predvars) #c("Rugged", "dist2WATER", "NDVI","elevation","slope"))
# sdvars <- paste(vars,"_sd",sep="")
# sigvars <- paste(vars,"_sig",sep="")
# 
temp <- lapply(vars, function(t) data[[t]] <<- NA)
# temp <- lapply(sdvars, function(t) data[[t]] <<- NA)
# temp <- lapply(sigvars, function(t) data[[t]] <<- NA)
# 
# data
# 

#######
# Put model info (e.g., coefficients) into the modeling data frame

# i=1
# for (i in 1:nrow(data)){
#   coefs<-mods2[[data$ID[i]]][2,]
#   data[i,vars] <- coefs
#   #data[i,sdvars] <- coefs[,'Std. Error']
#   #data[i,sigvars] <- ifelse(coefs[,'Pr(>|t|)']<0.05, 1, 0)
# 
# }
# 
# data


########
# Add mean covariate value of used sites to the modeling data frame

#toSummarize <- c("Rugged","NDVI","dist2WATER","elevation","SHRUBCVR","COSINE","SINE","slope","treecover","dist2MAIN","dist24x4")

temp <- lapply(allvars,function(t) data[[t]] <<- 
                 sapply(data$ID,function(y) {
                   sss <- subset(df,ID_YR==y)
                   return(mean(sss[[t]]))
                 }))
                 

hist(df$Rugged)

head(data)

#meancovcols <- paste("mean",toSummarize,sep="")




#######
#add recruitment success to modeling data frame

data <- dplyr::left_join(data,df2,by = c("ID" = "ID_YR"))

head(data)
nrow(data)

######## Simple testing of environment effect on recruit success

t=allvars[1]
allmods <- lapply(allvars,function(t) glm(as.formula(sprintf("Recruit_Success~%s",t)),family="binomial",data=data))
names(allmods) <- allvars

lapply(allmods,summary)     # nothing significant


# allmods_c <- lapply(vars,function(t) glm(as.formula(sprintf("Recruit_Success~%s",t)),family="binomial",data=data))
# names(allmods_c) <- vars
# 
# lapply(allmods_c,summary)  # NDVI coefficient and ruggedness coefficients are significant
                           # selection for ruggedness is negatively associated with survival
                           # selection for NDVI is positively associated with survival



plot(data$Recruit_Success~data$Rugged)    # the more rugged the terrain, the lower the prob of successful recruitment?
curve(plogis(allmods$Rugged$coefficients['(Intercept)']+allmods$Rugged$coefficients['Rugged']*x),add=T)


# plot(data$Recruit_Success~data$Rugged_c)    # the more rugged the terrain, the lower the prob of successful recruitment?  One high-influence point!
# curve(plogis(allmods_c$Rugged_c$coefficients['(Intercept)']+allmods_c$Rugged$coefficients['Rugged_c']*x),add=T)


test<-glm(Recruit_Success~Rugged+NDVI+slope,data=data,family="binomial")
summary(test)

test<-glm(Recruit_Success~Rugged+NDVI+slope,data=data,family="binomial")
summary(test)

test<-glm(Recruit_Success~Rugged+NDVI+slope,data=data,family="binomial")
summary(test)


test<-glm(Recruit_Success~NDVI+slope,data=data)      # nothing useful 
summary(test)
# no significance anywhere :(

test<-glm(Recruit_Success~Rugged+NDVI+dist2WATER,data=data,family="binomial")      # nothing useful 
summary(test)

# test<-glm(Recruit_Success~Rugged_c+NDVI_c+dist2WATER_c,data=data,family="binomial")      # direction opposite from that hypothesized
# summary(test)
# 
# test<-glm(Recruit_Success~Rugged_c+NDVI_c+slope_c,data=data,family="binomial")      # direction opposite from that hypothesized
# summary(test)

# look at more fine-grained patterns. Does the risk of mortality change by habitat type?












##########
# Survival analysis
##########


data2 <- data.frame(
  ID=data$ID
)

difftime(df2$Death.date2,df2$Birth.date2,units="days")   

data2$birth.date <- do.call("c", lapply(data2$ID,function(t) min(subset(df,ID_YR==t)$DATE2))) 
data2$birth.date2 <- data$Birth.date2
data2$YEAR <- year(data2$birth.date)
data2$final.date <- do.call("c", lapply(data2$ID,function(t) max(subset(df,ID_YR==t)$DATE2)))
data2$final.date2 <- data$Death.date2
data2$Success <- data$Recruit_Success
data2$time.monitored <- difftime(data2$final.date,data2$birth.date,units="days")
data2$time.monitored2 <- difftime(data2$final.date2,data2$birth.date2,units="days")
data2$time.monitored[is.na(data2$time.monitored)] <- ndays
data2$time.monitored2[is.na(data2$time.monitored2)] <- ndays

head(data2)

# simple K-M

with(data2,Surv(time.monitored2, Success == 0))
km1 <- survfit(Surv(time.monitored2, Success == 0) ~ 1, data=data2)
plot(km1, xlab = "Time (days)", ylab="Proportion surviving", main="Desert Deer")
summary(km1)

km2 <- survfit(Surv(time.monitored2, Success == 0) ~ YEAR, data=data2)
km2
plot(km2, xlab = "Time (days)", ylab="Survival",col=c("black", "red"), lty = 1:2, main="Kaplan-Meier Survival vs. YEAR")
survdiff(Surv(time.monitored2, Success == 0) ~ YEAR, data=data2)


#### cox proportional hazards model (not really appropriate here)

coxph.year <- coxph(Surv(time.monitored2, Success == 0) ~ YEAR, data=data2)
summary(coxph.year)


data2$Rugged <- data$Rugged
coxph.rugged <- coxph(Surv(time.monitored2, Success == 0) ~ YEAR+Rugged, data=data2)
summary(coxph.rugged)    # ruggedness coef not significant


########
# Break into time intervals
########

temp <- subset(df,ID_YR==data$ID[1])
sort(temp$DATE2)

?survSplit

data3 <- survSplit(Surv(time.monitored2, Success == 0) ~ ID + YEAR + birth.date,
                     cut=seq(1,ndays,1),zero=0,data=data2)


data3$DATE <- data3$birth.date+days(data3$tstart)
data3$Rugged <- NA
data3$dist2h2o <- NA
data3$ndvi <- NA
data3$elev <- NA
data3$slop <- NA

head(data3)

i=1
for(i in 1:nrow(data3)){
  temp <- subset(df,(ID_YR==data3$ID[i])&DATE2==data3$DATE[i])
  if(nrow(temp)>0){
    data3$Rugged[i] <- mean(temp$Rugged)
    data3$dist2h2o[i] <- mean(temp$dist2WATER)
    data3$ndvi[i] <- mean(temp$NDVI) 
    data3$elev[i] <- mean(temp$elevation)
    data3$slop[i] <- mean(temp$slope)
  }
}

head(data3)
summary(data3)

data3$ID2 <- as.numeric(as.factor(data3$ID))
data3 <- data3[order(data3$ID2),]
data3$ID2

data3$YEAR2 <- as.numeric(as.factor(data3$YEAR))


#### capture history

alive <- matrix(NA,nrow=nrow(data2),ncol=ndays2+1)
nobs <- numeric(nrow(data2))

i=1
for(i in 1:nrow(data2)){
  temp <- subset(data3[data3$tstart<=ndays2,],ID2==i)
  alive[i,1:nrow(temp)] <- 1-temp$event
  nobs[i] <- nrow(temp)
}


#### covars

rugged <- matrix(NA,nrow=nrow(data2),ncol=ndays2+1)
dist2h2o <- matrix(NA,nrow=nrow(data2),ncol=ndays2+1)
ndvi <- matrix(NA,nrow=nrow(data2),ncol=ndays2+1)
elev <- matrix(NA,nrow=nrow(data2),ncol=ndays2+1)
slop <- matrix(NA,nrow=nrow(data2),ncol=ndays2+1)

i=1
for(i in 1:nrow(data2)){
  temp <- subset(data3[data3$tstart<=ndays2,],ID2==i)
  rugged[i,1:nrow(temp)] <- temp$Rugged
  dist2h2o[i,1:nrow(temp)] <- temp$dist2h2o
  ndvi[i,1:nrow(temp)] <- temp$ndvi
  elev[i,1:nrow(temp)] <- temp$elev
  slop[i,1:nrow(temp)] <- temp$slop
}

######## try basic logistic regression model with time intervals?  Write out in JAGS 


cat("


model{

  ########## interpolation

  for(ind in 1:ninds){
    meanrugged[ind] ~ dnorm(0,0.1)
    precrugged[ind] ~ dgamma(0.1,0.1)
    meanelev[ind] ~ dnorm(0,0.1)
    precelev[ind] ~ dgamma(0.1,0.1)
    # meanslop[ind] ~ dnorm(0,0.1)
    # precslop[ind] ~ dgamma(0.1,0.1)
    meandist[ind] ~ dnorm(0,0.1)
    precdist[ind] ~ dgamma(0.1,0.1)
    meanndvi[ind] ~ dnorm(0,0.1)
    precndvi[ind] ~ dgamma(0.1,0.1)
    for(t in 1:ntimes[ind]){
      rugged[ind,t] ~ dnorm(meanrugged[ind],precrugged[ind])
      dist2h2o[ind,t] ~ dnorm(meandist[ind],precdist[ind])
      ndvi[ind,t] ~ dnorm(meanndvi[ind],precndvi[ind])
      elev[ind,t] ~ dnorm(meanelev[ind],precelev[ind])
      #slop[ind,t] ~ dnorm(meanslop[ind],precslop[ind])
    }
  }

  ########## SURVIVAL

  yearprec ~ dgamma(0.1,0.1)
  yearsd <- pow(1/yearprec,0.5)
  for(y in 1:nyears){
    yeareff[y] ~ dnorm(0,yearprec)
  }

  for(ind in 1:ninds){
    for(t in 1:(ntimes[ind]-1)){
      logit(psurviv[ind,t]) <- basesurv.l + yeareff[year[ind]] + ruggedeff * rugged[ind,t] +    # + ideff[id[ind]]
                                     disteff * dist2h2o[ind,t] + ndvieff * ndvi[ind,t] +
                                     eleveff * elev[ind,t] #+ slopeff * slop[ind,t]
    }
  }

  ########## LIKELIHOOD

  for(ind in 1:ninds){
    alive[ind,1] ~ dbern(1)
    for(t in 2:ntimes[ind]){
      palive[ind,t] <- alive[ind,t-1] * psurviv[ind,t-1]
      alive[ind,t] ~ dbern(palive[ind,t])
    }
  }

  ########## PRIORS

    basesurv ~ dunif(0,1)
    basesurv.l <- log(basesurv/(1-basesurv))
    ruggedeff ~ dnorm(0,0.1)
    disteff ~ dnorm(0,0.1)
    ndvieff ~ dnorm(0,0.1)
    eleveff ~ dnorm(0,0.1)
    #slopeff ~ dnorm(0,0.1)
    
  ######### PREDICTIONS
    # NDVI
    for(g in 1:npts){
      logit(psurviv_ndvi[g]) <- basesurv.l + yeareff[2] + ndvieff * ndvis[g]    # + ideff[id[ind]]
    }
    
    # Ruggedness
    for(g in 1:npts){
      logit(psurviv_rugged[g]) <- basesurv.l + yeareff[2] + ruggedeff * ruggeds[g]    # + ideff[id[ind]]
    }
    
    # Elevation
    for(g in 1:npts){
      logit(psurviv_elev[g]) <- basesurv.l + yeareff[2] + eleveff * elevs[g]    # + ideff[id[ind]]
    }
    
    # DistH2O
    for(g in 1:npts){
      logit(psurviv_dist[g]) <- basesurv.l + yeareff[2] + disteff * dist2h2os[g]    # + ideff[id[ind]]
    }
}
  
    " 
    ,file="mojavedeer_jags.txt"
)

npts = 20
ndvis <- seq(min(data3$ndvi,na.rm = T),max(data3$ndvi,na.rm = T),length=npts)
elevs <- seq(min(data3$elev,na.rm = T),max(data3$elev,na.rm = T),length=npts)
dist2h2os <- seq(min(data3$dist2h2o,na.rm = T),max(data3$dist2h2o,na.rm = T),length=npts)
ruggeds <- seq(min(data3$Rugged,na.rm = T),max(data3$Rugged,na.rm = T),length=npts)

data.for.jags <- list(
  alive = alive,
  year = tapply(data3$YEAR2,data3$ID2,function(t) t[1]),
  nyears = max(data3$YEAR2),
  ninds = nrow(data2),
  ntimes = nobs,
  rugged = rugged,
  dist2h2o = dist2h2o,
  ndvi = ndvi,
  elev = elev,
  #slop = slop,
  npts = npts,
  ndvis = ndvis,
  elevs = elevs,
  ruggeds = ruggeds,
  dist2h2os = dist2h2os
)

data.for.jags$ninds

initsfunc <- function(){
  list(
    basesurv = runif(1,0.9,0.95),
    ruggedeff = runif(1,-0.01,0.01),
    disteff = runif(1,-0.01,0.01),
    ndvieff = runif(1,-0.01,0.01),
    eleveff = runif(1,-0.01,0.01),
    #slopeff = runif(1,-0.01,0.01),
    yearprec = runif(1,9,10)
  )
}

initsfunc()

paramstosave <- c(
  "basesurv",
  "ruggedeff",
  "disteff",
  "ndvieff",
  "eleveff",
  #"slopeff",
  "yeareff",
  "meanrugged",
  "psurviv_rugged",
  "psurviv_ndvi",
  "psurviv_elev",
  "psurviv_dist"
)

  #?run.jags
mod <- run.jags("mojavedeer_jags.txt",data=data.for.jags, monitor=paramstosave,inits=initsfunc,n.chains = 3,
          burnin=5000,sample=10000,adapt=1000,method = 'parallel',summarise = F
         )

save(mod,file="JAGS_final.RData")
# failed.jags(c('model'))
# failed.jags(c('data'))
# failed.jags(c('inits'))

#?autorun.JAGS

library(coda)
mod2 <- as.mcmc.list(mod)

#varndx = 4
#traceplot(mod2)    #main=paramstosave[varndx]

#str(mod2)
hist(mod2[[1]][,"basesurv"])
hist(mod2[[1]][,"ruggedeff"])
hist(mod2[[1]][,"disteff"])
hist(mod2[[1]][,"ndvieff"],main="NDVI effect on survival",xlab="Coefficient")
hist(mod2[[1]][,"eleveff"],main="Elevation effect on survival",xlab="Coefficient")
hist(mod2[[1]][,"slopeff"],main="Slope effect on survival",xlab="Coefficient")

plot(mod2[,"basesurv"])
plot(mod2[,"eleveff"])
plot(mod2[,"ruggedeff"])
plot(mod2[,"ndvieff"])
plot(mod2[,"disteff"])

mod3 <- combine.mcmc(mod2)
mod3 <- as.matrix(mod3)
totsamps <- nrow(mod3) 
allnames <- colnames(mod3)



##### make survival figure (partial dependence plots)


df4 <- read.csv("120day_atheel.csv", header=T,stringsAsFactors = F) # recover original data

npts
nsamps <- 1000

svg("surv_fig.svg",6,5)

layout(matrix(1:4,nrow=2,byrow = T))
par(mai=c(0.7,0.7,0.1,0.1))


vars <- c("ndvi","elev","rugged","dist")
vars2 <- c("NDVI","elevation","Rugged","dist2WATER")
vars3 <- c("NDVI","Elevation","Ruggedness","Distance to Water")
vars4 <- c("ndvi","elev","Rugged","dist2h2o")

rnges <- list(ndvis,elevs,ruggeds,dist2h2os) 
rnds <- c(2,0,4,0)

var="elev"

panel=1
for(panel in 1:4){
  var = vars[panel]
  var2 = vars2[panel]
  var3 = vars3[panel]
  var4 = vars4[panel]
  thisrng = rnges[[panel]]
  
  
  plot(0.5,0.5,pch="",xlab=var3,xaxt="n",ylab="Survival rate (daily)",ylim=c(0.9,1),xlim=c(min(thisrng),max(thisrng)))
  
  samps <- sample(1:totsamps,nsamps,replace = F)
  p=1
  this <- matrix(0,nrow=nsamps,ncol=npts)
  for(p in 1:npts){
    this[,p] <- mod3[samps,sprintf("psurviv_%s[%s]",var,p)]
  }
  medians <- apply(this,2,median)
  maxs <- apply(this,2,quantile,probs=0.95)
  mins <- apply(this,2,quantile,probs=0.05)
  
  Hmisc::errbar(thisrng,medians,maxs,mins,add=T)
  axis(1,at=seq(min(thisrng),max(thisrng),length=5),labels=round(mean(df4[[var2]],na.rm=T)+sd(df4[[var2]],na.rm=T)*seq(min(thisrng),max(thisrng),length=5), rnds[panel]))
  rug(data3[[var4]][sample(1:nrow(data3),200)])
}
dev.off()


   # extract coefficients and credible intervals

library(coda)

mod4 <- combine.mcmc(mod2)

mean(mod3[,"ndvieff"])
HPDinterval(mod4[,"ndvieff"])

mean(mod3[,"eleveff"])
HPDinterval(mod4[,"eleveff"])

mean(mod3[,"disteff"])
HPDinterval(mod4[,"disteff"])

mean(mod3[,"ruggedeff"])
HPDinterval(mod4[,"ruggedeff"])


















###################
# SIDE ANALYSIS: functional response!
# 


# #Visualizations..
# #######functional response testing Holbrook et al.(2019) 
# library("plyr")
# 
# ?tdc
# 
# ####Method 1 plotting
# Used_Rug = ddply(used, .(ID), summarize, mean = mean(Rugged)) # calculate the mean value by individual 
# Avail_Rug = ddply(avail, .(ID), summarize, mean = mean(Rugged)) # calculate the mean value by individual 
# par(pty="s")
# plot(Used_Rug$mean~Avail_Rug$mean, xlim=c(-1,3),ylim=c(-1,3))
# abline(0,1)
# 
# Used_Wtr = ddply(used, .(ID), summarize, mean = mean(dist2WATER)) # calculate the mean value by individual 
# Avail_Wtr = ddply(avail, .(ID), summarize, mean = mean(dist2WATER)) # calculate the mean value by individual 
# par(pty="s")
# plot(Used_Wtr$mean~Avail_Wtr$mean, xlim=c(-1,3),ylim=c(-1,3))
# abline(0,1)
# 
# Used_ndvi = ddply(used, .(ID), summarize, mean = mean(NDVI)) # calculate the mean value by individual 
# Avail_ndvi = ddply(avail, .(ID), summarize, mean = mean(NDVI)) # calculate the mean value by individual 
# par(pty="s")
# plot(Used_ndvi$mean~Avail_ndvi$mean, xlim=c(-1,3),ylim=c(-1,3))
# abline(0,1)
# 
# ###Method 3A plotting
# plot(rugged$Estimate~Avail_Rug$mean, xlim=c(-1,2),ylim=c(-1,2))
# abline(0,1)
# abline(0,0, col="red")
# 
# plot(water$Estimate~Avail_Wtr$mean, xlim=c(-1,2),ylim=c(-1,2))
# abline(0,1)
# abline(0,0, col="red")
# 
# plot(ndvi$Estimate~Avail_ndvi$mean, xlim=c(-1,2),ylim=c(-1,2))
# abline(0,1)
# abline(0,0, col="red")
# 
# 


# 
# #############
# # try evaluating whether resource selection patterns change according to successful and unsuccessful?
# 
# # note: selection strength is proportional to the fraction of background points between the 5 and 95% quantiles of preferred habitat  
# 
# cat("
# 
# model{
# 
# #####priors
# 
# for(v in 1:nvars){
#   mean_wanted[v] ~ dnorm(0,0.01)I(-4,4)
#   var_wanted[v] ~ dgamma(1,1)     # previously set to 1
#   prec_wanted[v] <- 1/var_wanted[v] 
#   sd_wanted[v] <- pow(var_wanted[v],0.5)
# }
# 
# 
# meanp ~ dunif(0,1)    # assume selection propensity at global mean covar values does not vary across individuals
# logit.meanp <- log(meanp/(1-meanp))      # use propensity at global mean covariate values
# 
# 
# for(i in 1:ninds){
#   for(v in 1:nvars){
#     for(back in 1:nbacks){
#       prob.back[i,v,back] <- 1-pnorm(background[i,v,back],mean_wanted[v],prec_wanted[v])
#       q.back[i,v,back] <- log(prob.back[i,v,back]/(1-prob.back[i,v,back]))
#     }
#     coef[i,v] <- mean(q.back[i,v,1:nbacks])
#   }
# }
# 
# for(i in 1:ninds){
#   for(obs in 1:nobs[i]){
#     logit(pused[i,obs]) <- logit.meanp + coef[i,1]*ruggedness[i,obs] + coef[i,2]*ndvi[i,obs] + 
#                               coef[i,3]*disth20[i,obs]
#                               
#     used[i,obs] ~ dbern(pused[i,obs])
#     lik[i,obs] <- dbin(used[i,obs],pused[i,obs],1)
#     loglik[i,obs] <- log(lik[i,obs])
#   }
#   
# }
# 
# 
# }
#     "
# ,file="HSjags1_1.txt"
# )
# 
# #?loo::waic
# 
# 
# # #########
# # # playspace...
# # 
# backtest <- rnorm(100,0,1)
# wanted <- c(5,1)
# curve(dnorm(x,wanted[1],wanted[2]),-4,4,ylim=c(0,1))
# hist(backtest,add=t,lty="blank",col=gray(0.5),freq = F)
# qs <- 0.5-pnorm(backtest,wanted[1],wanted[2])
# mean(qs)
# 
# backtest <- rnorm(100,0,1)   # standard normal
# wanted <- c(-2,5)
# curve(dnorm(x,wanted[1],wanted[2]),-4,4,ylim=c(0,1))
# hist(backtest,add=t,lty="blank",col=gray(0.5),freq = F)
# qs <- 1-pnorm(backtest,wanted[1],wanted[2])
# log(mean(qs)/(1-mean(qs)))
# mean(log(qs/(1-qs)))
# 
# hist(data$Rugged_c)
# hist(data$dist2WATER_c)
# hist(data$NDVI_c)
# 
# a <- lapply(RSFdata2,function(t) glm(USED ~ NDVI,data=t))
# lapply(a,summary)
# allcoefs <- sapply(a,function(t) coef(t)[2])
# meanbacks <- sapply(RSFdata2,function(t) mean(t$NDVI[t$USED==0,]))
# plot(allcoefs~meanbacks)   
# summary(lm(allcoefs~meanbacks))    # weak negative relationship but not significant
# #sapply(1:length(a),function(t) hist(background[t,2,],xlim=c(-3,3),freq=F,col=rgb(t/250,t/250,t/250,0.2),breaks=10,main=allcoefs[[t]],lty="blank"))
# 
# a <- lapply(RSFdata2,function(t) glm(Rugged ~ NDVI,data=t))
# #lapply(a,summary)
# allcoefs <- sapply(a,function(t) coef(t)[2])
# meanbacks <- sapply(RSFdata2,function(t) mean(t$Rugged[t$USED==0,]))
# plot(allcoefs~meanbacks)   
# summary(lm(allcoefs~meanbacks))    # positive relationship- the more rugged habitat available the more they select for it!
# sapply(1:length(a),function(t) hist(background[t,1,],xlim=c(-3,3),freq=F,col=rgb(t/250,t/250,t/250,0.2),breaks=10,main=allcoefs[[t]],lty="blank"))
# # where available habitat is not rugged, they seem to avoid the rugged areas!
# 
# str(background)
# 
# #########
# # set up data for jags
# 
# 
# sapply(RSFdata,nrow)
# 
# RSFdata2 <- list()
# i=IDs[1]
# for(i in IDs){
#   tempused <- subset(RSFdata[[i]],USED==1)
#   tempback <- subset(RSFdata[[i]],USED==0)
#   dayz <- unique(yday(tempused$DATE2))
#   tempused2 <- tempused[unlist(lapply(dayz,function(t) sample(which(yday(tempused$DATE2)==t),2,replace=FALSE))),]
#   tempback2 <- tempback[sample(1:nrow(tempback),nrow(tempused2)*2,replace=F),]
#   RSFdata2[[i]] <- rbind(tempused2,tempback2)
#   #nrow(tempused2)
# }
# 
# RSFdata2[[60]]$ID_YR
# 
# ninds <- length(RSFdata2)
# maxobs <- max(sapply(RSFdata2,nrow))
# nobs <- sapply(RSFdata2,nrow)
# nvars <- 3
# maxback <- 250  #max(sapply(RSFdata2,function(t) nrow(subset(t,USED==0)) ))
# 
# nbacks <- 250 # sapply(RSFdata,function(t) nrow(subset(t,USED==0)) )
# 
# varnames <- c("Rugged","NDVI","dist2WATER")
# 
# used <- array(NA,dim=c(ninds,maxobs))
# ruggedness <- array(NA,dim=c(ninds,maxobs))
# ndvi <- array(NA,dim=c(ninds,maxobs))
# disth20 <- array(NA,dim=c(ninds,maxobs))
# succeeded <- numeric(ninds)
# 
# i=24
# for(i in 1:ninds){
#   used[i,1:nrow(RSFdata2[[i]])] <- RSFdata2[[i]]$USED
#   ruggedness[i,1:nrow(RSFdata2[[i]])] <- RSFdata2[[i]]$Rugged
#   ndvi[i,1:nrow(RSFdata2[[i]])] <- RSFdata2[[i]]$NDVI
#   disth20[i,1:nrow(RSFdata2[[i]])] <- RSFdata2[[i]]$dist2WATER
#   ndx <- which(data$ID==RSFdata2[[i]]$ID_YR[1])
#   succeeded[i] <- data$Recruit_Success[ndx]+1
# }
# used[3,]
# ruggedness[3,]
# 
# succeeded[24] <- 1   # temporary
# 
# background <- array(NA,dim=c(ninds,nvars,maxback))
# 
# i=1
# for(i in 1:ninds){
#   thisback <- subset(RSFdata[[i]],USED==0)[sample(1:nrow(subset(RSFdata[[i]],USED==0)),maxback,replace = T),]
#   v=1
#   for(v in 1:nvars){
#     background[i,v,1:nrow(thisback)] <- thisback[[varnames[v]]]
#   }
# }
# background[24,3,250]
# str(background)
# 
# data.for.jags <- list(
#   used = used,
#   background=background,
#   nvars=nvars,
#   nbacks=nbacks,
#   nobs=nobs,
#   ninds = ninds,
#   ruggedness = ruggedness,
#   ndvi = ndvi,
#   disth20 = disth20
# )
# 
# #data.for.jags
# 
# initsfunc <- function(){
#   list(
#     mean_wanted = runif(nvars,-0.1,0.1),
#     var_wanted = runif(nvars,1,2),
#     meanp = runif(1,0.4,0.5)
#   )
# }
# 
# initsfunc()
# 
# paramstosave <- c(
#   "mean_wanted",
#   "sd_wanted",
#   "meanp",
#   "coef",
#   "loglik"
# )
# 
# #?run.jags (takes a long time)
# rsfmod <- run.jags("HSjags1_1.txt",data=data.for.jags, monitor=paramstosave,inits=initsfunc,n.chains = 2,
#                    burnin=500,sample=500,adapt=500
# )
# 
# # failed.jags(c('model'))
# # failed.jags(c('data'))
# # failed.jags(c('inits'))
# 
# #?autorun.JAGS
# 
# rsfmod2 <- coda::as.mcmc.list(rsfmod)
# 
# 
# plot(rsfmod2[,'mean_wanted[2]'])
# plot(rsfmod2[,'sd_wanted[2]'])
# plot(rsfmod2[,'meanp'])
# 
# 
# curve(dnorm(x,rsfmod2[[1]][,'mean_wanted[1]'][1],rsfmod2[[1]][,'sd_wanted[1]'][1]),-4,4,ylim=c(0,1),xlab=varnames[1],ylab="Optimality")   #rsfmod2[[1]][,'sd_wanted[1]'][1]
# sapply(1:200,function(t) curve(dnorm(x,rsfmod2[[1]][,'mean_wanted[1]'][t],rsfmod2[[t]][,'sd_wanted[1]'][1]),-4,4,add=T) ) #rsfmod2[[1]][,'sd_wanted[1]'][t]
# ndxs <- c(which.min(apply(background[,1,],1,mean)),which.max(apply(background[,1,],1,mean)))
# sapply(ndxs,function(t) hist(background[t,1,],freq=F,col=rgb(t/250,t/250,t/250,0.2),add=T,breaks=10,lty="blank"))
# 
# 
# curve(dnorm(x,rsfmod2[[1]][,'mean_wanted[2]'][1],rsfmod2[[1]][,'sd_wanted[2]'][1]),-4,4,ylim=c(0,1),xlab=varnames[2],ylab="Optimality")   #rsfmod2[[1]][,'sd_wanted[1]'][1]
# sapply(1:200,function(t) curve(dnorm(x,rsfmod2[[1]][,'mean_wanted[2]'][t],rsfmod2[[1]][,'sd_wanted[2]'][t]),-4,4,add=T) ) #rsfmod2[[1]][,'sd_wanted[1]'][t]
# ndxs <- c(which.min(apply(background[,2,],1,mean)),which.max(apply(background[,2,],1,mean)))
# sapply(ndxs,function(t) hist(background[t,2,],freq=F,col=rgb(t/100,t/100,t/100,0.5),add=T,breaks=10,lty="blank"))
# 
# curve(dnorm(x,rsfmod2[[1]][,'mean_wanted[3]'][1],rsfmod2[[1]][,'sd_wanted[3]'][1]),-4,4,ylim=c(0,1),xlab=varnames[3],ylab="Optimality")   #rsfmod2[[1]][,'sd_wanted[1]'][1]
# sapply(1:200,function(t) curve(dnorm(x,rsfmod2[[1]][,'mean_wanted[3]'][t],rsfmod2[[1]][,'sd_wanted[3]'][t]),-4,4,add=T) ) #rsfmod2[[1]][,'sd_wanted[1]'][t]
# ndxs <- c(which.min(apply(background[,3,],1,mean)),which.max(apply(background[,3,],1,mean)))
# sapply(ndxs,function(t) hist(background[t,2,],freq=F,col=rgb(t/250,t/250,t/250,0.5),add=T,breaks=10,lty="blank"))
# 
# 
# # Q: is the sign of the logistic regression relations related to where the background is in relation to the desired range?
# 
# ## note: not too much difference in background percentages among individuals... so this method isn't necessarily that necessary here... 
# 
# loglikmat <- t(as.matrix(rsfmod2[,grep("loglik",colnames(rsfmod2[[1]]))]))
# #str(loglikmat)
# waic(loglikmat)
# 
# str(rsfmod2[,'loglik'])
# 
# 
# 
# ################
# # does desired habitat differ between those who were or were not successful in recruiting?
# 
# 
# cat("
# 
# model{
# 
# #####priors
# 
# #prec_wanted <- 1   #~ dgamma(0.1,0.1)
# 
# for(v in 1:nvars){
#   for(r in 1:2){
#     mean_wanted[v,r] ~ dnorm(0,0.01)I(-4,4)
#     var_wanted[v,r] ~ dgamma(1,1)     # previously set to 1
#     prec_wanted[v,r] <- 1/var_wanted[v,r] 
#     sd_wanted[v,r] <- pow(var_wanted[v,r],0.5)
#   }
# }
# 
# # p.succeeded[1] ~ dunif(0,1)
# # p.succeeded[2] <- 1-p.succeeded[1]
# # 
# # for(i in 1:ninds){
# #   succeeded[i] ~ dcat(p.succeeded[])
# # }
# 
# meanp ~ dunif(0,1)
# logit.meanp <- log(meanp/(1-meanp))      # use propensity at global mean covariate values
# 
# for(i in 1:ninds){
#   for(v in 1:nvars){
#     for(back in 1:nbacks){
#       prob.back[i,v,back] <- 1-pnorm(background[i,v,back],mean_wanted[v,succeeded[i]],prec_wanted[v,succeeded[i]])
#       q.back[i,v,back] <- log(prob.back[i,v,back]/(1-prob.back[i,v,back]))
#     }
#     coef[i,v] <- mean(q.back[i,v,1:nbacks])
#   }
# }
# 
# for(i in 1:ninds){
#   for(obs in 1:nobs[i]){
#     logit(pused[i,obs]) <- logit.meanp + coef[i,1]*ruggedness[i,obs] + coef[i,2]*ndvi[i,obs] + 
#                               coef[i,3]*disth20[i,obs]
#                               
#     used[i,obs] ~ dbern(pused[i,obs])
#     lik[i,obs] <- dbin(used[i,obs],pused[i,obs],1)
#     loglik[i,obs] <- log(lik[i,obs])
#   }
#   
# }
# 
# 
# }
#     "
# ,file="HSjags2_2.txt"
# )
# 
# 
# data.for.jags <- list(
#   used = used,
#   background=background,
#   nvars=nvars,
#   nbacks=nbacks,
#   nobs=nobs,
#   ninds = ninds,
#   ruggedness = ruggedness,
#   succeeded = succeeded,
#   ndvi = ndvi,
#   disth20 = disth20
# )
# 
# #data.for.jags
# 
# initsfunc <- function(){
#   list(
#     mean_wanted = matrix(runif(nvars*2,-0.1,0.1),ncol=2),
#     var_wanted = matrix(runif(nvars*2,2,3),ncol=2),
#     meanp = runif(1,0.4,0.5)
#   )
# }
# 
# initsfunc()
# 
# paramstosave <- c(
#   "mean_wanted",
#   "sd_wanted",
#   "meanp",
#   "coef",
#   "loglik"
# )
# 
# #?run.jags (takes a long time)
# rsfmodr <- run.jags("HSjags2_2.txt",data=data.for.jags, monitor=paramstosave,inits=initsfunc,n.chains = 2,
#                     burnin=500,sample=500,adapt=500
# )
# 
# 
# 
# rsfmodr2 <- coda::as.mcmc.list(rsfmodr)
# 
# 
# plot(rsfmodr2[,'mean_wanted[1,1]'])
# plot(rsfmodr2[,'mean_wanted[1,2]'])
# plot(rsfmodr2[,'sd_wanted[1,1]'])
# plot(rsfmodr2[,'sd_wanted[1,2]'])
# 
# plot(rsfmodr2[,'mean_wanted[2,1]'])
# plot(rsfmodr2[,'mean_wanted[2,2]'])
# plot(rsfmodr2[,'sd_wanted[2,1]'])
# plot(rsfmodr2[,'sd_wanted[2,2]'])
# 
# plot(rsfmodr2[,'mean_wanted[3,1]'])
# plot(rsfmodr2[,'mean_wanted[3,2]'])
# plot(rsfmodr2[,'sd_wanted[3,1]'])
# plot(rsfmodr2[,'sd_wanted[3,2]'])
# 
# 
# layout(matrix(c(1:3),nrow=3))
# par(mai=c(0.8,0.8,0.1,0.1))
# 
# curve(dnorm(x,rsfmodr2[[1]][,'mean_wanted[1,1]'][1],rsfmodr2[[1]][,'sd_wanted[1,1]'][1]),-4,4,ylim=c(0,1),main="",col="red",ylab="Optimality",xlab=varnames[1])   #rsfmod2[[1]][,'sd_wanted[1]'][1]
# sapply(1:200,function(t) curve(dnorm(x,rsfmodr2[[1]][,'mean_wanted[1,1]'][t],rsfmodr2[[1]][,'sd_wanted[1,1]'][t]),-4,4,add=T,col="red") ) #rsfmod2[[1]][,'sd_wanted[1]'][t]
# curve(dnorm(x,rsfmodr2[[1]][,'mean_wanted[1,2]'][1],rsfmodr2[[1]][,'sd_wanted[1,2]'][1]),-4,4,ylim=c(0,1),main=varnames[1],col="green",add=T)   #rsfmod2[[1]][,'sd_wanted[1]'][1]
# sapply(1:200,function(t) curve(dnorm(x,rsfmodr2[[1]][,'mean_wanted[1,2]'][t],rsfmodr2[[1]][,'sd_wanted[1,2]'][t]),-4,4,add=T,col="green") ) #rsfmod2[[1]][,'sd_wanted[1]'][t]
# legend("topleft",bty="n",col=c("red","green"),lty=c(1,1),legend=c("Unsuccessful","Successful"))
# ndxs <- data$Recruit_Success==1
# hist(background[ndxs,1,],freq=F,col=rgb(1,0,0,0.2),add=T,breaks=15,lty="blank")
# hist(background[!ndxs,1,],freq=F,col=rgb(0,1,0,0.2),add=T,breaks=15,lty="blank")
# 
# 
# curve(dnorm(x,rsfmodr2[[1]][,'mean_wanted[2,1]'][1],rsfmodr2[[1]][,'sd_wanted[2,1]'][1]),-4,4,ylim=c(0,1),main="",col="red",ylab="Optimality",xlab=varnames[2])   #rsfmod2[[1]][,'sd_wanted[1]'][1]
# sapply(1:200,function(t) curve(dnorm(x,rsfmodr2[[1]][,'mean_wanted[2,1]'][t],rsfmodr2[[1]][,'sd_wanted[2,1]'][t]),-4,4,add=T,col="red") ) #rsfmod2[[1]][,'sd_wanted[1]'][t]
# curve(dnorm(x,rsfmodr2[[1]][,'mean_wanted[2,2]'][1],rsfmodr2[[1]][,'sd_wanted[2,2]'][1]),-4,4,ylim=c(0,1),main=varnames[2],col="green",add=T)   #rsfmod2[[1]][,'sd_wanted[1]'][1]
# sapply(1:200,function(t) curve(dnorm(x,rsfmodr2[[1]][,'mean_wanted[2,2]'][t],rsfmodr2[[1]][,'sd_wanted[2,2]'][t]),-4,4,add=T,col="green") ) #rsfmod2[[1]][,'sd_wanted[1]'][t]
# legend("topleft",bty="n",col=c("red","green"),lty=c(1,1),legend=c("Unsuccessful","Successful"))
# ndxs <- data$Recruit_Success==1
# hist(background[ndxs,2,],freq=F,col=rgb(1,0,0,0.2),add=T,breaks=15,lty="blank")
# hist(background[!ndxs,2,],freq=F,col=rgb(0,1,0,0.2),add=T,breaks=15,lty="blank")
# 
# 
# curve(dnorm(x,rsfmodr2[[1]][,'mean_wanted[3,1]'][1],rsfmodr2[[1]][,'sd_wanted[3,1]'][1]),-4,4,ylim=c(0,1),main="",col="red",ylab="Optimality",xlab=varnames[3])   #rsfmod2[[1]][,'sd_wanted[1]'][1]
# sapply(1:200,function(t) curve(dnorm(x,rsfmodr2[[1]][,'mean_wanted[3,1]'][t],rsfmodr2[[1]][,'sd_wanted[3,1]'][t]),-4,4,add=T,col="red") ) #rsfmod2[[1]][,'sd_wanted[1]'][t]
# curve(dnorm(x,rsfmodr2[[1]][,'mean_wanted[3,2]'][1],rsfmodr2[[1]][,'sd_wanted[3,2]'][1]),-4,4,ylim=c(0,1),main=varnames[3],col="green",add=T)   #rsfmod2[[1]][,'sd_wanted[1]'][1]
# sapply(1:200,function(t) curve(dnorm(x,rsfmodr2[[1]][,'mean_wanted[3,2]'][t],rsfmodr2[[1]][,'sd_wanted[3,2]'][t]),-4,4,add=T,col="green") ) #rsfmod2[[1]][,'sd_wanted[1]'][t]
# legend("topleft",bty="n",col=c("red","green"),lty=c(1,1),legend=c("Unsuccessful","Successful"))
# ndxs <- data$Recruit_Success==1
# hist(background[ndxs,3,],freq=F,col=rgb(1,0,0,0.2),add=T,breaks=15,lty="blank")
# hist(background[!ndxs,3,],freq=F,col=rgb(0,1,0,0.2),add=T,breaks=15,lty="blank")
# 
# 
# # Q: is the sign of the logistic regression relations related to where the background is in relation to the desired range?
# 
# ## note: not too much difference in background percentages among individuals... so this method isn't necessarily that necessary here... 
# 
# loglikmatr <- t(as.matrix(rsfmodr2[,grep("loglik",colnames(rsfmodr2[[1]]))]))
# #str(loglikmat)
# waic(loglikmatr)
# 
# str(rsfmod2[,'loglik'])
# 


##############
# OLD CODE
##############

#####################
# Make datasets for at heel (30 days or until death of lamb) vs post loss
#####################

# rownames(alldat) = 1:nrow(alldat)
# atheel = subset(alldat,USED==1)
# atheel = atheel[rep(FALSE,times=nrow(atheel)),]
# lostlamb = atheel[rep(FALSE,times=nrow(atheel)),]
# 
# 
# id=IDs[5]
# for(id in IDs){
#   deathdate = df2$Death.date2[df2$ID_YR==id]
#   if(!is.na(deathdate)){
#     temp = subset(alldat,ID_YR==id&DATE2<deathdate&USED==1)
#     toadd <- rownames(alldat)%in%rownames(temp)
#     atheel = rbind(atheel,alldat[toadd,])
#     temp = subset(alldat,ID_YR==id&DATE2>=deathdate&USED==1)
#     toadd <- rownames(alldat)%in%rownames(temp)
#     lostlamb = rbind(lostlamb,alldat[toadd,])
#   }else{
#     temp = subset(alldat,ID_YR==id&USED==1)
#     toadd <- rownames(alldat)%in%rownames(temp)
#     atheel = rbind(atheel,alldat[toadd,])
#   } 
# }
# atheel = subset(atheel,dayssince<=30)
# atheel = subset(atheel,ID_YR%in%IDs)
# lostlamb = subset(lostlamb,dayssince<=90)
# lostlamb = subset(lostlamb,ID_YR%in%IDs)
# 
# nrow(atheel)
# nrow(lostlamb)
# 
# hasstages <- list()
# 
# atheel$hasstages = 0
# lostlamb$hasstages = 0
# alldat$hasstages = 0
# 
# id=IDs[1]
# for(id in IDs){
#   temp1 = subset(atheel,ID_YR==id&USED==1)
#   temp2 = subset(alldat,ID_YR==id&USED==0)
#   temp3 = temp2[sample(1:nrow(temp2),min(nrow(temp1)*5,nrow(temp2)),replace = F),]
#   atheel <- rbind(atheel,temp3)
#   
#   temp4 = subset(lostlamb,ID_YR==id&USED==1)
#   temp5 = subset(alldat,ID_YR==id&USED==0)
#   temp6 = temp2[sample(1:nrow(temp5),min(nrow(temp4)*5,nrow(temp5)),replace = F),]
#   lostlamb <- rbind(lostlamb,temp6)
#   
#   hasstages[[id]] = as.numeric(c(nrow(temp1)>0 , nrow(temp4)>0))
#   hasstages[[id]] = c(hasstages[[id]],prod(hasstages[[id]][1:2]))
#   
#   ndx = atheel$ID_YR==id
#   atheel$hasstages[ndx] <- hasstages[[id]][3]
#   ndx = lostlamb$ID_YR==id
#   lostlamb$hasstages[ndx] <- hasstages[[id]][3]
# }
# 
# nrow(atheel)
# nrow(lostlamb)
# nrow(alldat)
# 
# atheel$STATUS = "atheel"
# lostlamb$STATUS = "lostlamb"
# 
# alldat2 <- rbind(atheel,lostlamb)
# nrow(alldat2)
# 
# atheel$weights <- ifelse(atheel$USED == 1, 1, 1000)
# lostlamb$weights <- ifelse(lostlamb$USED == 1, 1, 1000)
# alldat2$weights <- ifelse(alldat2$USED == 1, 1, 1000)




