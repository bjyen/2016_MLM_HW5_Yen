#MLM HW5
# April 16th, 2016
library(LukeMLM)
library(lattice)
library(lme4)
library(arm)
data(nlsy)
set.seed(1)
nlsysub<-nlsy[nlsy$pubid %in% sample(unique(nlsy$pubid),2000),]
nlsysub$age12
nlsysub$satot

# Without predictor
smod0 <- lmer (satot ~ 1 + (1|pubid), data=nlsysub)
# Unconditional Linear Growth Model: i.e. level-1 predictor is tied to time 
smod1 <- lmer (satot ~ ageint + (1|pubid), data=nlsysub)
# Determine random effect
smod1a <- lmer (satot ~ ageint + (ageint|pubid), data=nlsysub)
smod1b<-lmer(satot ~ age12 + (age12|pubid), data=nlsysub)
anova(smod0,smod1,smod1a,smod1b) # smod1a has smaller AIC, wins

# Center time at age 12 since over age 12 is defined as youth. The advantage is for easier interpretation for the parameter


# Modeling non-linear change patterns
# If you have more than a few timepoints in your longitudinal data, you can move beyond a linear form of
# change and explore curvilinear change. The easiest way is to build models with polynomial time predictors.

# quadratic and cubic models
smod1c <- lmer (satot ~  age12 +  I (age12^2)+(1|pubid), data=nlsysub) 
smod1_2 <- lmer (satot  ~ age12 + I (age12^2) + (age12|pubid), data=nlsysub)
smod1_3 <- lmer (satot  ~ age12 + I (age12^2) + I (age12^3) + (1|pubid), data=nlsysub)
smod1_4 <- lmer (satot  ~ age12 + I (age12^2) + I (age12^3) + (age12|pubid), data=nlsysub)

anova(smod1c,smod1_2,smod1_3,smod1_4) # smod1_4 wins, i.e. cubic model with age as rndom effect wins

#https://cran.r-project.org/web/packages/gamm4/gamm4.pdf
# GAMM models
#summary(nlsysub$satot)
#library(gamm4)
#satot.gm4 <- gamm4 (satot ~ s (ageint),random= ~ (ageint|pubid),data=nlsysub)

#xdum <- seq (0,90,by=1)
#satot.pred <- predict (satot.gm4$gam, data.frame (ageint=xdum),se=TRUE)
#plot (xdum,satot.pred$fit,type='l',col=2,
#      xlab="Age at Interview",ylab="Predicted number of subtance use days ")
# Overdispersion


overdisp_fun <- function(model) {
  ## number of variance parameters in
  ## an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow (m)*( nrow (m)+1)/2
  }
  model.df <- sum ( sapply ( VarCorr (model),vpars))+ length ( fixef (model))
  rdf <- nrow ( model.frame (model))-model.df
  rp <- residuals (model,type="pearson")
  Pearson.chisq <- sum (rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq (Pearson.chisq, df=rdf, lower.tail=FALSE)
  c (chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun (smod1_4) # ration=80.08275 >>1

# GAMM models
summary(nlsysub$satot)
library(gamm4)
smod1_4.gm4 <- gamm4 (satot ~ s (ageint + I (ageint^2) + I (ageint^3)),random= ~ (ageint|pubid),data=nlsysub)
xdum <- seq (0,90,by=1)
satot.pred <- predict (smod1_4.gm4$gam, data.frame (ageint=xdum),se=TRUE)
plot(xdum,satot.pred$fit,type='l',col=2,
      xlab="Age at Interview",ylab="Predicted number of subtance use days",main="Predicted number of substance use days in cubic change model without level-2 predictors")


# Solve overdispersion by zero inflation and negative bionomial model

install.packages ("R2admb")
install.package ("reshape")
install.packages ("glmmADMB",
                  repos= c ("http://glmmadmb.r-forge.r-project.org/repos",getOption ("repos")),
                  type="source")
install.packages ("coefplot2",repos="http://www.math.mcmaster.ca/bolker/R", type="source")
library(glmmADMB)
nlsysub$pubid<-as.factor(nlsysub$pubid)

smod1_4b<-glmmadmb(satot  ~ age12 + I (age12^2) + I (age12^3) + (age12|pubid),
                   data=nlsysub,
                   zeroInflation=TRUE,
                   family="nbinom")
anova(smod1_4,smod1_4b)


#############################
##     Question2     ########
#############################
smod2 <- lmer (satot ~ age12 + sex97 + nonwhite + (age12|pubid), data=nlsysub)
smod3 <- lmer (satot ~ age12*sex97 + age12*nonwhite + sex97*nonwhite+(age12|pubid), data=nlsysub)
smod3a <- lmer (satot ~ age12*sex97 + age12*nonwhite+(age12|pubid), data=nlsysub)
anova(smod2,smod3,smod3a) # ignore gender*race since smod3 and somd3a share same AIC


smod3a <- lmer (satot ~ age12*sex97 + age12*nonwhite+(age12|pubid), data=nlsysub)
smod4 <- lmer (satot ~  age12*nonwhite+(age12|pubid), data=nlsysub)
anova(smod3a,smod4) # smod3a wins, it means sex has influence

smod3a <- lmer (satot ~ age12*sex97 + age12*nonwhite+(age12|pubid), data=nlsysub)
# the quesiton says No include the interaction quadratic or cubic time variables
smod5 <- lmer (satot ~  age12*sex97+age12*nonwhite+ I(age12^2)+I(age12^3)+ (age12|pubid), data=nlsysub)
anova(smod3a,smod5) #smod5 wins
#The prediction graph should include four growth curves: white male,
#white female, nonwhite male, and nonwhite female. Interpret your results.
smod5a <- lmer (satot ~  ageint*sex97+ageint*nonwhite+ I(ageint^2)+I(ageint^3)+ (ageint|pubid), data=nlsysub)


satot_wf <- data.frame(satot = seq (0,30,30/25),ageint= seq(0,25,1),nonwhite = "White",sex97="Female")
satot_wm <- data.frame(satot = seq (0,30,30/25),ageint= seq(0,25,1),nonwhite = "White",sex97="Male")
satot_nwf <- data.frame (satot = seq(0,30,30/25),ageint= seq(0,25,1),nonwhite = "Non-white",sex97="Female")
satot_nwm <- data.frame (satot = seq (0,30,30/25),ageint= seq(0,25,1), nonwhite = "Non-white",sex97="Male")

y_wf <-predict(smod5a,satot_wf,re.form=NA)
y_wm <-predict(smod5a,satot_wm,re.form=NA)
y_nwf<-predict(smod5a,satot_nwf,re.form=NA)
y_nwm<-predict(smod5a,satot_nwm,re.form=NA)

plot (satot_wf$satot,y_wf,type="o",col="blue",lty=1,
      xlim= c(12,25),ylim= c(0,30),
      main="Change in total number of substance use days over time",xlab="Noncentering age",
      ylab=" Predicted number of substance")
points (satot_wm$satot,y_wm,type="l",lty=1,col="brown1")
points (satot_nwf$satot,y_nwf,type="o",lty=1,col="olivedrab")
points (satot_nwm$satot,y_nwm,type="l",lty=1,col="cyan")
legend("topleft",c("White Female","White Male","Non-white Female","Non-white Male"),col=c("blue","brown1","olivedrab","cyan"),lty=c(1,1,1,1))

# "bottomright", "bottom", "bottomleft", 

# "left", "topleft", "top", "topright", 

# "right" and "center"

library(stargazer)

stargazer(smod1b,smod1_4,smod5, title="Three growth models for Substance Use Days",align=TRUE)


#####################
####  Question3  #####
#####################
# Divide the scholl level

nlsysub$mdtx <- c ( lapply ( split (nlsysub, nlsysub$pubid),
                             function(x) cummax (x$schlmdl)),recursive=T)
nlsysub$hstx <- c ( lapply ( split (nlsysub, nlsysub$pubid),
                             function(x) cummax (x$schlhigh)),recursive=T)
nlsysub$cltx <- c ( lapply ( split (nlsysub, nlsysub$pubid),
                             function(x) cummax (x$schlclg)),recursive=T)
nlsysub2<-data.frame(nlsysub)
nlsysub2$mdtx


smod8a <- lmer (satot ~ age12+hstx + (age12|pubid), data=nlsysub2)
smod8b <- lmer (satot ~ age12+cltx + (age12|pubid), data=nlsysub2)
smod8c <- lmer (satot ~ age12+hstx+cltx + (age12|pubid), data=nlsysub2)

anova(smod8a,smod8b,smod8c) #it turns out the transition to college accounts for a larger shift in expected substance use
# smod8b wins

smod9a <- lmer (satot ~ age12*hstx + (age12|pubid), data=nlsysub2)
smod9b <- lmer (satot ~ age12*cltx + (age12|pubid), data=nlsysub2)
smod9c <- lmer (satot ~ age12*hstx+age12*cltx + (age12|pubid), data=nlsysub2)

anova(smod9a,smod9b,smod9c) 
stargazer(smod9a,smod9b,smod9c, title="Three growth models for Substance Use Days with education levels",align=TRUE)

## Use smod9a to predict satot for a person who enters high school at age 15
## Use smod9b to predict satot for a person who enters high school at age 18

summary(smod9b)

xbase <- seq (12,18,.5)
xclg <- seq (18,25,.5)
ybase <-  -0.02673  + 1.65035 *(xbase-12)
yclg <- ( -0.02673 + 1.65035 *3 +5.55040) + (1.65035 -0.59114)*(xclg-18)
yclgb <- ( -0.02673 + 1.65035 *3) + 1.65035*(xclg-18)
plot (xbase,ybase,type="n",xlab="Age",ylab="Predict Substance",
      main="Shift of Substance Use Days at college Transition",xlim= c (12,25),ylim= c (0,25))
points (xbase,ybase,type="o",col="red",lty=1,lwd=3)
points (xclg,yclg,type="l",col="black",lty=1,lwd=3)
points (xclg,yclgb,type="l",col="blue",lty=2,lwd=3)
segments (18,(-0.02673 + 1.65035 *3 ),18,(-0.02673 + 1.65035 *3 +5.55040),lty=1,lwd=3)
abline(v=c(15,18),lty=2)
#legend("topleft",c("in High school transition","W/ college transition","w/o college transition & w/ high school"),col=c("red","black","blue"))



## ## Use smod9a to predict satot for a person who high school


summary(smod9a)

xbase <- seq (12,15,.5)
xclg <- seq (15,25,.5)
ybase <- 0.17282  + 1.59726 *(xbase-12)
yclg <- (0.17282 + 1.59726*3-0.11756) + (1.59726+ 0.12559)*(xclg-15)
yclgb <- (0.17282 + 1.59726*3) + 1.59726*(xclg-15)
plot (xbase,ybase,type="n",xlab="Age",ylab="Predict Substance",
      main="Shift of Substance Use Days at High School Transition",xlim= c (11,25),ylim= c (0,25))
points (xbase,ybase,type="o",col="red",lty=1,lwd=3)
points (xclg,yclg,type="l",col="black",lty=1,lwd=3)
points (xclg,yclgb,type="l",col="blue",lty=2,lwd=3)
segments (15,(0.17282 + 1.59726*3),15,(0.17282 + 1.59726*3-0.11756),lty=1,lwd=3)
abline (v=15,lty=2)
#legend("topleft",c("w/o High school transition","W/ high school transition","w/o college transition & w/ high school"),col=c("red","black","blue"))


print ( xyplot (satot ~ ageint | pubid, type= c ("p","r"),
                data=nlsysub[nlsysub$pubid %in% sample ( unique (nlsysub$pubid),20),]))
