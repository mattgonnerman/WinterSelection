require(ggplot2)
require(dplyr)
require(tidyverse)
require(hexbin)

######################
## Prepare the data ##
######################
# Raw data
data <- turkeyData.zm
step <- data$step
angle <- data$angle

# MomentuHMM model
m <- hmm.top.model

# Animal Index
nbAnimals <- length(unique(m$data$ID))
animalsInd <- 1:nbAnimals
ID <- unique(m$data$ID)[animalsInd]

# States decoding with Viterbi 
states <- viterbi(m)

#######################################
## Plot the Step Length Density/Hist ##
#######################################
turkey_states <- turkey_states %>%
  filter(!is.na(location_lat)) %>%
  mutate(State = ifelse(State == 1, "Roosting", 
                 ifelse(State == 2, "Stationary",
                 ifelse(State == 3, "Mobile", NA)))) %>%
  mutate(State = factor(State, levels = c("Roosting", "Stationary", "Mobile"))) 


ggplot(turkey_states, aes(x = step, group = State)) +
  geom_density(aes(color = State, fill = State), alpha = .7, size = 1.2) +
  theme_classic(base_size = 25) +
  scale_y_continuous(limits=c(0,.04), oob = scales::squish) + 
  scale_x_continuous(limits=c(0,500)) +
  scale_color_manual(name = "Behavioral\nState",
                      # labels = c("Roost", "Stationary", "Mobile"),
                      values = c("#1cade4", "#f1a806", "#46e300")) +
  scale_fill_manual(name = "Behavioral\nState",
                      # labels = c("Roost", "Stationary", "Mobile"),
                      values = c("#1cade4", "#f1a806", "#46e300")) + 
  theme(legend.position = c(0.77, 0.77)) +
  xlab("Step Length") +
  ylab("Density")
ggsave("Results/StepLengthDensity.jpeg", width = 7, height = 7, units = "in")

ggplot(turkey_states, aes(x = angle, group = State)) +
  geom_density(aes(color = State), alpha = .7, size = 1.2) +
  theme_classic(base_size = 25) +
  # scale_y_continuous(limits=c(0,.04), oob = scales::squish) + 
  scale_x_continuous(limits=c(-pi,pi)) +
  scale_color_manual(name = "Behavioral\nState",
                     # labels = c("Roost", "Stationary", "Mobile"),
                     values = c("#f1a806", "#46e300")) +
  scale_fill_manual(name = "Behavioral\nState",
                    # labels = c("Roost", "Stationary", "Mobile"),
                    values = c("#f1a806", "#46e300")) + 
  theme(legend.position = c(0.5, 0.3)) +
  xlab("Turning Angle") +
  ylab("Density")
ggsave("Results/TurningAngleDensity.jpeg", width = 7, height = 7, units = "in")


##################################################
## Plot the t.p. as functions of the covariates ##
##################################################
#Load MLE of Transition Porbability betas
TPbetas <- read.csv('Results/HMM - MLE of betas.csv')
row.names(TPbetas) <- TPbetas$X
TPbetas <- TPbetas %>% select(-X)
TPbetas <- as.data.frame(t(TPbetas)) %>%
  dplyr::select('(Intercept)', WC.Z, SD.Z, 'cosinorCos(hour, period = 24)', 'cosinorSin(hour, period = 24)')

X.wc <- seq(-2,.5,.1)
X.sd <- seq(-.5,6,.1)
X.coshour <- cos(2*pi*seq(0,24, .25)/24)
X.sinhour <- sin(2*pi*seq(0,24, .25)/24)

x <-1

eq <- TPbetas[x,1] + TPbetas[x,2]*-1 + TPbetas[x,3]*X.sd + TPbetas[x,4]*cos(2*pi*11.54/24) + TPbetas[x,5]*sin(2*pi*11.54/24)
plot((exp(eq)/(1+exp(eq)))~X.sd)


require(VGAM)
pneumo <- transform(pneumo, let = log(exposure.time))
fit <- vglm(cbind(normal, mild, severe) ~ let,
            multinomial, trace = TRUE, data = pneumo)  # For illustration only!
fitted(fit)
predict(fit)

mlogit(fitted(fit))
mlogit(fitted(fit)) - predict(fit)  # Should be all 0s

multilogit(predict(fit), inverse = TRUE)  # rowSums() add to unity
multilogit(predict(fit), inverse = TRUE, refLevel = 1)  # For illustration only
multilogit(predict(fit), inverse = TRUE) - fitted(fit)  # Should be all 0s

mlogit(fitted(fit), deriv = 1)
mlogit(fitted(fit), deriv = 2)


##################################################
## Plot the t.p. as functions of the covariates ##
##################################################

covs=NULL
meansList<-c("matrix","numeric","integer","logical","Date","POSIXlt","POSIXct","difftime")
if(is.null(covs)){
    covs <- m$data[which(m$data$ID %in% ID),][1,]
    for(j in names(m$data)[which(unlist(lapply(m$data,function(x) any(class(x) %in% meansList))))]){
      if(inherits(m$data[[j]],"angle")) covs[[j]] <- CircStats::circ.mean(m$data[[j]][which(m$data$ID %in% ID)][!is.na(m$data[[j]][which(m$data$ID %in% ID)])])
      else covs[[j]]<-mean(m$data[[j]][which(m$data$ID %in% ID)],na.rm=TRUE)
    }
}

gamInd<-(length(m$mod$estimate)-(nbCovs+1)*nbStates*(nbStates-1)*mixtures+1):(length(m$mod$estimate))-(ncol(m$covsPi)*(mixtures-1))-ifelse(nbRecovs,nbRecovs+1+nbG0covs+1,0)-ncol(m$covsDelta)*(nbStates-1)*(!m$conditions$stationary)*mixtures
quantSup<-qnorm(1-(1-alpha)/2)
  
# values of each covariate
rawCovs <- m$rawCovs
if(is.null(covs)) {
  rawCovs <- m$rawCovs
  meanCovs <- colSums(rawCovs)/nrow(rawCovs)
  } else {
    rawCovs <- m$data[,names(covs),drop=FALSE]
    meanCovs <- as.numeric(covs)
  }
covIndex <- 1:ncol(rawCovs)

covNames <- getCovNames(m,p,i)
DMterms<-covNames$DMterms
 
for(cov in covIndex) {
  if(!is.factor(rawCovs[,cov])){
     gridLength <- 101
     hGridLength <- gridLength*ifelse(inherits(m,"hierarchical"),nlevels(m$data$level),1)
        
        inf <- min(rawCovs[,cov],na.rm=TRUE)
        sup <- max(rawCovs[,cov],na.rm=TRUE)
        
        # set all covariates to their mean, except for "cov"
        # (which takes a grid of values from inf to sup)
        tempCovs <- data.frame(matrix(covs[names(rawCovs)][[1]],nrow=hGridLength,ncol=1))
        if(ncol(rawCovs)>1)
          for(i in 2:ncol(rawCovs))
            tempCovs <- cbind(tempCovs,rep(covs[names(rawCovs)][[i]],gridLength))
        
        tempCovs[,cov] <- rep(seq(inf,sup,length=gridLength),each=hGridLength/gridLength)
      } else {
        gridLength<- nlevels(rawCovs[,cov])
        hGridLength <- gridLength*ifelse(inherits(m,"hierarchical"),nlevels(m$data$level),1)
        # set all covariates to their mean, except for "cov"
        tempCovs <- data.frame(matrix(covs[names(rawCovs)][[1]],nrow=hGridLength,ncol=1))
        if(ncol(rawCovs)>1)
          for(i in 2:ncol(rawCovs))
            tempCovs <- cbind(tempCovs,rep(covs[names(rawCovs)][[i]],gridLength))
        
        tempCovs[,cov] <- as.factor(rep(levels(rawCovs[,cov]),each=hGridLength/gridLength))
      }
      
      names(tempCovs) <- names(rawCovs)
      tmpcovs<-covs[names(rawCovs)]
      for(i in which(unlist(lapply(rawCovs,is.factor)))){
        tempCovs[[i]] <- factor(tempCovs[[i]],levels=levels(rawCovs[,i]))
        tmpcovs[i] <- as.character(tmpcovs[[i]])
      }
      for(i in which(!unlist(lapply(rawCovs,is.factor)))){
        tmpcovs[i]<-round(covs[names(rawCovs)][i],2)
      }
      if(!is.null(recharge)){
        tmprecovs<-covs[names(m$reCovs)]
        for(i in which(unlist(lapply(m$reCovs,is.factor)))){
          tmprecovs[i] <- as.character(tmprecovs[[i]])
        }
        for(i in which(!unlist(lapply(m$reCovs,is.factor)))){
          tmprecovs[i]<-round(recovs[names(m$reCovs)][i],2)
        }
      }
      
      if(inherits(m$data,"hierarchical")) class(tempCovs) <- append("hierarchical",class(tempCovs))
      
      tmpSplineInputs<-getSplineFormula(newformula,m$data,tempCovs)
      desMat <- stats::model.matrix(tmpSplineInputs$formula,data=tmpSplineInputs$covs)
      
      for(mix in 1:mixtures){
        
        if(is.null(recharge)){
          trMat <- trMatrix_rcpp(nbStates,beta$beta[(mix-1)*(nbCovs+1)+1:(nbCovs+1),,drop=FALSE],desMat,m$conditions$betaRef)
        } else {
          trMat <- array(unlist(lapply(split(tmpSplineInputs$covs,1:nrow(desMat)),function(x) tryCatch(get_gamma_recharge(m$mod$estimate[c(gamInd[unique(c(m$conditions$betaCons))],length(m$mod$estimate)-nbRecovs:0)],covs=x,formula=tmpSplineInputs$formula,hierRecharge=hierRecharge,nbStates=nbStates,betaRef=m$conditions$betaRef,betaCons=m$conditions$betaCons,workBounds=rbind(m$conditions$workBounds$beta,m$conditions$workBounds$theta),mixture = mix),error=function(e) NA))),dim=c(nbStates,nbStates,nrow(desMat)))
        }
        
        if(!inherits(m,"hierarchical")){
          
          plotTPM(nbStates,cov,ref=1:nbStates,tempCovs,trMat,rawCovs,lwd,arg,plotCI,Sigma,gamInd,m,desMat,nbRecovs,tmpSplineInputs$formula,hierRecharge,mix,muffWarn,quantSup,tmpSplineInputs$covs,stateNames=1:nbStates)
          
          txt <- paste(names(rawCovs)[-cov],"=",tmpcovs[-cov],collapse=", ")
          if(nbRecovs & names(rawCovs)[cov]=="recharge"){
            tmpNames <- c(names(rawCovs)[-cov],colnames(m$reCovs))
            txt <- paste(tmpNames[!duplicated(tmpNames)],"=",c(tmpcovs[-cov],tmprecovs)[!duplicated(tmpNames)],collapse=", ")
          }
          if(ncol(rawCovs)>1 | nbRecovs) do.call(mtext,c(list(paste0(ifelse(mixtures>1,paste0("Mixture ",mix," t"),"T"),"ransition probabilities",ifelse(nbRecovs," at next time step: ",": "),txt),side=3,outer=TRUE,padj=2,cex=cex.main),marg))
          else do.call(mtext,c(list(paste0(ifelse(mixtures>1,paste0("Mixture ",mix," t"),"T"),"ransition probabilities"),side=3,outer=TRUE,padj=2,cex=cex.main),marg))
          
        } else {
          for(j in 1:(m$conditions$hierStates$height-1)){
            
            txt <- paste(names(rawCovs)[-cov][which(names(rawCovs)[-cov]!="level")],"=",tmpcovs[which(tmpcovs$level==j),-cov][which(names(rawCovs)[-cov]!="level")],collapse=", ")
            if(nbRecovs & grepl("recharge",names(rawCovs)[cov])){
              tmpNames <- c(names(rawCovs)[-cov][which(names(rawCovs)[-cov]!="level" & !grepl("recharge",names(rawCovs)[-cov]))],colnames(m$reCovs)[which(colnames(m$reCovs)!="level" & !grepl("recharge",colnames(m$reCovs)))])
              tmprecovs <- tmprecovs[,which(colnames(tmprecovs)!="level" & !grepl("recharge",colnames(tmprecovs))),drop=FALSE]
              txt <- paste(tmpNames[!duplicated(tmpNames)],"=",c(tmpcovs[which(tmpcovs$level==j),-cov][which(names(rawCovs)[-cov]!="level" & !grepl("recharge",names(rawCovs)[-cov]))],tmprecovs)[!duplicated(tmpNames)],collapse=", ")
            }
            
            if(j==1) {
              
              ref <- m$conditions$hierStates$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==j+1)
              
              # only plot if there is variation in stationary state proabilities
              if(!all(apply(trMat[ref,ref,which(tempCovs$level==j)],1:2,function(x) all( abs(x - mean(x)) < 1.e-6 )))){
                
                plotTPM(nbStates,cov,ref,tempCovs[which(tempCovs$level==j),],trMat[,,which(tempCovs$level==j)],rawCovs,lwd,arg,plotCI,Sigma,gamInd,m,desMat[which(tempCovs$level==j),],nbRecovs,tmpSplineInputs$formula,hierRecharge,mix,muffWarn,quantSup,tmpSplineInputs$covs[which(tempCovs$level==j),],stateNames=names(ref))
                
                if(ncol(rawCovs[-cov])>1 | nbRecovs) do.call(mtext,c(list(paste0(ifelse(mixtures>1,paste0("Mixture ",mix," t"),"T"),"ransition probabilities for level",j,ifelse(nbRecovs," at next time step: ",": "),txt),side=3,outer=TRUE,padj=2,cex=cex.main),marg))
                else do.call(mtext,c(list(paste0(ifelse(mixtures>1,paste0("Mixture ",mix," t"),"T"),"ransition probabilities for level",j),side=3,outer=TRUE,padj=2,cex=cex.main),marg))
                
                #if(length(covnames[-cov])>1) do.call(mtext,c(list(paste0(ifelse(mixtures>1,paste0("Mixture ",mix," s"),"S"),"tationary state probabilities for level",j,": ",paste(covnames[-cov][which(covnames[-cov]!="level")]," = ",tmpcovs[which(tmpcovs$level==j),-cov][which(covnames[-cov]!="level")],collapse=", ")),side=3,outer=TRUE,padj=2,cex=cex.main),marg))
                #else do.call(mtext,c(list(paste0(ifelse(mixtures>1,paste0("Mixture ",mix," s"),"S"),"tationary state probabilities for level",j),side=3,outer=TRUE,padj=2,cex=cex.main),marg))
              } 
            } else {
              t <- data.tree::Traverse(m$conditions$hierStates,filterFun=function(x) x$level==j)
              names(t) <- m$conditions$hierStates$Get("name",filterFun=function(x) x$level==j)
              for(k in names(t)){
                ref <- t[[k]]$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==j+1)#t[[k]]$Get("state",filterFun = data.tree::isLeaf)
                # only plot if jth node has children and there is variation in stationary state proabilities
                if(!is.null(ref) && !all(apply(trMat[ref,ref,which(tempCovs$level==j)],1:2,function(x) all( abs(x - mean(x)) < 1.e-6 )))){
                  
                  plotTPM(nbStates,cov,ref,tempCovs[which(tempCovs$level==j),],trMat[,,which(tempCovs$level==j)],rawCovs,lwd,arg,plotCI,Sigma,gamInd,m,desMat[which(tempCovs$level==j),],nbRecovs,tmpSplineInputs$formula,hierRecharge,mix,muffWarn,quantSup,tmpSplineInputs$covs[which(tempCovs$level==j),],stateNames=names(ref))
                  
                  if(ncol(rawCovs[-cov])>1 | nbRecovs) do.call(mtext,c(list(paste0(ifelse(mixtures>1,paste0("Mixture ",mix," t"),"T"),"ransition probabilities for level",j," ",k,ifelse(nbRecovs," at next time step: ",": "),txt),side=3,outer=TRUE,padj=2,cex=cex.main),marg))
                  else do.call(mtext,c(list(paste0(ifelse(mixtures>1,paste0("Mixture ",mix," t"),"T"),"ransition probabilities for level",j," ",k),side=3,outer=TRUE,padj=2,cex=cex.main),marg))
                  
                }
              }
            }
          }
        }
      }
    }



#####################################################
###Trying to think of cool ways to visualize data ###
#####################################################
#Weather Data
TPdata <- read.csv("Results/HMMBehavioralStates_output.csv") %>%
  dplyr::select(ID, Timestamp, step, angle, WC.Z, SD.Z, hour, State) %>%
  mutate(State = as.factor(State))

Sdata <- TPdata %>% filter(State == 2)
Mdata <- TPdata %>% filter(State == 3)

ggplot(Sdata, aes(x = step, y = angle)) +
  geom_hex(bins = 50) +
  scale_fill_viridis_c(option = "magma") +
  theme_classic()

ggplot(Mdata, aes(x = step, y = angle)) +
  geom_hex(bins = 50) +
  scale_fill_viridis_c(option = "magma") +
  theme_classic() + 
  xlim(0, 1000)


ggplot(TPdata, aes(x = log(step), group = State)) +
  geom_density(aes(color = State, fill = State), alpha = .4) +
  theme_classic() 