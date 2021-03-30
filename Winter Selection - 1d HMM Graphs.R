require(ggplot2)
require(dplyr)
require(tidyverse)
require(hexbin)

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

x <-6

eq <- TPbetas[x,1] + TPbetas[x,2]*-1 + TPbetas[x,3]*X.sd + TPbetas[x,4]*cos(2*pi*11.54/24) + TPbetas[x,5]*sin(2*pi*11.54/24)
plot((exp(eq)/(1+exp(eq)))~X.sd)

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




##################################################
## Plot the t.p. as functions of the covariates ##
##################################################
m <- turk_m1.zm
reForm <- formatRecharge(nbStates,m$conditions$formula,m$conditions$betaRef,m$data,par=m$mle)
recharge <- reForm$recharge
hierRecharge <- reForm$hierRecharge
newformula <- reForm$newformula



if(nbStates>1) {
  par(mar=c(5,4,4,2)-c(0,0,1.5,1)) # bottom, left, top, right
  
  gamInd<-(length(m$mod$estimate)-(nbCovs+1)*nbStates*(nbStates-1)*mixtures+1):(length(m$mod$estimate))-(ncol(m$covsPi)*(mixtures-1))-ifelse(nbRecovs,nbRecovs+1+nbG0covs+1,0)-ncol(m$covsDelta)*(nbStates-1)*(!m$conditions$stationary)*mixtures
  quantSup<-qnorm(1-(1-alpha)/2)
  
  if(nbCovs>0) {
    
    # values of each covariate
    #rawCovs <- m$rawCovs[which(m$data$ID %in% ID),,drop=FALSE]
    #if(is.null(covs)) {
    #  rawCovs <- m$rawCovs
    #  meanCovs <- colSums(rawCovs)/nrow(rawCovs)
    #} else {
    #  rawCovs <- m$data[,names(covs),drop=FALSE]
    #  meanCovs <- as.numeric(covs)
    #}
    
    if(inherits(m,"hierarchical")) {
      covIndex <- which(!(names(rawCovs)=="level"))
      covs$level <- NULL
      covs <- data.frame(covs[rep(1:nrow(covs),nlevels(m$data$level)),,drop=FALSE],level=rep(levels(m$data$level),each=nrow(covs)))
    } else covIndex <- 1:ncol(rawCovs)
    
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
      
      if(plotStationary) {
        par(mfrow=c(1,1))
        if(inherits(m,"hierarchical")){
          if(is.null(recharge)){
            tmpSplineInputs$covs <- tempCovs
          } else {
            tmpSplineInputs$covs <- tmpSplineInputs$covs[which(tmpSplineInputs$covs$level==levels(m$data$level)[1]),]
          }
        }
        statPlot(m,Sigma,nbStates,tmpSplineInputs$formula,tmpSplineInputs$covs,tempCovs,tmpcovs,cov,hierRecharge,alpha,gridLength,gamInd,names(rawCovs),col,plotCI,...)
      }
    }
  }
}