# To run, use the command:
# sbatch job.sbatch


library(data.table)
library(SuperLearner)
library(cvAUC)
library(doMPI)

set.seed(1)

# load the data
data = as.data.table(read.csv(
	'Zoster_Merck_Finaldata_protocol022_March2012_mock.csv',
	stringsAsFactors=FALSE))

# Load the full data set (not just case-cohort)
data.pp = as.data.table(read.csv(
	'p022_v2_mock.csv',
	stringsAsFactors=FALSE))

data.pp[,immunocoh:=as.numeric(usubjid%in%data$USUBJID[data$SUBCOHFL=='Y']),]
data.pp[,wgt.init:=immunocoh/predict(glm(immunocoh~factor(plan_arm)*sex*age,data=subset(data.pp,hz_confirm==0),family=binomial),newdata=data.pp,type='response'),]

# Inverse weight
data[,wgt:=F_ITT + (1-F_ITT) * data.pp$wgt.init[match(USUBJID,data.pp$usubjid)],]

rm(data.pp)

# Transform pre/post titer using log10
data[,LBPRE.log10:=log10(LBPRE),]
data[,LBPOST.log10:=log10(LBPOST),]

# Add average PRE/POST log10 titer
data[,LBAVG.log10:=(LBPRE.log10+LBPOST.log10)/2,]

# Add fold-rise log10 titer
data[,LB.fr.log10:=LBPOST.log10-LBPRE.log10,]

# Remove unlogged titers
data[,c('LBPRE','LBPOST'):=NULL]

# Rename treatment VACC
data[,VACC:=ACTLARM,]
data[,ACTLARM:=NULL,]

# Make sex numeric (1 => Female)
data[,SEXbin:=as.numeric(SEX=='F'),]
data[,SEX:=NULL,]

# Named vector of candidate predictor columns in data
predictors = c('AGE','LBPRE.log10','LBPOST.log10','LB.fr.log10','LBAVG.log10','SEXbin')

# super-learner library (prediction algorithms only -- screening algorithms to be included separately)
# NOTE: did not include SL.xgboost in this mock analysis because fails in this fake data set when titer information not included (worked on the real data set)
SL.pred.algs = c('SL.bayesglm','SL.glm','SL.glm.interaction','SL.mean','SL.gam','SL.cforest')#,'SL.xgboost')

# Remove average titers
screen.normal <- function(Y, X, family, obsWeights, id, ...) {
  vars <- rep(TRUE, ncol(X))
  vars[names(X) %in% c('LBAVG.log10','LB.fr.log10')] <- FALSE
  vars
}

# Screening algorithm that removes PRE/AVG titers
screen.noPRE <- function(Y, X, family, obsWeights, id, ...) {
  vars <- rep(TRUE, ncol(X))
  vars[names(X) %in% c('LBPRE.log10','LBAVG.log10','LB.fr.log10')] <- FALSE
  vars
}

# Screening algorithm that removes POST titer
screen.noPOST <- function(Y, X, family, obsWeights, id, ...) {
  vars <- rep(TRUE, ncol(X))
  vars[names(X) %in% c('LBPOST.log10','LBAVG.log10','LB.fr.log10')] <- FALSE
  vars
}

# Screening algorithm that removes AVG, POST titers
screen.PREfr <- function(Y, X, family, obsWeights, id, ...) {
  vars <- rep(TRUE, ncol(X))
  vars[names(X) %in% c('LBPOST.log10','LBAVG.log10')] <- FALSE
  vars
}

# Screening algorithm that removes PRE/POST titers (only uses average)
screen.AVG <- function(Y, X, family, obsWeights, id, ...) {
  vars <- rep(TRUE, ncol(X))
  vars[names(X) %in% c('LBPRE.log10','LBPOST.log10','LB.fr.log10')] <- FALSE
  vars
}

# Screening algorithm that removes all titers (reference)
screen.noTiter <- function(Y, X, family, obsWeights, id, ...) {
  vars <- rep(TRUE, ncol(X))
  vars[names(X) %in% c('LBPRE.log10','LBPOST.log10','LBAVG.log10','LB.fr.log10')] <- FALSE
  vars
}

# Function that takes as input SL.pred.algs and a vector of screening algorithms and outputs
# a list of vectors that can be input into SuperLearner as SL.library

make.SL.library = function(pred.algs,screens){lapply(pred.algs,function(x){if(x=='SL.mean'){ c('SL.mean','screen.normal') } else {c(x,screens)}})}

# Screening options to run
screen.options = list(
	noTiter=c('screen.noTiter'), # No titers
	noPOST=c('screen.noPOST'), # No POST titers
	noPRE=c('screen.noPRE'), # No PRE titers
	All=c('screen.normal','screen.noPRE','screen.noPOST','screen.PREfr','screen.AVG') # All titers
	)

# Uncomment below if want to save SuperLearner fit on full data set
# SL.out = list()
cvAUC.out = list()
AUCplot = list()
predByCase = list()

for(i in seq(screen.options)){
	print(i)
	# Current screening library
	screen.lib = screen.options[[i]]

	cvAUC.out[[i]] <- AUCplot[[i]] <- predByCase[[i]] <- list()
	# Uncomment below if want to save SuperLearner fit on full data set
	# SL.out[[i]] <- list()

	for(j in seq(unique(data$VACC))){
		A = unique(data$VACC)[j]
		data.curr = subset(data,VACC==A)
		print(A)
		####
		# Uncomment below if want to save SuperLearner fit on full data set
		# SL fits on the full data set (run ten times to minimize fold dependence)
		# tmp = lapply(1:10,function(k){
		# 	SL = with(data.curr,
		# 		SuperLearner(
		# 			Y = F_ITT,
		# 			X = as.data.frame(subset(data.curr,select=predictors)),
		# 			obsWeights=wgt,
		# 			SL.library = make.SL.library(SL.pred.algs,screen.lib),
		# 			family = binomial(),
		# 			verbose=FALSE,
		# 			control=list(saveFitLibrary = FALSE)))
		# 	list(preds=SL$SL.predict[,1],coef=SL$coef,risk=SL$cvRisk)
		# })
		# save(tmp,file='results/tmp.Rdata')
		# SL.out[[i]][[j]] = list(preds=rowMeans(sapply(tmp,function(x){x$preds})),coef=rowMeans(sapply(tmp,function(x){x$coef})),risk=rowMeans(sapply(tmp,function(x){x$risk})))
		# rm(tmp)

		####
		# Stratify V-fold cross-validation so that all validation samples have roughly the same number of events
		.cvFolds <- function(Y, V){ # Create CV folds (stratify by outcome) -- code from cvAUC documentation
			Y0 <- split(sample(which(Y==0)), rep(1:V, length=length(which(Y==0))))
			Y1 <- split(sample(which(Y==1)), rep(1:V, length=length(which(Y==1))))
			folds <- vector("list", length=V)
			for (v in seq(V)) {folds[[v]] <- c(Y0[[v]], Y1[[v]])}
			return(folds)
		}
		V = 15
		V.folds = .cvFolds(data.curr$F_ITT,V) # What we'd use if we were doing V-fold CV
		V2.folds = apply(combn(seq(V),2),2,function(x){c(V.folds[[x[1]]],V.folds[[x[2]]])}) # V choose 2 fold
		
		sapply(V2.folds,function(x){sum(data.curr$F_ITT[x])})

		# Start a cluster with 20 cores
		cl <- startMPIcluster(20)
		registerDoMPI(cl)
		exportDoMPI(cl,varlist=c('data.curr','V2.folds','predictors','make.SL.library','SL.pred.algs','screen.lib','screen.normal','screen.noPRE','screen.noPOST','screen.PREfr','screen.AVG','screen.noTiter'))
		# CV.SL fits
		CV.SL = foreach(i=1:length(V2.folds),.packages=c('SuperLearner','data.table'), .verbose=TRUE, .errorhandling='pass') %dopar% {
			fold = V2.folds[[i]]
			SL.tmp = with(data.curr[-fold,],
				SuperLearner(
					Y = F_ITT,
					X = as.data.frame(subset(data.curr[-fold,],select=predictors)),
					obsWeights=wgt,
					SL.library = make.SL.library(SL.pred.algs,screen.lib),
					family = binomial(),
					verbose=FALSE,
					newX = as.data.frame(subset(data.curr[fold,],select=predictors)),
					control=list(saveFitLibrary = FALSE)))
			out = cbind(SL.tmp$SL.predict,SL.tmp$library.predict)
			colnames(out)[1] = 'Super Learner'
			return(out)
		}
		closeCluster(cl)
		
		# data frame of cross-validated predictions and case status (j stratifies this data frame by vaccination status)
		predByCase[[i]][[j]] = data.frame(
			unlist(lapply(CV.SL,function(x){x[,1]})),
			unlist(lapply(V2.folds,function(fold){data.curr$F_ITT[fold]})))

		AUCplot[[i]][[j]] = cvAUC(lapply(CV.SL,function(x){x[,1]}),lapply(V2.folds,function(fold){data.curr$F_ITT[fold]}))$perf
		
		cvAUC.out[[i]][[j]] = t(sapply(1:ncol(CV.SL[[1]]),function(ii){
			tmp = ci.cvAUC(lapply(CV.SL,function(x){x[,ii]}),lapply(V2.folds,function(fold){data.curr$F_ITT[fold]}))
			truese = sqrt(length(unlist(V2.folds))/length(data.curr$F_ITT)) * tmp$se # Adjust se to account for V2.fold cross-validation
			c(tmp$cvAUC,tmp$cvAUC+c(-1,1)*qnorm(0.975)*truese)
		}))
		rownames(cvAUC.out[[i]][[j]]) = colnames(CV.SL[[1]])
		rm(CV.SL)
	}
	names(cvAUC.out[[i]]) <- names(AUCplot[[i]]) <- names(predByCase[[i]]) <- unique(data$VACC)
	# Uncomment below if want to save SuperLearner fit on full data set
	# names(SL.out[[i]]) <- unique(data$VACC)
}

names(cvAUC.out) <- names(AUCplot) <- names(predByCase) <- names(screen.options)
# Uncomment below if want to save SuperLearner fit on full data set
# names(SL.out) <- names(screen.options)

save.image(file='results/outs.Rdata')

mpi.finalize()
