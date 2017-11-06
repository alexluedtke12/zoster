
load('results/outs.Rdata')
rm(cl)

library(ggplot2)
library(plyr)

my.colors=c('red','blue','green','black')

num.folds = length(AUCplot[[1]][['Placebo']]@alpha.values)
xseq = seq(0,1,length=1e4)

############################################################
# ROC Plot
dat = do.call(rbind,
	lapply(c('Placebo','ZOSTAVAX'),function(VACC){
		Vacc = 	paste0(toupper(substr(VACC, 1, 1)), tolower(substr(VACC, 2, nchar(VACC))))
		do.call(rbind,lapply(1:length(AUCplot),function(j){
			sfs = lapply(1:num.folds,function(i){
				approxfun(AUCplot[[j]][[VACC]]@x.values[[i]],c(AUCplot[[j]][[VACC]]@y.values[[i]]))
			})
			sf = function(xx){rowMeans(sapply(sfs,function(ff){ff(xx)}))}
			data.frame(x=xseq,y=sf(xseq),scrn=names(AUCplot)[j],arm=Vacc)
		}))
	})
)

dat$scrn = revalue(dat$scrn,c('noTiter'='None','noPOST'='Day 1 Only','noPRE'='Week 6 Only','All'='Day 1, Week 6'))

dat$arm = factor(dat$arm)

p = ggplot(dat,aes(x=x,y=y,colour=scrn)) + theme_bw() + geom_line() + facet_grid(arm~.) + geom_abline(colour='grey') + xlab('False Positive Rate') + ylab('Cross-Validated True Positive Rate') + labs(colour='Titers Included') + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme(panel.spacing = unit(1, "lines"))
ggsave('figures/ROC.pdf',plot=p)

############################################################
# CV-AUC

dat2 = do.call(rbind,
	lapply(c('Placebo','ZOSTAVAX'),function(VACC){
		Vacc = 	paste0(toupper(substr(VACC, 1, 1)), tolower(substr(VACC, 2, nchar(VACC))))
		do.call(rbind,lapply(1:length(cvAUC.out),function(j){
			out = data.frame(cvAUC.out[[j]][[VACC]])
			colnames(out) = c('CVAUC','xmin','xmax')
			alg = c(rownames(out)[1],
				sapply(strsplit(rownames(out)[-1],'SL.'),function(x){strsplit(x,'_')[[2]]})[1,])
			VarRed = c('None',
				sapply(strsplit(rownames(out)[-1],'SL.'),function(x){strsplit(x,'screen.')[[2]]})[2,])
			data.frame(out,Algorithm=alg,VarReduction=VarRed,scrn=names(AUCplot)[j],arm=Vacc)
		}))
	})
)

dat2$Algorithm = revalue(dat2$Algorithm,c('Super Learner'='Super Learner','bayesglm'='Bayes GLM','gam'='GAM','glm'='GLM','glm.interaction'='GLM, Interactions','knn'='k-Nearest Neighbors','mean'='Sample Mean','polymars'='Polynomial Splines','step'='Stepwise Regression','xgboost'='Gradient Boosting','cforest'='Random Forest','nnet'='Neural Network'))
dat2$Algorithm = factor(dat2$Algorithm,levels=c(levels(dat2$Algorithm)[levels(dat2$Algorithm)!='Super Learner'],'Super Learner'))
dat2$VarReduction = revalue(dat2$VarReduction,c('None'='N/A (Super Learner)','noTiter'='No Titers','noPOST'='Day 1 Only','noPRE'='Week 6 Only','PREfr'='Fold-Rise Only','normal'='Both Titers','AVG'='Average Only'))
dat2$scrn = revalue(dat2$scrn,c('noTiter'='None','noPOST'='Day 1 Only','noPRE'='Week 6 Only','All'='Day 1, Week 6','AVG'='Average Only'))

p2 = ggplot(dat2,aes(x=CVAUC,y=Algorithm,colour=VarReduction)) + theme_bw() + facet_grid(scrn~arm) + geom_point(size=0.7,alpha=0.85,position=position_jitter(w = 0, h = 0.05)) + geom_errorbarh(aes(xmax=xmax,xmin=xmin),alpha=0.85) + geom_vline(xintercept=0.5,colour='grey') + xlab('Cross-Validated Area Under the Curve') + labs(colour='Variable Screening')
ggsave('figures/AUC.pdf',plot=p2)

############################################################
# Boxplot

dat3 = do.call(rbind,
	lapply(c('Placebo','ZOSTAVAX'),function(VACC){
		Vacc = 	paste0(toupper(substr(VACC, 1, 1)), tolower(substr(VACC, 2, nchar(VACC))))
		do.call(rbind,lapply(1:length(predByCase),function(j){
			out = data.frame(predByCase[[j]][[VACC]])
			colnames(out) = c('PredProb','CaseStatBin')
			tmp = rep('Zoster\nControl',nrow(out))
			tmp[out$CaseStatBin==1] = 'Zoster\nCase'
			data.frame(PredProb=out$PredProb,CaseStat=tmp,scrn=names(predByCase)[j],arm=Vacc)
		}))
	})
)

dat3$scrn = revalue(dat3$scrn,c('noTiter'='None','noPOST'='Day 1 Only','noPRE'='Week 6 Only','All'='Day 1, Week 6'))

p3 = ggplot(dat3,aes(x=CaseStat,y=PredProb)) + theme_bw() + facet_grid(scrn~arm) + geom_boxplot() + scale_y_log10() + theme(axis.title.x=element_blank()) + ylab('Predicted Probability of Zoster')
ggsave('figures/predProbs.pdf',plot=p3)

