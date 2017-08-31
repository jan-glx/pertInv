library("pertInv")

covariates.dt<-fread("results/covariates.dt.csv")
load(file="results/count_matrix.RData")

Y <- stabilize_Anscombes(count_matrix) #count_matrix
X <- model.matrix(~I(CDR^2)+I(CDR^3)+CDR+total_counts_scaled+batch,data=covariates.dt) #+guide

fit <- lm(Y ~ .,as.data.table(X))
Y_adj <- residuals(fit)

# development_mode <- TRUE
subset <- seq_len(if(is_dev()) 8 else ncol(Y_adj))

glmboost_selection <- data.table(target_gene = colnames(Y_adj)[subset], i=subset)
glmboost_selection <-
  glmboost_selection[
    ,{cat(".");.(predictor_gene = names(coef(suppressWarnings(mboost::glmboost(Y_adj[,-i[1]],Y_adj[,i[1]], center = TRUE)))))
    },by=target_gene]


corr_pval <- psych::corr.test(Y_adj,Y_adj[,subset,drop=F],adjust="none")$p

corr_pval <- melt(data.table(corr_pval, keep.rownames=TRUE), variable.factor = FALSE, id.vars = "rn")
setnames(corr_pval, c("predictor_gene", "target_gene", "p.value"))
corr_pval <- corr_pval[`target_gene` != `predictor_gene`]

corr_pval[,p.val.adj := p.adjust(p.value, metho="BH")]

corr_pval[,selected:=FALSE]
corr_pval[glmboost_selection, on=c("target_gene", "predictor_gene"), selected:= TRUE]

significant_at_FDR.05_threshold <- corr_pval[p.val.adj<0.05,max(p.value)]



figure("Comparison of feature selection with glmboost vs cor.test",
ggplot(corr_pval, aes(x=p.value, fill=factor(selected,levels=c(TRUE,FALSE))))+geom_histogram(binwidth=0.1,boundary=0)+
  scale_fill_discrete(name="Selected by glmboost",
                          breaks=c(TRUE, FALSE),
                          labels=c("yes", "no"))+xlab("correlation test p-value")
)

setorder(corr_pval,p.value)
corr_pval[,fraction_selected := cumsum(selected)]

library(viridis)
figure("glmboost slection vs. correllation test",
ggplot(corr_pval, aes(x=order(p.val.adj),color=-log(p.val.adj+.Machine$double.eps),y=fraction_selected))+geom_line(size=2)+
  xlab("rank correlation p.value")+ylab("cumulative selected by glmboost method")+scale_color_viridis(name="neg. log \ncorrellation p.value\n(BH adjusted)")+
  theme(legend.position = c(0.8, 0.2))
)

