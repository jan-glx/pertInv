library("InvariantCausalPrediction")
library("pertInv")

load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")

YY = cbind(rowSums(count_matrix), count_matrix)


i <- 53 #"m_Nfkb1_2"
dtrain <- xgb.DMatrix(YY[isTrain,], label =  guide_matrix[isTrain, i])
xgb.cv(list(booster="gbtree",objective="binary:logistic", eta=0.02,eval_metric="logloss",
            subsample=0.5,colsample_bytree=0.8,max_depth=1), dtrain, nfold=3, nrounds=300)

xgb.cv(list(booster="gblinear",objective="binary:logistic", eta=0.01,eval_metric="logloss",lambda=1,alpha=100), dtrain, nfold=3, nrounds=200)
fitted_booster <- xgb.train(list(booster="gbtree",objective="binary:logistic", eta=0.02,eval_metric="logloss",
                                 subsample=0.8,colsample_bytree=0.8,max_depth=1,eval_metric="logloss"), dtrain, nrounds=300)
fitted_booster <- xgb.train(list(booster="gblinear",objective="binary:logistic", eta=0.02,eval_metric="logloss",
                                 subsample=0.8,colsample_bytree=0.8,max_depth=1,eval_metric="logloss"), dtrain, nrounds=300)

dt2 <- data.table(guide=i,
                  p_ko=predict(fitted_booster, YY[isTest,]),
                  guide_detected=guide_matrix[isTest, i],
                  LL_method="xgboost")
dt2[, .(x_entropy_loss=cross_entropy_lossf(p_ko,guide_detected)),
   by=.(LL_method,guide)]
