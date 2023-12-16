library(Matrix)
library(glmnet)
load("/gpfs/alpine/proj-shared/csc489/R/neg_log10.p.adjusted_1500.Rdata")
load("/gpfs/alpine/proj-shared/csc489/R/cosine_1000.Rdata")
sum(rownames(neg_log10.p.adjusted)==rownames(cosine_similarity))
idx = grep("PheCode:",rownames(neg_log10.p.adjusted))
length(idx)
neg_log10.p.adjusted = neg_log10.p.adjusted[idx,]
#PheName = rownames(neg_log10.p.adjusted)
PheName = c("PheCode:411.4", "PheCode:555.1", "PheCode:714.1", "PheCode:555.2", 'PheCode:428', 'PheCode:250.1', 'PheCode:250.2', 'PheCode:296.2')
head(PheName)
load("/gpfs/alpine/proj-shared/csc489/R/embedding_1000.Rdata")
embed = embed[which(!is.na(rownames(embed))),]
sum(rownames(embed)==colnames(neg_log10.p.adjusted))
sum(rownames(embed)==rownames(cosine_similarity))
dim(neg_log10.p.adjusted)
ls()

rm_family = function(codeselect, phecode){
  idx = grep("PheCode:", codeselect)
  if(length(idx)==0){
    return(codeselect)
  }else{
    idx1 = grep(phecode, codeselect)
    idx2 = which(sapply(codeselect, function(code) grepl(code, phecode)))
    idx = union(idx1, idx2)
    if(length(idx)==0){
      return(codeselect)
    }else{
      return(codeselect[setdiff(1:length(codeselect),idx)])
    }
  }
}

t1 = Sys.time()
lambda_lst = list(c(seq(1,51,5)*1e-6,seq(60,1000,50)*1e-6),
                  c(seq(1,30,2)*1e-6,seq(35,350,20)*1e-6),
                  c(seq(1,100,4)*1e-7,seq(1,5,0.5)*1e-5),
                  c(seq(1,20,1)*1e-7,seq(21,40,2)*1e-7,seq(5,10,2)*1e-6))
lambda_lst = lapply(lambda_lst, sort, decreasing = TRUE)
alpha_lst = seq(0.25,1,0.25)
names(lambda_lst) = alpha_lst
lasso_fit_lambda = list()
length(lasso_fit_lambda) = length(alpha_lst)
names(lasso_fit_lambda) = alpha_lst
for(j in 1:length(lambda_lst)){
  best_lambda = lambda_lst[[j]]
  lasso_phe_lst = lapply(PheName, function(phecode){
    idx = which(rownames(embed)==phecode)
    idselect = which(rownames(neg_log10.p.adjusted)==phecode)
    i = idselect
    if(i%%10==0) cat(idselect,", ")
    if(i%%100==0) cat("\n")
    idselect = which(neg_log10.p.adjusted[idselect,]!=0)
    idselect = setdiff(idselect, idx)
    codeselect = colnames(neg_log10.p.adjusted)[idselect]
    codeselect = rm_family(codeselect, phecode)
    idselect = which(colnames(neg_log10.p.adjusted)%in%codeselect)
    if(length(idselect)==0){
      return(numeric())
    }
    x = t(embed[idselect,]*abs(cosine_similarity[idx,idselect]))
    x = x[,order(cosine_similarity[idx,idselect], decreasing = TRUE)]
    best_model <- glmnet(x = x, 
                         y = embed[idx,], alpha = alpha_lst[j], lambda = best_lambda,
                         intercept = FALSE,
                         standardize = FALSE)
    tmp = coef(best_model)
    colnames(tmp) = best_lambda
    tmp = tmp*abs(cosine_similarity[idx,match(rownames(tmp),colnames(cosine_similarity))])
    return(tmp)
  })
  names(lasso_phe_lst) = PheName
  lasso_fit_lambda[[j]] = lasso_phe_lst
  names(lasso_fit_lambda) = alpha_lst
  save(lasso_fit_lambda, file = "/gpfs/alpine/proj-shared/csc489/LASSO/ElasticNet_fit_alpha_lambda_dim1000_RmFamily_cos_220828.Rdata")
  #change the file here
}
