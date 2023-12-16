GetRes1 = function(lasso_fit_lambda, embed, phecode_lst){
  res = lapply(1:length(lasso_fit_lambda), function(i){
    x = lasso_fit_lambda[[i]]
    res_part = lapply(phecode_lst, function(phecode){
      beta_list = x[[which(names(x)==phecode)]]
      tmp = lapply(1:ncol(beta_list), function(j){
        beta = beta_list[,j]
        names(beta) = rownames(beta_list)
        beta = beta[which(beta!=0)]
        beta1 = beta[which(names(beta)%in%rownames(embed))]
        if(length(beta)!=length(beta1)){
          cat("wrong!",phecode)
        }
        stopifnot(length(beta)==length(beta1))
        if(length(beta1)==0){
          y = embed[which(rownames(embed)==phecode),]
          return(norm(y,'2'))
        }else{
          if(length(beta1)==1){
            x = t(t(embed[match(names(beta1), rownames(embed)),]))
          }else{
            x = t(embed[match(names(beta1), rownames(embed)),])
          }
          beta1 = as.matrix(beta1, ncol = 1)
          y = embed[which(rownames(embed)==phecode),]
          return(norm(y-x%*%beta1,'2'))
        }
      })
      tmp = unlist(tmp)
      names(tmp) = colnames(beta_list)
      return(tmp)
    })
    names(res_part) = phecode_lst
    return(res_part)
  })
  names(res) = names(lasso_fit_lambda)
  return(res)
}

GetRes0 = function(lasso_fit_lambda, embed, phecode_lst){
  res = lapply(1:length(lasso_fit_lambda), function(i){
    x = lasso_fit_lambda[[i]]
    res_part = lapply(phecode_lst, function(phecode){
      beta_list = x[[which(names(x)==phecode)]]
      idx = which(rownames(beta_list)%in%rownames(embed))
      if(length(idx)==0){
        tmp = rep(norm(embed[which(rownames(embed)==phecode),],'2'), ncol(beta_list))
        names(tmp) = colnames(beta_list)
        return(tmp)
      }else if(length(idx)==1){
        beta_list = beta_list[idx,]
        tmp = t(t(embed[match(rownames(beta_list),rownames(embed)),]))%*%t(beta_list)
      }else{
        beta_list = beta_list[idx,]
        tmp = t(embed[match(rownames(beta_list),rownames(embed)),])%*%beta_list
      }
      tmp = tmp-embed[which(rownames(embed)==phecode),]
      tmp = apply(tmp, 2, norm, '2')
      names(tmp) = colnames(beta_list)
      return(tmp)
    })
    names(res_part) = phecode_lst
    return(res_part)
  })
  names(res) = names(lasso_fit_lambda)
  return(res)
}

GetRes = function(lasso_fit_lambda, embed, phecode_lst, cosd = NULL){
  dd = 0
  if(length(cosd)!=0){
    idx = match(rownames(cosd),rownames(embed))
    embed = embed[idx,]
    cosd = abs(cosd)
    dd = 1
  }
  res = lapply(1:length(lasso_fit_lambda), function(i){
    x = lasso_fit_lambda[[i]]
    res_part = lapply(phecode_lst, function(phecode){
      beta_list = x[[which(names(x)==phecode)]]
      idx = which(rownames(beta_list)%in%rownames(embed))
      idy = which(rownames(embed)==phecode)
      if(length(idx)==0){
        tmp = rep(norm(embed[which(rownames(embed)==phecode),],'2'), ncol(beta_list))
        names(tmp) = colnames(beta_list)
        return(tmp)
      }else if(length(idx)==1){
        beta_list = beta_list[idx,]
        idbeta = which(rownames(embed)==names(beta_list))
        if(dd == 0){
          tmp = t(t(embed[idbeta,]))%*%t(beta_list)
        }else{
          tmp = t(t(embed[idbeta,]))*cosd[idy,idbeta]
          tmp = tmp %*% t(beta_list)
        }
      }else{
        beta_list = beta_list[idx,]
        idbeta = match(rownames(beta_list),rownames(embed))
        if(dd == 0){
          tmp = t(embed[idbeta,])%*%beta_list
        }else{
          tmp = t(embed[idbeta,]*cosd[idy,idbeta])%*%beta_list
        }
      }
      tmp = tmp-embed[which(rownames(embed)==phecode),]
      tmp = apply(tmp, 2, norm, '2')
      names(tmp) = colnames(beta_list)
      return(tmp)
    })
    names(res_part) = phecode_lst
    return(res_part)
  })
  names(res) = names(lasso_fit_lambda)
  return(res)
}

load("/gpfs/alpine/proj-shared/csc489/LASSO/ElasticNet_fit_alpha_lambda_dim1000_RmFamily_cos_220819.Rdata")
load("/gpfs/alpine/proj-shared/csc489/R/embedding_1000.Rdata")
#load("/gpfs/alpine/proj-shared/csc489/R/cosine_1000.Rdata")

phecode_lst <- c("PheCode:411.4", "PheCode:555.1", "PheCode:714.1", "PheCode:555.2", 'PheCode:428', 'PheCode:250.1', 'PheCode:250.2', 'PheCode:296.2')
embed = as.matrix(embed)

trainres = GetRes(lasso_fit_lambda, embed, phecode_lst)

Nlist = lapply(1:length(lasso_fit_lambda), function(i){
  x = lasso_fit_lambda[[i]]
  N_part = lapply(phecode_lst, function(phecode){
    beta = x[[which(names(x)==phecode)]]
    N_beta = apply(beta, 2, function(x) sum(x!=0))
    names(N_beta) = colnames(beta)
    return(N_beta)
  })
  names(N_part) = phecode_lst
  return(N_part)
})
names(Nlist) = names(lasso_fit_lambda)

rm(list = setdiff(ls(),c('trainres','Nlist','phecode_lst','GetRes')))
save(trainres, Nlist, file = '/gpfs/alpine/proj-shared/csc489/LASSO/res/ElasticNet_res_RmFamily_Missing_cos_220820.Rdata')

load("/gpfs/alpine/proj-shared/csc489/R/MGB_devideby_VA_cosine_1000.Rdata")
load("/gpfs/alpine/proj-shared/csc489/LASSO/ElasticNet_fit_alpha_lambda_MGB_dim1000_RmFamily_cos_220819.Rdata")
load("/gpfs/alpine/proj-shared/csc489/R/MGB_embedding_normed_for_lasso.Rdata")
MGB_embed_map = as.matrix(MGB_embed_map)

testres = GetRes(lasso_fit_lambda_MGB, MGB_embed_map, phecode_lst, MGBdVA)

if(1==0){
  testres = lapply(1:length(lasso_fit_lambda_MGB), function(i){
    x = lasso_fit_lambda_MGB[[i]]
    cat("i =", i,"\n")
    testres_part = lapply(phecode_lst, function(phecode){
      beta_list = x[[which(names(x)==phecode)]]
      tmp = lapply(beta_list, function(beta){
        beta = beta[which(beta!=0)]
        if(length(beta)==0){
          y = MGB_embed_map[which(rownames(MGB_embed_map)==phecode),]
          return(norm(y,'2'))
        }
        beta1 = beta[which(names(beta)%in%rownames(MGB_embed_map))]
        if(length(beta)!=length(beta1)){
          cat("wrong!",phecode)
        }
        stopifnot(length(beta)==length(beta1))
        if(length(beta1)==0){
          y = MGB_embed_map[which(rownames(MGB_embed_map)==phecode),]
          return(norm(y,'2'))
        }else{
          if(length(beta1)==1){
            x = t(t(MGB_embed_map[match(names(beta1), rownames(MGB_embed_map)),]))
          }else{
            x = t(MGB_embed_map[match(names(beta1), rownames(MGB_embed_map)),])
          }
          beta1 = as.matrix(beta1, ncol = 1)
          y = MGB_embed_map[which(rownames(MGB_embed_map)==phecode),]
          return(norm(y-x%*%beta1,'2'))
        }
      })
      tmp = unlist(tmp)
      names(tmp) = names(beta_list)
      return(tmp)
    })
    names(testres_part) = phecode_lst
    return(testres_part)
  })
  names(testres) = names(lasso_fit_lambda_MGB)
}


save(trainres, testres, Nlist, file = '/gpfs/alpine/proj-shared/csc489/LASSO/res/ElasticNet_res_RmFamily_Missing_cos_220820.Rdata')

rm(list = setdiff(ls(),c('MGB_embed_map','phecode_lst','GetRes','MGBdVA')))

load('/gpfs/alpine/proj-shared/csc489/LASSO/res/ElasticNet_res_RmFamily_Missing_cos_220820.Rdata')
load("/gpfs/alpine/proj-shared/csc489/LASSO/ElasticNet_fit_alpha_lambda_dim1000_RmFamily_cos_220819.Rdata")

testres2 = GetRes(lasso_fit_lambda, MGB_embed_map, phecode_lst, MGBdVA)


readme = 'testres: test residual after handling with the missing data in MGB embedding. testres2: without handling the missing data in MGB embedding.'
save(trainres, testres, testres2, Nlist, readme, file = '/gpfs/alpine/proj-shared/csc489/LASSO/res/ElasticNet_res_RmFamily_Missing_cos_220820.Rdata')
