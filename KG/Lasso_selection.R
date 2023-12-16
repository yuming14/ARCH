load("/gpfs/alpine/proj-shared/csc489/R/MGB_embedding_mapped_for_lasso.Rdata")
phecode_lst <- c("PheCode:411.4", "PheCode:555.1", "PheCode:714.1", "PheCode:555.2", 'PheCode:428', 'PheCode:250.1', 'PheCode:250.2', 'PheCode:296.2')
MGB_embed_map = as.matrix(MGB_embed_map)
MGB_norm = apply(MGB_embed_map[match(phecode_lst, rownames(MGB_embed_map)),], 1 , norm, '2')

#PheCode:411.4 PheCode:555.1 PheCode:714.1 PheCode:555.2   PheCode:428 
#0.5919238     0.7910263     0.5764494     0.7870555     0.6030546 
#PheCode:250.1 PheCode:250.2 PheCode:296.2 
#0.5746983     0.5154348     0.4983847 

load("/gpfs/alpine/proj-shared/csc489/LASSO/lasso_fit_lambda_dim1000_RmFamily.Rdata")
phecode_lst <- c("PheCode:411.4", "PheCode:555.1", "PheCode:714.1", "PheCode:555.2", 'PheCode:428', 'PheCode:250.1', 'PheCode:250.2', 'PheCode:296.2')

Nlist2 = lapply(1:length(lasso_fit_lambda), function(i){
  x = lasso_fit_lambda[[i]]
  N_part = lapply(phecode_lst, function(phecode){
    beta = x[[which(names(x)==phecode)]]
    return(length(beta)-1)
  })
  N_part = unlist(N_part)
  names(N_part) = phecode_lst
  return(N_part)
})
Nlist = do.call('rbind', Nlist)
rownames(Nlist) = names(lasso_fit_lambda)
colnames(Nlist) = phecode_lst

load('EightLasso_reslist_RmFamily.Rdata')