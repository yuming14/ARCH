library(Matrix)
library(readr)
library(dplyr)
library(pROC)
library(knitr)
library(PRROC)

get_supervise_dict = function(Side_effects_1, Side_effects_01, Side_effects_02,
                              train_prop = 0.3){
  stopifnot(train_prop>0 & train_prop<1)
  druglist = table(Side_effects_1$Drug_RXNORM)
  set.seed(214)
  druglist = sample(druglist)
  drugprop = cumsum(druglist)/sum(druglist)
  idtrain = min(which(drugprop>=train_prop))
  drug_train = names(druglist)[1:idtrain]
  drug_test = names(druglist)[-c(1:idtrain)]
  dict = list(train = Side_effects_1[which(Side_effects_1$Drug_RXNORM%in%drug_train),],
              test = Side_effects_1[which(Side_effects_1$Drug_RXNORM%in%drug_test),])
  dict01 = list(train = Side_effects_01[which(Side_effects_01$Drug%in%drug_train),],
                test = Side_effects_01[which(Side_effects_01$Drug%in%drug_test),])
  dict02 = list(train = Side_effects_02[which(Side_effects_02$Drug%in%drug_train),],
                test = Side_effects_02[which(Side_effects_02$Drug%in%drug_test),])
  return(list(posi_dict = dict, nega_dict1 = dict01, nega_dict2 = dict02))
}

get_supervise_dict_unfair = function(Side_effects, train_prop = 0.3){
  stopifnot(train_prop>0 & train_prop<1)
  set.seed(214)
  n = nrow(Side_effects)
  id1 = sample(1:n, floor(n*train_prop))
  id2 = setdiff(1:n, id1)
  return(list(train = Side_effects[id1,],
              test = Side_effects[id2,]))
}


get_grad_incorrect = function(M, Other, LOINC, id, coef){
  alpha = coef[1]
  beta = coef[2]
  lambda = coef[3]
  
  S = LOINC %*% M %*% t(Other)
  
  Salpha = exp(-alpha * S)
  Salphas = matrix(0, nrow = nrow(LOINC), ncol = nrow(Other))
  ida = cbind(id$idL[which(id$ans==1)],id$idO[which(id$ans==1)])
  Salphas[ida] = Salpha[ida]
  Salphas = rowSums(Salphas) + exp(-lambda*alpha)
  rm(ida)
  
  Sbeta = exp(beta * S)
  Sbetas = matrix(0, nrow = nrow(LOINC), ncol = nrow(Other))
  idb = cbind(id$idL[which(id$ans==0)],id$idO[which(id$ans==0)])
  Sbetas[idb] = Sbeta[idb]
  Sbetas = rowSums(Sbetas) + exp(lambda*beta)
  rm(idb)
  
  gradM = lapply(1:nrow(id), function(i){
    if(id$ans[i]==1){
      coef = - Salpha[id$idL[i],id$idO[i]]/Salphas[id$idL[i]]
    }else{
      coef = Sbeta[id$idL[i],id$idO[i]]/Sbetas[id$idL[i]]
    }
    return(t(t(LOINC[id$idL[i],])) %*% t(Other[id$idO[i],]) * coef)
  })
  gradM = Reduce("+", gradM)
  return(gradM)
}

get_grad = function(M, Other, LOINC, id, coef){
  alpha = coef[1]
  beta = coef[2]
  lambda = coef[3]
  
  id$eMe = unlist(lapply(1:nrow(id), function(i){
    return(c(t(LOINC[id$idL[i],]) %*% M %*% t(t(Other[id$idO[i],]))))
  }))
  id$eMe = exp(id$eMe * c(-alpha, beta)[2-id$ans])
  idj = id %>%
    group_by(idL,ans) %>%
    summarise(sumj = sum(eMe))
  idj$sumj = idj$sumj  + exp(lambda*c(-alpha,beta)[2-idj$ans])
  
  gradM = lapply(1:nrow(idj), function(i){
    idi = which(id$idL==idj$idL[i] & id$ans==idj$ans[i])
    if(length(idi)==1){
      ej = t(t(Other[id$idO[idi],])) * id$eMe[idi]
    }else{
      ej = t(Other[id$idO[idi],]) %*% t(t(id$eMe[idi]))
    }
    gMsub = (-1)^(idj$ans[i]) * t(t(LOINC[idj$idL[i],])) %*% t(ej) / idj$sumj[i]
    return(gMsub)
  })
  
  gradM = Reduce("+", gradM)
  return(gradM)
}

get_supervied = function(Other, LOINC, dict_label, type = 1, 
                         coef = c(2,50,0.5),
                         maxstep = 100, epsilon = 1e-2, stepsize = 1e-1){
  ## obtain supervised Other_embed
  ## dict_label: three columns: 'other', 'loinc' and 'ans'
  #### ans = 1 indicates positive pairs while ans = 0 indicates negative pairs.
  ## type: 
  #### type = 0: use positive pairs only; 
  #### type = 1: use both positive pairs and negative pairs.
  ## coef = c(alpha, beta, lambda)
  ## L = 1/alpha \sum log(1+\sum exp(-alpha(S_{ij}-lambda))) + 
  #### 1/beta \sum log(1+\sum exp(beta(S_{ij}-lambda)))
  maxstep = max(1, maxstep)
  epsilon = max(epsilon, 1e-100)
  stepsize = max(stepsize, 1e-100)
  
  alpha = coef[1]
  beta = coef[2]
  lambda = coef[3]
  
  stopifnot(c("other","loinc","ans")%in%colnames(dict_label))
  dict_label = dict_label %>%
    filter(other%in%rownames(Other) & loinc%in%rownames(LOINC))
  if(type == 0){
    dict_label = dict_label[which(dict_label$ans==1),]
  }else if(type == 1){
    dict_label = dict_label[which(dict_label$ans%in%c(0,1)),]
  }
  stopifnot(nrow(dict_label)>2)
  
  olist = unique(dict_label$other)
  llist = unique(dict_label$loinc)
  Other1 = Other[match(olist,rownames(Other)),]
  LOINC = LOINC[match(llist,rownames(LOINC)),]
  id = data.frame(idO = match(dict_label$other, olist),
                  idL = match(dict_label$loinc, llist),
                  ans = dict_label$ans)
  
  M = diag(1, nrow = ncol(LOINC), ncol = ncol(Other))
  delta = numeric(maxstep)
  for(step in 1:maxstep){
    gradM = get_grad(M, Other1, LOINC, id, coef)
    M = M - stepsize * gradM
    delta[step] = norm(gradM)
    if(delta[step] < epsilon) break
  }
  delta = delta[1:step]
  if(step == maxstep) cat("not converge!","\n")
  newOther = Other %*% t(M)
  return(list(Other = newOther, M=M, delta=delta))
}


Procrustes <- function(X1,X2){
  #Return Omega = arg min||X1 - X2 Omgea||_F
  H = t(X2)%*%X1
  mod = svd(H)
  return(mod$u%*%t(mod$v))
}

get_supervise_embedding = function(embed_Drug, embed_side, label){
  idx = data.frame(drug = c(label$Drug_RXNORM, label$Drug_CUI,
                            label$Drug_RXNORM, label$Drug_CUI),
                   side = c(label$Side_PheCode, label$Side_PheCode,
                            label$Side_CUI, label$Side_CUI))
  idx = na.omit(idx)
  idx$drug = match(idx$drug, rownames(embed_Drug))
  idx$side = match(idx$side, rownames(embed_side))
  idx = na.omit(idx)
  idx = idx[!duplicated(idx),]
  Omega = Procrustes(embed_side[idx$side,],embed_Drug[idx$drug,])
  embed_Drug = embed_Drug%*%Omega
  return(embed_side)
}

get_AUC = function(embed_drug, embed_side, 
                   Drug_Side_effects_1, Drug_Side_effects_01,
                   align = TRUE){
  drug_code_list = unique(c(Side_effects_1$Drug_RXNORM,
                            Side_effects_1$Drug_CUI,
                            Side_effects_01$Drug))
  drug_code_list = intersect(drug_code_list, rownames(embed_drug))
  side_code_list = unique(c(Side_effects_1$Side_PheCode,
                            Side_effects_1$Side_CUI,
                            Side_effects_01$Side))
  side_code_list = intersect(side_code_list, rownames(embed_side))
  cos = embed_drug[match(drug_code_list,rownames(embed_drug)),] %*% 
    t(embed_side[match(side_code_list,rownames(embed_side)),])
  
  idx_posi = data.frame(drug1 = match(Drug_Side_effects_1$Drug_RXNORM,drug_code_list),
                        drug2 = match(Drug_Side_effects_1$Drug_CUI,drug_code_list),
                        side1 = match(Drug_Side_effects_1$Side_PheCode,side_code_list),
                        side2 = match(Drug_Side_effects_1$Side_CUI,side_code_list))
  idx_posi = as.matrix(idx_posi)
  idx_posi = cbind(cos[idx_posi[,c(1,3)]],
                   cos[idx_posi[,c(1,4)]],
                   cos[idx_posi[,c(2,3)]],
                   cos[idx_posi[,c(2,4)]])
  idx_posi = apply(idx_posi, 1, max, na.rm = TRUE)
  idx_posi = idx_posi[which(!is.na(idx_posi))]
  # idx_posi = idx_posi[which(!is.infinite(idx_posi))]
  
  idx_nega = cbind(drug = match(Drug_Side_effects_01$Drug, drug_code_list),
                   side = match(Drug_Side_effects_01$Side, side_code_list))
  idx_nega = cos[idx_nega]
  idx_nega = idx_nega[which(!is.na(idx_nega))]
  
  if(align){
    set.seed(214)
    n = min(length(idx_posi), length(idx_nega))
    idx_posi = sample(idx_posi, n)
    idx_nega = sample(idx_nega, n)
  }
  
  pr<-pr.curve(scores.class1 = idx_nega, scores.class0 = idx_posi)
  roc<-roc.curve(scores.class1 = idx_nega, scores.class0 = idx_posi)
  ans = c(ROC = roc$auc, PRC1 = pr$auc.integral, PRC2 = pr$auc.davis.goadrich,
          n_posi = length(idx_posi), n_nega = length(idx_nega))
  return(ans)
}

get_dict_to_label = function(Side_effects_1, Side_effects_0, 
                             drug_code_list, side_code_list){
  dict_posi = data.frame(other = c(Side_effects_1$Drug_RXNORM,Side_effects_1$Drug_RXNORM,
                                   Side_effects_1$Drug_CUI,Side_effects_1$Drug_CUI),
                         loinc = c(Side_effects_1$Side_PheCode,Side_effects_1$Side_CUI,
                                   Side_effects_1$Side_PheCode,Side_effects_1$Side_CUI),
                         ans = 1)
  dict_posi = dict_posi %>%
    filter(other%in%drug_code_list, loinc%in%side_code_list)
  
  dict_nega = data.frame(other = Side_effects_0$Drug,
                         loinc = Side_effects_0$Side,
                         ans = 0)
  dict_nega = dict_nega %>%
    filter(other%in%drug_code_list, loinc%in%side_code_list)
  
  return(rbind(dict_posi, dict_nega))
}
