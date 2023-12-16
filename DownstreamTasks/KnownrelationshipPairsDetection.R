get_type = function(AllRelationPairs){
  pairs = lapply(AllRelationPairs, function(x){ 
    x = x[,c('type',"pair","RELA","id")]
    x = x[!duplicated(x),]
  })
  pairs = do.call("rbind",pairs)
  pairs = pairs[!duplicated(pairs),]
  tn = lapply(1:28, function(i){
    if(sum(pairs$id==i)==0) return(c(type = "NA", name = "NA"))
    if(i<=9){
      name = unique(pairs$pair[which(pairs$id==i)])
      type = ifelse(i<9,unique(pairs$type[which(pairs$id==i)]),"rm")
      name = paste(name,"(",ifelse(i<=3,"sim",ifelse(i<=8,"rela","")),")",sep="")
    }else{
      name = ifelse(i<=27, pairs$RELA[which(pairs$id==i)], "CUI-code")
      type = unique(pairs$type[which(pairs$id==i)])
      name = paste(name,"(",ifelse(type=="similarity","sim","rela"),")",sep="")
    }
    stopifnot(length(type)==1&length(name)==1)
    return(c(type = type, name = name))
  })
  tn = do.call("rbind", tn)
  return(list(type = tn[,1], name = tn[,2]))
}


get_power_pmitest_sub = function(log_p_matrix, cosine_similarity, pairs, target_FDR){
  target_FDR = sort(unique(target_FDR), decreasing = TRUE)
  if(nrow(pairs)==0) return(list(ans = numeric(length(target_FDR)*3+2), 
                                 ptable = NULL))
  idx = which(pairs$code1%in%rownames(log_p_matrix) &
                pairs$code2%in%rownames(log_p_matrix))
  if(length(idx)==0) return(list(ans = numeric(length(target_FDR)*3+2), 
                                 ptable = NULL))
  pairs = pairs[idx,]
  grouplist = data.frame(g1 = dict$group[match(pairs$code1, dict$code)],
                         g2 = dict$group[match(pairs$code2, dict$code)])
  group0 = grouplist[!duplicated(grouplist),]
  group = data.frame(g1 = apply(group0, 1, min),
                     g2 = apply(group0, 1, max))
  rm(group0)
  group = group[!duplicated(group),]
  group = na.omit(group)
  codegroup = dict$group[match(rownames(log_p_matrix),dict$code)]
  pvalue_list = lapply(1:nrow(group), function(i){
    if(group[i,1]==group[i,2]){
      idx = which(codegroup==group[i,1])
      idlist = as.data.frame(t(combn(idx, 2)))
      colnames(idlist) = c("id1", "id2")
      idgrouplist = which(grouplist$g1==group[i,1] & grouplist$g2==group[i,2])
    }else{
      idx1 = which(codegroup==group[i,1])
      idx2 = which(codegroup==group[i,2])
      idlist = data.frame(id1 = rep(idx1, each = length(idx2)),
                          id2 = rep(idx2, length(idx1)))
      idgrouplist = which(grouplist$g1==group[i,1] & grouplist$g2==group[i,2])
      idgrouplist = unique(c(idgrouplist,
                             which(grouplist$g1==group[i,2] & grouplist$g2==group[i,1])))
    }
    grouplist_i = data.frame(id1 = match(pairs$code1[idgrouplist], rownames(log_p_matrix)),
                             id2 = match(pairs$code2[idgrouplist], rownames(log_p_matrix)))
    grouplist_i = data.frame(id1 = c(grouplist_i$id1, grouplist_i$id2),
                             id2 = c(grouplist_i$id2, grouplist_i$id1))
    idlist0 = setdiff(idlist, grouplist_i)
    idlist0$relation = 0
    idlist1 = intersect(idlist, grouplist_i)
    idlist1$relation = 1
    idlist = rbind(idlist1, idlist0)
    rm(idlist0, idlist1)
    #grouplist_i$relation = 1
    #idlist = merge(idlist, grouplist_i, all.x = TRUE)
    idlist$logp = log_p_matrix[cbind(idlist$id1, idlist$id2)]
    return(idlist)
  })
  pvalue_list = do.call("rbind", pvalue_list)
  pvalue_list = pvalue_list[order(pvalue_list$logp),]
  pvalue_list$fdr = 1
  nmax = n = nrow(pvalue_list)
  ans = numeric()
  NRelation = sum(pvalue_list$relation)
  order_cos = order(cosine_similarity[cbind(pvalue_list$id1,pvalue_list$id2)], 
                    decreasing = TRUE)
  for(fdr in target_FDR){
    nmax = max(c(which(pvalue_list$logp[1:nmax]<= log(fdr) + log(c(1:nmax) / n / (log(n)+1)))),0)
    pvalue_list$fdr[1:nmax] = fdr
    ans = c(ans, 
            sum(pvalue_list$relation[1:nmax])/NRelation, 
            sum(pvalue_list$relation[order_cos[1:nmax]])/NRelation,
            nmax)
  }
  ans = c(ans, NRelation, n)
  names(ans) = c(paste(rep(c("power_test", "power_cos", "NSel"), length(target_FDR)),
                       rep(target_FDR, each = 3)),
                 "NPairs", "NHypo")
  return(list(ans = ans, ptable = pvalue_list))
}


get_power_compare_sub = function(cosine_similarity, N, pairs){
  if(nrow(pairs)==0) return(list(ans = numeric(length(N)*2+2), 
                                 ptable = NULL))
  idx = which(pairs$code1%in%rownames(cosine_similarity) &
                pairs$code2%in%rownames(cosine_similarity))
  if(length(idx)==0) return(list(ans = numeric(length(N)*2+2), 
                                 ptable = NULL))
  pairs = pairs[idx,]
  grouplist = data.frame(g1 = dict$group[match(pairs$code1, dict$code)],
                         g2 = dict$group[match(pairs$code2, dict$code)])
  group0 = grouplist[!duplicated(grouplist),]
  group = data.frame(g1 = apply(group0, 1, min),
                     g2 = apply(group0, 1, max))
  rm(group0)
  group = group[!duplicated(group),]
  group = na.omit(group)
  codegroup = dict$group[match(rownames(cosine_similarity),dict$code)]
  pvalue_list = lapply(1:nrow(group), function(i){
    if(group[i,1]==group[i,2]){
      idx = which(codegroup==group[i,1])
      idlist = as.data.frame(t(combn(idx, 2)))
      colnames(idlist) = c("id1", "id2")
      idgrouplist = which(grouplist$g1==group[i,1] & grouplist$g2==group[i,2])
    }else{
      idx1 = which(codegroup==group[i,1])
      idx2 = which(codegroup==group[i,2])
      idlist = data.frame(id1 = rep(idx1, each = length(idx2)),
                          id2 = rep(idx2, length(idx1)))
      idgrouplist = which(grouplist$g1==group[i,1] & grouplist$g2==group[i,2])
      idgrouplist = unique(c(idgrouplist,
                             which(grouplist$g1==group[i,2] & grouplist$g2==group[i,1])))
    }
    grouplist_i = data.frame(id1 = match(pairs$code1[idgrouplist], rownames(cosine_similarity)),
                             id2 = match(pairs$code2[idgrouplist], rownames(cosine_similarity)))
    grouplist_i = data.frame(id1 = c(grouplist_i$id1, grouplist_i$id2),
                             id2 = c(grouplist_i$id2, grouplist_i$id1))
    idlist0 = setdiff(idlist, grouplist_i)
    idlist0$relation = 0
    idlist1 = intersect(idlist, grouplist_i)
    idlist1$relation = 1
    idlist = rbind(idlist1, idlist0)
    rm(idlist0, idlist1)
    #grouplist_i$relation = 1
    #idlist = merge(idlist, grouplist_i, all.x = TRUE)
    idlist$cos = cosine_similarity[cbind(idlist$id1,idlist$id2)]
    return(idlist)
  })
  pvalue_list = do.call("rbind", pvalue_list)
  pvalue_list = pvalue_list[order(pvalue_list$cos,decreasing = TRUE),]
  NRelation = sum(pvalue_list$relation)
  ans = numeric()
  for(n in N){
    ans = c(ans,
            sum(pvalue_list$relation[1:n])/NRelation,
            n)
  }
  ans = c(ans, NRelation, nrow(pvalue_list))
  names(ans) = c(paste(rep(c("power", "NSel"), length(N)),
                       rep(N, each = 2)),
                 "NPairs", "NHypo")
  return(list(ans = ans, ptable = pvalue_list))
}

Get_Power_PMItest = function(log_p_matrix, cosine_similarity, AllRelationPairs, target_FDR, 
                             echo = TRUE, runall = FALSE){
  tn = get_type(AllRelationPairs)
  fulltable = list()
  length(fulltable) = 3
  if(runall){
    for(k in 3:1){
      anslist = list()
      length(anslist) = 28
      for(i in 1:28){
        if(echo) cat("pair",k,"type",i,"\n")
        pairs = AllRelationPairs[[k]]
        pairs = pairs[which(pairs$id==i),]
        ans = get_power_pmitest_sub(log_p_matrix, cosine_similarity, pairs, target_FDR)
        if(echo) cat(ans$ans,"\n")
        anslist[[i]] = ans$ans
      }
      #ptable = lapply(anslist, function(x) x$ptable)
      #ans = do.call("rbind", lapply(anslist, function(x) x$ans))
      ans = do.call("rbind", anslist)
      ans = as.data.frame(ans)
      rownames(ans) = tn$name
      pairtype = tn$type
      id = which(pairtype=="similarity")
      ans = rbind(ans, apply(ans[id,],2,weighted.mean,ans[id,ncol(ans)-1]))
      ans[nrow(ans),ncol(ans)-1] = sum(ans[id,ncol(ans)-1])
      ans[nrow(ans),ncol(ans)] = sum(ans[id,ncol(ans)])
      rownames(ans)[nrow(ans)] = "weighted.sim"
      id = which(pairtype=="related")
      ans = rbind(ans, apply(ans[id,],2,weighted.mean,ans[id,ncol(ans)-1]))
      ans[nrow(ans),ncol(ans)-1] = sum(ans[id,ncol(ans)-1])
      ans[nrow(ans),ncol(ans)] = sum(ans[id,ncol(ans)])
      rownames(ans)[nrow(ans)] = "weighted.rela"
      ans = ans[c(which(pairtype=="similarity"),
                  which(pairtype=="related"),
                  which(rownames(ans)=="weighted.sim"),
                  which(rownames(ans)=="weighted.rela")),]
      #ans[,3*(1:length(target_FDR))] = apply(ans[,3*(1:length(target_FDR))], 2, round)
      if(echo) kable(ans,"latex",3)
      #return(list(ans=ans, ptable = ptable))
      fulltable[[k]] = ans
    }
    names(fulltable) = c("CUI-CUI","CUI-codi","codi-codi")
    return(fulltable)
  }else{
    anslist = list()
    length(anslist) = 28
    for(i in 1:28){
      k = ifelse(i<=9, 3, ifelse(i==28, 2, 1))
      if(echo) cat("pair",k,"type",i,"\n")
      pairs = AllRelationPairs[[k]]
      pairs = pairs[which(pairs$id==i),]
      ans = get_power_pmitest_sub(log_p_matrix, cosine_similarity, pairs, target_FDR)
      if(echo) cat(ans$ans,"\n")
      anslist[[i]] = ans$ans
    }
    ans = do.call("rbind", anslist)
    ans = as.data.frame(ans)
    rownames(ans) = tn$name
    pairtype = tn$type
    id = which(pairtype=="similarity")
    ans = rbind(ans, apply(ans[id,],2,weighted.mean,ans[id,ncol(ans)-1]))
    ans[nrow(ans),ncol(ans)-1] = sum(ans[id,ncol(ans)-1])
    ans[nrow(ans),ncol(ans)] = sum(ans[id,ncol(ans)])
    rownames(ans)[nrow(ans)] = "weighted.sim"
    id = which(pairtype=="related")
    ans = rbind(ans, apply(ans[id,],2,weighted.mean,ans[id,ncol(ans)-1]))
    ans[nrow(ans),ncol(ans)-1] = sum(ans[id,ncol(ans)-1])
    ans[nrow(ans),ncol(ans)] = sum(ans[id,ncol(ans)])
    rownames(ans)[nrow(ans)] = "weighted.rela"
    ans = ans[c(which(pairtype=="similarity"),
                which(pairtype=="related"),
                which(rownames(ans)=="weighted.sim"),
                which(rownames(ans)=="weighted.rela")),]
    if(echo) kable(ans,"latex",3)
    return(ans)
  }
}


Get_Power_compare = function(cosine_similarity, refans, x = "comp", AllRelationPairs, 
                             echo = TRUE, runall = FALSE){
  tn = get_type(AllRelationPairs)
  fulltable = list()
  length(fulltable) = 3
  if(runall){
    for(k in 3:1){
      anslist = list()
      length(anslist) = 28
      for(i in 1:28){
        if(echo) cat("pair",k,"type",i,"\n")
        pairs = AllRelationPairs[[k]]
        pairs = pairs[which(pairs$id==i),]
        idsel = grep("NSel",colnames(refans[[k]]))
        N = unlist(refans[[k]][which(rownames(refans[[k]])==tn$name[i]),idsel])
        if(length(N)==0) N = rep(0, length(idsel))
        ans = get_power_compare_sub(cosine_similarity, N, pairs)
        if(echo) cat(ans$ans,"\n")
        anslist[[i]] = ans$ans
      }
      #ptable = lapply(anslist, function(x) x$ptable)
      #ans = do.call("rbind", lapply(anslist, function(x) x$ans))
      ans = do.call("rbind", anslist)
      ans = as.data.frame(ans)
      rownames(ans) = tn$name
      pairtype = tn$type
      id = which(pairtype=="similarity")
      ans = rbind(ans, apply(ans[id,],2,weighted.mean,ans[id,ncol(ans)-1]))
      ans[nrow(ans),ncol(ans)-1] = sum(ans[id,ncol(ans)-1])
      ans[nrow(ans),ncol(ans)] = sum(ans[id,ncol(ans)])
      rownames(ans)[nrow(ans)] = "weighted.sim"
      id = which(pairtype=="related")
      ans = rbind(ans, apply(ans[id,],2,weighted.mean,ans[id,ncol(ans)-1]))
      ans[nrow(ans),ncol(ans)-1] = sum(ans[id,ncol(ans)-1])
      ans[nrow(ans),ncol(ans)] = sum(ans[id,ncol(ans)])
      rownames(ans)[nrow(ans)] = "weighted.rela"
      ans = ans[c(which(pairtype=="similarity"),
                  which(pairtype=="related"),
                  which(rownames(ans)=="weighted.sim"),
                  which(rownames(ans)=="weighted.rela")),]
      #ans[,3*(1:length(target_FDR))] = apply(ans[,3*(1:length(target_FDR))], 2, round)
      if(echo) kable(ans,"latex",3)
      #return(list(ans=ans, ptable = ptable))
      fulltable[[k]] = ans
    }
    names(fulltable) = c("CUI-CUI","CUI-codi","codi-codi")
    return(fulltable)
  }else{
    anslist = list()
    length(anslist) = 28
    for(i in 1:28){
      k = ifelse(i<=9, 3, ifelse(i==28, 2, 1))
      if(echo) cat("pair",k,"type",i,"\n")
      pairs = AllRelationPairs[[k]]
      pairs = pairs[which(pairs$id==i),]
      idsel = grep("NSel",colnames(refans))
      N = unlist(refans[which(rownames(refans)==tn$name[i]),idsel])
      if(length(N)==0) N = rep(0, length(idsel))
      ans = get_power_compare_sub(cosine_similarity, N, pairs)
      if(echo) cat(ans$ans,"\n")
      anslist[[i]] = ans$ans
    }
    ans = do.call("rbind", anslist)
    ans = as.data.frame(ans)
    rownames(ans) = tn$name
    pairtype = tn$type
    id = which(pairtype=="similarity")
    ans = rbind(ans, apply(ans[id,],2,weighted.mean,ans[id,ncol(ans)-1]))
    ans[nrow(ans),ncol(ans)-1] = sum(ans[id,ncol(ans)-1])
    ans[nrow(ans),ncol(ans)] = sum(ans[id,ncol(ans)])
    rownames(ans)[nrow(ans)] = "weighted.sim"
    id = which(pairtype=="related")
    ans = rbind(ans, apply(ans[id,],2,weighted.mean,ans[id,ncol(ans)-1]))
    ans[nrow(ans),ncol(ans)-1] = sum(ans[id,ncol(ans)-1])
    ans[nrow(ans),ncol(ans)] = sum(ans[id,ncol(ans)])
    rownames(ans)[nrow(ans)] = "weighted.rela"
    ans = ans[c(which(pairtype=="similarity"),
                which(pairtype=="related"),
                which(rownames(ans)=="weighted.sim"),
                which(rownames(ans)=="weighted.rela")),]
    if(echo) kable(ans,"latex",3)
    return(ans)
  }
}




