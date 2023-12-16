library(ggplot2)

load('/gpfs/alpine/proj-shared/csc489/LASSO/res/ElasticNet_res_RmFamily_Missing_cos_220820.Rdata')
load('/gpfs/alpine/proj-shared/csc489/R/dict.Rdata')
load('testres1_0826.RData')

pdf('plot/res_alpha_1.pdf')
i = which(names(trainres)=="0.25")
for(j in 1:length(trainres[[i]])){
  code = names(trainres[[i]])[j]
  lambdalist = as.numeric(names(trainres[[i]][[j]]))
  trainres_tmp = log(trainres[[i]][[j]])
  testres_tmp = log(testres[[i]][[j]])
  testres1_tmp = log(testres1[[i]][[j]])
  testres2_tmp = log(testres2[[i]][[j]])
  N_tmp = Nlist[[i]][[j]]
  aic = trainres_tmp + N_tmp/1000
  bic = trainres_tmp + N_tmp*log(1000)/2000
  df = data.frame(lambda = rep(lambdalist, 6),
                  logres = c(testres_tmp, testres1_tmp, testres2_tmp, 
                             testres_tmp+aic, testres1_tmp+aic, testres2_tmp+aic),
                  type = rep(c("test(MapToAll)","test1(MapToNoneZero)","test2(NHM)",
                               "test+AIC","test1+AIC","test2+AIC"),
                             each=length(N_tmp)))
  dt = data.frame(x = lambdalist, 
                  y = rep(c(-1.9,-2.1),length(N_tmp))[1:length(N_tmp)], 
                  text = as.character(N_tmp))
  #p1 = ggplot(data = df, mapping = aes(x = lambda, y = logres, group = type, colour = type)) + 
  #  geom_line() + 
  #  geom_point() + 
  #  geom_hline(aes(yintercept=min(bic)), colour="#A3A500", linetype="dashed") +
  #  geom_vline(aes(xintercept=lambdalist[which.min(bic)]), colour="#A3A500", linetype="dashed") +
  #  annotate('text', x = dt$x, y = dt$y, label = dt$text, size = 2) +
  #  ggtitle(paste(code,": ",dict$desc[which(dict$code==code)]))
  p1 = ggplot(data = df, mapping = aes(x = log(lambda), y = logres, 
                                       group = type, colour = type, 
                                       linetype = type, shape = type)) + 
    geom_line() + 
    geom_point() + 
    scale_linetype_manual(values = c(2,1,2,1,2,1)) +
    scale_shape_manual(values = c(39,16,39,16,39,16)) +
    scale_color_manual(values = c("#C65146","#EC6A5C",
                                  "#4FB0C6","#4F86C6",
                                  "#8FBC94","#548687")) +
    geom_hline(aes(yintercept=min(testres_tmp+aic)), 
               colour="#EC6A5C", linetype="dotted") +
    geom_vline(aes(xintercept=log(lambdalist[which.min(testres_tmp+aic)])), 
               colour="#EC6A5C", linetype="dotted") +
    geom_hline(aes(yintercept=min(testres1_tmp+aic)), 
               colour="#4F86C6", linetype="dotted") +
    geom_vline(aes(xintercept=log(lambdalist[which.min(testres1_tmp+aic)])), 
               colour="#4F86C6", linetype="dotted") +
    geom_hline(aes(yintercept=min(testres2_tmp+aic)), 
               colour="#548687", linetype="dotted") +
    geom_vline(aes(xintercept=log(lambdalist[which.min(testres2_tmp+aic)])), 
               colour="#548687", linetype="dotted") +
    annotate('text', x = log(dt$x), y = dt$y, label = dt$text, size = 1.75) +
    ggtitle(paste(code,": ",dict$desc[which(dict$code==code)]))
  plot(p1)
}
dev.off()


load('/gpfs/alpine/proj-shared/csc489/LASSO/res/ElasticNet_res_RmFamily_Missing_cos_220819.Rdata')
load('/gpfs/alpine/proj-shared/csc489/R/dict.Rdata')
pdf('plot/res_alpha_0.25.pdf')
i = which(names(trainres)=="0.25")
for(j in 1:length(trainres[[i]])){
  code = names(trainres[[i]])[j]
  lambdalist = as.numeric(names(trainres[[i]][[j]]))
  trainres_tmp = log(trainres[[i]][[j]])
  testres_tmp = log(testres[[i]][[j]])
  testres2_tmp = log(testres2[[i]][[j]])
  N_tmp = Nlist[[i]][[j]]
  aic = trainres_tmp + N_tmp/1000
  bic = trainres_tmp + N_tmp*log(1000)/2000
  df = data.frame(lambda = rep(lambdalist, 6),
                  logres = c(trainres_tmp, testres_tmp, testres2_tmp, aic, bic,testres_tmp+aic),
                  type = rep(c('train',"test(HMiss)","test(NHMiss)","AIC","BIC","AIC+test"),
                             each=length(N_tmp)))
  dt = data.frame(x = lambdalist, 
                  y = rep(c(-1.9,-2.1),length(N_tmp))[1:length(N_tmp)], 
                  text = as.character(N_tmp))
  #p1 = ggplot(data = df, mapping = aes(x = lambda, y = logres, group = type, colour = type)) + 
  #  geom_line() + 
  #  geom_point() + 
  #  geom_hline(aes(yintercept=min(bic)), colour="#A3A500", linetype="dashed") +
  #  geom_vline(aes(xintercept=lambdalist[which.min(bic)]), colour="#A3A500", linetype="dashed") +
  #  annotate('text', x = dt$x, y = dt$y, label = dt$text, size = 2) +
  #  ggtitle(paste(code,": ",dict$desc[which(dict$code==code)]))
  p1 = ggplot(data = df, mapping = aes(x = log(lambda), y = logres, group = type, colour = type)) + 
    geom_line() + 
    geom_point() + 
    geom_hline(aes(yintercept=min(testres_tmp+aic)), 
               colour="#A3A500", linetype="dashed") +
    geom_vline(aes(xintercept=log(lambdalist[which.min(testres_tmp+aic)])), 
               colour="#A3A500", linetype="dashed") +
    annotate('text', x = log(dt$x), y = dt$y, label = dt$text, size = 1.75) +
    ggtitle(paste(code,": ",dict$desc[which(dict$code==code)]))
  plot(p1)
}
dev.off()



load('EightLasso_res_RmFamily_Missing.Rdata')
load('/gpfs/alpine/proj-shared/csc489/R/dict.Rdata')

pdf('res.pdf')

phecode_lst <- c("PheCode:411.4", "PheCode:555.1", "PheCode:714.1", "PheCode:555.2", 'PheCode:428', 'PheCode:250.1', 'PheCode:250.2', 'PheCode:296.2')
par(mar = c(7,4,6,6))

for(phecode in phecode_lst){
  idx = which(phecode_lst==phecode)
  mymain = paste(phecode,": ",dict$desc[which(dict$code==phecode)])
  plot(as.numeric(rownames(trainres)),
       trainres[,idx],
       type ="b", xlab = "lambda", ylab = "Train Residual", pch = 3,
       ylim = c(min(trainres[,idx])-0.05, max(trainres[,idx])+0.15), 
       main = mymain, lwd = 2, col = "black")
  text(rownames(Nlist), rep(c(0.08,0.12), ceiling(nrow(Nlist)/2))[1:nrow(Nlist)], 
       as.character(Nlist[,idx]), 
       cex=0.75, pos = 1, col="blue")
  par(new = TRUE)
  tmp = max(testres[,idx])-min(testres[,idx])
  tmp = mean(c(tmp,0.1))
  plot(as.numeric(rownames(testres)), testres[,idx], 
       lwd = 2, col = "red", pch = 3,
       ylim = c(min(testres[,idx])-tmp*0.4, max(testres[,idx])+tmp*0.1),
       type ="b", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "")
  axis(side = 4)
  mtext("Test Residual", side = 4, line = 3)
  legend("bottomright",col = c('black','red'), lwd = 2, legend = c('train','test'))
}
dev.off()




res.list = list(trainres = trainres, testres = testres, Nlist = Nlist)
trainres = trainres[-(1:5),]
testres = testres[-(1:5),]
Nlist = Nlist[-(1:5),]

pdf('res2.pdf')
par(mar = c(7,4,6,6))

for(phecode in phecode_lst){
  idx = which(phecode_lst==phecode)
  mymain = paste(phecode,": ",dict$desc[which(dict$code==phecode)])
  plot(rep(as.numeric(rownames(trainres)),2),
       c(trainres[,idx], testres[,idx]),
       col = c(rep('black',nrow(trainres)),rep('red',nrow(testres))),
       type = 'p', ylim = c(0.1,0.5),
       xlab = "lambda", ylab = "Residual", lwd = 2,
       pch = 3, main = mymain)
  abline(h = min(testres[,idx]), col = "red", lwd = 0.7, lty = 2)
  abline(v = as.numeric(rownames(testres)[which.min(testres[,idx])]), 
         col = "red", lwd = 0.7, lty = 2)
  text(rownames(Nlist), 
       rep(c(0.115,0.13), ceiling(nrow(Nlist)/2))[1:nrow(Nlist)], 
       as.character(Nlist[,idx]), 
       cex=0.75, pos = 1, col="blue")
  legend(0.0013,0.25,col = c('black','red'), pch = 3, lwd = 2,
         legend = c('train','test'))
}
dev.off()