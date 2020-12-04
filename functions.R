#install.packages("ggpubr")
#remotes::install_version("cowplot", version = "0.9.2")
#devtools::install_github("wilkelab/cowplot")

library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(readr)
#library(phyloseq)
library(ggpubr)
#library(gbm)
library(lawstat)
library(Metrics)



FeatureIDset <- function(DF, FeatureID){
  DF1 <- DF[FeatureID,]
  rownames(DF1) <- FeatureID
  DF1[is.na(DF1)] <- 0
  return(DF1)
}

rmdup <- function(mt, FUN = 'sum'){
  # library(Matrix.utils)
  # #mt <- LXmt
  # c <- Matrix.utils::aggregate.Matrix(mt, row.names(mt), fun = FUN)
  # MT <- matrix(c, nrow = length(c@Dimnames[[1]]))
  # rownames(MT) <- c@Dimnames[[1]]
  # colnames(MT) <- c@Dimnames[[2]]
  # return(MT)

  if(FUN == 'sum'){
    return(t(sapply(by(mt,rownames(mt),colSums),identity)))
  }
  if(FUN == 'mean'){
    return(t(sapply(by(mt,rownames(mt),colMeans),identity)))
  }
}



RanSpl <- function(DT, label=0.5, ratio=0.3, seed=3492, labeled=F){ # column: sample ID ## labeled Ture: last row is Group info
  #DT = D_con
  #ratio = 0.3
  if(!labeled){
    set.seed(seed)
    test_i <- sample(1:ncol(DT), replace = F, size = ncol(DT)*ratio)
    train_i <- setdiff(1:ncol(DT), test_i)
    train = DT[,train_i]; train['Group',] = label
    test = DT[,test_i]; test['Group',] = label
  }else{
    #DT = umaster_mt
    grps <- names(table(as.numeric(DT[nrow(DT),])))
    dt1 <- DT[,which(DT[nrow(DT),] == grps[1])]
    dt2 <- DT[,which(DT[nrow(DT),] == grps[2])]
    
    set.seed(seed)
    test_i <- sample(1:ncol(dt1), replace = F, size = ncol(dt1)*ratio)
    train_i <- setdiff(1:ncol(dt1), test_i)
    train1 = dt1[,train_i]; test1 = dt1[,test_i]
    
    set.seed(seed)
    test_i <- sample(1:ncol(dt2), replace = F, size = ncol(dt2)*ratio)
    train_i <- setdiff(1:ncol(dt2), test_i)
    train2 = dt2[,train_i]; test2 = dt2[,test_i]
    
    train = cbind(train1, train2)
    test = cbind(test1, test2)
  }
  
  return(list(
    TRAIN = train,
    TEST = test
  ))
}


boostp <- function(DF, Xx=5){ ### row: feature(last row: group label) / column: sample ID
  #DF = conspl$TEST 
  #Xx=3
  Gpi <- DF[nrow(DF),1]
  DF <- DF[-nrow(DF),]
  btDF <- apply(DF, 1, FUN = function(r){sample(r, replace = T, size = length(r)*(Xx-1) )}) %>% t() %>% as.data.frame()
  colnames(btDF) <- paste('Bts', 1:(ncol(DF)*(Xx-1)), sep = "" )
  DF <- cbind(DF, btDF)
  DF['Group',] <- Gpi
  return(DF)
}

EnsemblValid <- function(testPredfls){
  predDataList <- lapply(testPredfls, read.csv)
  predVals <- lapply(predDataList, function(df){df$prd})
  ensVal <- do.call(args = predVals, what = cbind) %>% rowMeans()
  resTbl <- predDataList[[1]]
  resTbl$Predict_Ens <- ensVal
  resTbl$Group <- ifelse(resTbl$y==0.1, 'Control', 'Case')
  head(resTbl)
  
  XX <- Epi::ROC(resTbl$Predict_Ens, resTbl$y)
  
  ggplot(XX$res, aes(rev(1-spec), rev(sens) )) +
    theme_bw() +
    geom_line(color = '#88144F') +
    geom_abline(slope = 1, color='#69539F') +
    xlab("1 - Specificity") +
    ylab("Sensitivity") +
    annotate('text', x = 0.75, y = 0.25, label = paste('AUC:', round(XX$AUC, 3)), size = 10) +
    scale_y_continuous(breaks=seq(0.0, 1, 0.1)) +
    scale_x_continuous(breaks=seq(0.0, 1, 0.1)) -> ROC
    
  
  ggplot(resTbl, aes(Predict_Ens, 1, color = Group)) +
    geom_jitter(alpha = 0.5) +
    theme_classic() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank()) +
    xlab("Prediction Value") +
    scale_color_manual(values = c("#E12C72", "#ACDE67")) -> ScatterPlot
  
  x = list('roc'=ROC, 'scatter'=ScatterPlot)
}


LoadOTUtxt <- function(file=NA,
                       Genus=T, TaxaSelc="", Filter=1000, levelsep = " s__.*", SKip = 1, LogScale = T){
  #file = "Profiling/SR_Control/output_closed/otu_table_w.txt"; Filter = 1000
  #file = "~/HDD7T_B/MDnaviProfile_5_ALLforMODEL/BLOOD_LungCancer/Profile/output_close/otu_table_mc2_w.txt"
  #Genus = T
  #TaxaSelc = ""
  #Filter = 1000
  #levelsep = " s__.*"
  
  scon <- read_delim(file, delim = "\t", skip = SKip)
  scon <- as.data.frame(scon)
  scon_mt <- scon[,-c(1, ncol(scon))] %>% as.matrix()
  rownames(scon_mt) <- scon$taxonomy
  
  if(Genus){
    rownames(scon_mt) <- gsub(pattern = levelsep, replacement = "", rownames(scon_mt))
    scon_mt <- rmdup(scon_mt) %>% as.data.frame()
    if(length(TaxaSelc)>1){
      scon_mt <- scon_mt[TaxaSelc,]
      rownames(scon_mt) <- TaxaSelc
      scon_mt[is.na(scon_mt)] <- 0
    }
  }else{
    scon_mt <- rmdup(scon_mt) %>% as.data.frame() 
  }
  
  scon_mt <- scon_mt[,which(colSums(scon_mt) > Filter)]
  ## D_0__ 제거
  scon_mt <- scon_mt[grep(pattern = "D_0__", x = rownames(scon_mt), invert = T),]
  ## Mitochondria/Chloroplast 제거
  scon_mt <- scon_mt[grep(pattern = "Mitochondria", rownames(scon_mt), invert = T),]
  scon_mt <- scon_mt[grep(pattern = "Chloroplast", rownames(scon_mt), invert = T),]
  
  ## g__unculture/unidentified ~~ -> unidentified 통일
  reID <- gsub("g__uncultured.*", "g__unidentified", rownames(scon_mt))
  reID <- gsub("g__unident.*", "g__unidentified", reID)
  scon_mt2 <- as.matrix(scon_mt)
  rownames(scon_mt2) <- reID
  scon_mt2 <- rmdup(scon_mt2) %>% as.data.frame()
  
  if(LogScale){
    return(log2(scon_mt2+1))
  }else{
    return(scon_mt2)
  }
}


dflog2comp <- function(DF){ ## colnames: sampleID, rownames: features
  DF <- log2(DF+1)
  DF <- apply(DF, 2, function(c){ (c/sum(c))*100 }) %>% as.data.frame()
  return(DF)
}


createfdumpInput <- function(DF, label='세부진단명', extrainfo = c("성별", "나이", "출처(기관)")){
  #label='세부진단명'
  #extrainfo = c("성별", "나이", "출처(기관)")
  #DF <- sampledatalist_uc
  #DF <- COfile_UR
  DF <- DF[,c("검체번호", "시퀀싱 일자", "분석 기기", "(재)시퀀싱 일자", "(재)분석 기기", label, extrainfo)]
  #DF <- DF[,c("검체번호", "시퀀싱 일자", "분석 기기", "(재)시퀀싱 일자", "(재)분석 기기")]
  
  ID <- DF$검체번호
  Machine <- ifelse(is.na(DF$`(재)분석 기기`), DF$`분석 기기`, DF$`(재)분석 기기`)
  DATE <- ifelse(is.na(DF$`(재)시퀀싱 일자`), format(DF$`시퀀싱 일자`), format(DF$`(재)시퀀싱 일자`)) %>%
    gsub(pattern = "-", replacement = "\\.")
  Label <- DF[,label]
  
  return(cbind(ID, DATE, Machine, Label, DF[,extrainfo]))
  }




# loadOtuWinfo <- function(file, Group='None'){
#   #file = "data/Profile/"
#   #Group = "Control"
#   #file = 'data_200309/Profile_CD/Profile_S'
#   #file = "data/URINE/Profiling_LC/"
#   #Group = "Control"
#   
#   dG <- LoadOTUtxt(dir(file, recursive = T, full.names = T) %>% grep(pattern = 'otu_table_mc2_w.txt', value = T), Genus = T)
#   dS <- LoadOTUtxt(dir(file, recursive = T, full.names = T) %>% grep(pattern = 'otu_table_mc2_w.txt', value = T), Genus = F)
#   
#   meta <- read.csv(dir(file, recursive = T, pattern = "fileList.csv", full.names = T), row.names = 1) %>% as.data.frame() %>% filter(exist==1)
#   rownames(meta) <- meta$ID
#   sIDs <-  colsplit(string = colnames(dS), pattern = "-", names = c(1,2,3,4))[,1]
#   meta <- meta[sIDs,]
# 
#   ### Chao1
#   #dir(file, recursive = T, full.names = T) %>% grep(pattern = 'output_close/otu_table_mc2_w.biom', value = T)
#   
#   Biom <- phyloseq::import_biom(dir(file, recursive = T, full.names = T) %>% grep(pattern = 'otu_table_mc2_w.biom', value = T))
#   rf <- phyloseq::rarefy_even_depth(Biom, sample.size = 1000, verbose = F)
#   er <- estimate_richness(rf)
#   er$ID <- rownames(er)
#   er <- gather(er[,c(1,2,4,6,7,10)], "Method", "Value", -"ID")
#   er$Method <- factor(er$Method, levels = c("Observed", "Chao1", "Shannon", "Simpson", "ACE"))
#   er$Group <- Group
#   
#   return(list(metaInfo=meta, Genus=dG, Species=dS, Rich=er))
# }

# DiveR <- function(biomfile){
#   Biom <- phyloseq::import_biom(BIOMfilename = biomfile)
#   rf <- phyloseq::rarefy_even_depth(Biom, sample.size = 1000, verbose = F)
#   er <- estimate_richness(rf)
#   er$ID <- rownames(er)
#   er <- gather(er[,c(1,2,4,6,7,10)], "Method", "Value", -"ID")
#   er$Method <- factor(er$Method, levels = c("Observed", "Chao1", "Shannon", "Simpson", "ACE"))
#   return(er)
# }

OTUmerge <- function(df1, df2, method = 'all'){ # method = ['all', 'intersect']
  #df1 = OTUlist$CDs$Genus
  #df2 = OTUlist$CONs$Genus
  if(method == 'all'){
    IDs <- unique(c(rownames(df1), rownames(df2)))
  }
  if(method == 'intersect'){
    IDs <- intersect(rownames(df1), rownames(df2)) 
  }
  DFx <- cbind(df1[IDs,], df2[IDs,])
  rownames(DFx) <- IDs
  DFx[is.na(DFx)] <- 0
  return(DFx)
}

normComp <- function(Df, log=T){
  if(log){
    apply(log2(Df+1), 2, function(c){c/sum(c)}) %>% as.data.frame() %>% return()
  }else{
    apply(Df, 2, function(c){c/sum(c)}) %>% as.data.frame() %>% return()
  }
}

SSAcalc <- function(df, cut=0.5){ ## sensitivity, specificity, accuracy calculator # higher(0.9/1.0) is Positive
  #library(dplyr)
  #df <- read.csv("data/URINE/AI_Modeling/QT_200320_4/RandomSplit_1/RLN_itr1_AUC0.966/test_pred.csv")
  #cut=0.5
  tY <- ifelse(df$y > cut, 1, 0)
  pY <- ifelse(df$prd > cut, 1, 0)
  
  matrix(c(
    sum(tY==1 & pY==1), sum(tY==0 & pY==1),#TP / FP
    sum(tY==1 & pY==0), sum(tY==0 & pY==0) #FN / TN
  ), nrow = 2, byrow = T) %>% as.data.frame() -> cfm
  
  rownames(cfm) <- c("Positive", "Negative")
  colnames(cfm) <- c("Ture", "False")
  
  
  return(list(confusion = cfm,
              accuracy = (cfm[1,1]+cfm[2,2])/sum(cfm),
              sensitivity = sum(cfm[1,1])/sum(cfm[,1]),
              specificity = sum(cfm[2,2])/sum(cfm[,2])
              )
         )
}



genderValancer4Control <- function(control_meta, case_meta){
  #control_meta = cont_meta; case_meta = Cancer$meta
  FMr <- table(case_meta$성별)['F']/table(case_meta$성별)['M'] # 케이스 성비 f/m
  MFr <- table(case_meta$성별)['M']/table(case_meta$성별)['F'] # 케이스 성비 m/f
  cont_FMr <- table(control_meta$성별)['F']/table(control_meta$성별)['M'] # 컨트롤 성비
  wc <- ceiling(FMr*table(control_meta$성별)['M'])
  
  if(FMr < cont_FMr){ ## 남자가 더 많을때
    midx <- which(cont_meta$성별 == 'M')
    set.seed(3492)
    widx <- sample(which(cont_meta$성별 == 'F'), wc, replace = F)
  }else{ ## 여자가 더 많을때
    midx <- which(cont_meta$성별 == 'F')
    wc <- ceiling(MFr*table(control_meta$성별)['F']) # 컨트롤 여자수
    set.seed(3492)
    widx <- sample(which(cont_meta$성별 == 'M'), wc, replace = F)
  }
  
  return(c(widx, midx))
}



otuAccumulator <- function(otuDt, RelativeAdundance = T, K = 1e-7){  #K=5e-6
  #options(warn=-1)
  #otuDt = otu_acc
  #otuDt = FoodOtu
  #K = 1e-7; RelativeAdundance = T
  #rm(otuDt)
  #otuDt <- LoadOTUtxt(otu1)
  
  if(sum(rownames(otuDt) == 'Unassigned') == 1){ #denovo Unassigned 제거
    otuDt <- otuDt[which(rownames(otuDt) != 'Unassigned'),]
  }
  
  if(RelativeAdundance){
    otuDt <- otuDt %>% apply(2, function(c){c/sum(c)}) %>% as.data.frame()
  }
  
  taxalev <- colsplit(rownames(otuDt), pattern = "; ", names = c('k', 'p', 'c', 'o', 'f', 'g'))
  apply(otuDt, 2, function(c){
    #c = otuDt[,1]
    lapply(1:5, function(i){
      #i = 5
      tmpa <- data.frame(tx = as.character(taxalev[,i]))
      tmpa$tx <- as.character(tmpa$tx)
      tmpd <- data.frame(tx = as.character(taxalev[,i]), value = as.numeric(c)) %>% group_by(tx) %>% summarise(sum = sum(value), .groups = 'drop') %>% ungroup()
      tmpm <- merge(tmpa, tmpd) # %>% as.matrix; rownames(tmpm) <- as.character(tmpm[,1])
      tmpx <- matrix(tmpm$sum); rownames(tmpx) <- tmpm$tx
      return( (K)*(10^i)*tmpx[as.character(tmpa$tx),] )
    }) %>% as.data.frame() -> acmtx
    acmtx$g <- as.numeric(c)
    return(rowSums(acmtx))
  }) %>% as.data.frame() -> cumDt
  rownames(cumDt) <- rownames(otuDt)
  return(cumDt)
}


ENSvalues <- function(d1){
  # rln 
  rlnbs <- dir(d1, pattern = "RLN") %>% gsub(pattern = ".*AUC", replacement = "") %>% sort(decreasing = T); rlnbs <- rlnbs[1]
  rdd <- paste(d1, "/", dir(d1, pattern = "RLN") %>% grep(pattern = rlnbs, value = T), sep = "")[1]
  rfi <- read.csv(dir(rdd, pattern = "FeatureImportance", full.names = T), row.names = 2)
  rtp <- read.csv(dir(rdd, pattern = "test_pred.csv", full.names = T), row.names = 1)
  # gbm
  gdd <- paste(d1, "/", dir(d1, pattern = "GBM"), sep = "") #paste(d1, "/", dir(d1, pattern = "RLN") %>% grep(pattern = rlnbs, value = T), sep = "")
  gfi <- read.csv(dir(gdd, pattern = "FeatureImportance", full.names = T), row.names = 2)
  gtp <- read.csv(dir(gdd, pattern = "test_pred.csv", full.names = T), row.names = 1)
  etp <- (rtp + gtp)/2
  return(etp)
}



#Dir <- "/home/jwkang/MNT/HDD3T_JWK/Project_T/X200327_MDnavi_v3/Modeling/output//ColonCancer"
MultiTestValidation <- function(AIresDir, Threshold = 0.5, Range = 'ALL'){
  #AIresDir = "Results3_PSmatchTrial/AImodel/ColonCancer/Output_P_ColonCancer/"; Threshold = 0.5
  
  Dir <- AIresDir
  Rans <- dir(Dir, pattern = "RandomSplit")
  
  # if(exists('RES_DF')){
  #   rm(RES_DF)
  # }
  
  if(Range != 'ALL'){
    Rans <- Rans[1:Range]
  }
  
  xi = 1
  for(rd in Rans){
    #rd = Rans[1]
    print(rd)
    d1 <- paste(Dir, rd, sep = "/")
    
    # rln 
    rlnbs <- dir(d1, pattern = "RLN") %>% gsub(pattern = ".*AUC", replacement = "") %>% sort(decreasing = T); rlnbs <- rlnbs[1]
    rdd <- paste(d1, "/", dir(d1, pattern = "RLN") %>% grep(pattern = rlnbs, value = T), sep = "")[1]
    rfi <- read.csv(dir(rdd, pattern = "FeatureImportance", full.names = T), row.names = 1)
    rtp <- read.csv(dir(rdd, pattern = "test_pred.csv", full.names = T), row.names = 1)
    
    # gbm
    gdd <- paste(d1, "/", dir(d1, pattern = "GBM"), sep = "") #paste(d1, "/", dir(d1, pattern = "RLN") %>% grep(pattern = rlnbs, value = T), sep = "")
    gfi <- read.csv(dir(gdd, pattern = "FeatureImportance", full.names = T), row.names = 1)
    gtp <- read.csv(dir(gdd, pattern = "test_pred.csv", full.names = T), row.names = 1)
    
    etp <- (rtp + gtp)/2
    
    # rln <- Epi::ROC(rtp$prd, rtp$y, plot = F)$AUC
    # gbm <- Epi::ROC(gtp$prd, gtp$y, plot = F)$AUC
    # ens <- Epi::ROC(etp$prd, etp$y, plot = F)$AUC
    
    rln <- auc(ifelse(rtp$y > 0.5, 1, 0), rtp$prd)
    gbm <- auc(ifelse(gtp$y > 0.5, 1, 0), gtp$prd)
    ens <- auc(ifelse(etp$y > 0.5, 1, 0), etp$prd)
    
    ssa <- SSAcalc(rtp, cut = Threshold)
    
    
    if(xi == 1){
      RES_DF <- data.frame(AUC_ANN = rln, AUC_GBM = gbm, AUC_Ensemble = ens, 
                           Accuracy = ssa$accuracy, Sensitivity = ssa$sensitivity, specificity = ssa$specificity)
      annFI_DF <- rfi
      gbmFI_DF <- gfi
      xi = 0
    }else{
      RES_DF <- rbind(RES_DF, data.frame(AUC_ANN = rln, AUC_GBM = gbm, AUC_Ensemble = ens,
                                         Accuracy = ssa$accuracy, Sensitivity = ssa$sensitivity, specificity = ssa$specificity))
      annFI_DF <- rbind(annFI_DF, rfi)
      gbmFI_DF <- rbind(gbmFI_DF, gfi)
    }
    #data.frame(AUC_ANN = rln, AUC_GBM = gbm, AUC_Ensemble = ens)
  }
  
  gather(RES_DF[,1:3], "Method") %>% mutate(Method = factor(Method, levels = c("AUC_GBM", "AUC_ANN", "AUC_Ensemble"))) %>% 
    ggplot(aes(Method, value, color = Method)) +
    geom_boxplot() +
    ylim(c(0.6, 1)) +
    geom_jitter(alpha = 0.3, width = 0.2) +
    theme_bw() -> Img
  
  RES_DF_t <- RES_DF
  RES_DF['mean',] = colMeans(RES_DF_t)
  RES_DF['median',] = apply(RES_DF_t, 2, function(c){median(c)})
  RES_DF['sd',] = apply(RES_DF_t, 2, function(c){sd(c)})
  RES_DF['max',] = apply(RES_DF_t, 2, function(c){max(c)})
  RES_DF['min',] = apply(RES_DF_t, 2, function(c){min(c)})
  
  FI_DF <- annFI_DF %>% group_by(feature) %>% summarise(ann_mean_weight = mean(weight)) %>% ungroup() %>% arrange(feature)
  tmpA <- gbmFI_DF %>% group_by(feature) %>% summarise(gbm_mean_weight = mean(weight)) %>% ungroup() %>% arrange(feature)
  FI_DF$gbm_mean_weight <- tmpA$gbm_mean_weight
  FI_DF$mean_weight <- ((FI_DF$ann_mean_weight + FI_DF$gbm_mean_weight)/2)
  FI_DF <- FI_DF %>% arrange(desc(mean_weight))
  
  return(list(DF=RES_DF, Plot=Img, FeatureImportance = FI_DF))
}




FoldChg <- function(DF, Labeldf){
  #DF = df
  #Labeldf = GRP
  apply(DF, 1, function(r){
    g1n = names(table(Labeldf))[1]
    g2n = names(table(Labeldf))[2]
    return(log2(mean(r[Labeldf == g1n]) / mean(r[Labeldf == g2n])))
  }) -> lc2
  
  FC <- data.frame(taxa = rownames(DF), fc2 = lc2)
  return(FC)
}

BiomTTest <- function(DF, Labeldf){
  #DF = df
  #Labeldf = GRP
  apply(DF, 1, function(r){
    #r = as.numeric(DF[1,])
    ftest_pval <- var.test(r ~ Labeldf)$p.value
    eqv <- ifelse(ftest_pval > 0.05, T, F)
    x <- t.test(r ~ Labeldf, var.equal = eqv, alternative = "two.sided")
    return(t.test(r ~ Labeldf, var.equal = eqv, alternative = "two.sided")$p.value)
  }) -> Pvals
  
  g1n = names(table(Labeldf))[1]
  g2n = names(table(Labeldf))[2]
  
  lg2fc = log2(as.numeric(rowMeans(DF[,which(Labeldf == g1n)]) / rowMeans(DF[,which(Labeldf == g2n)])))
  
  PV <- data.frame(taxa = rownames(DF), 
                   ttest_pvalue = Pvals, 
                   fdr = p.adjust(Pvals, method = 'fdr'), 
                   bonferroni = p.adjust(Pvals, method = 'bonferroni')) %>% mutate(adjP001 = ifelse( (fdr < 0.05 & bonferroni < 0.05), 1, 0))
  PV[,paste('log2FC_', g1n, "_", g2n, sep = "")] = lg2fc
  return(PV)
}

CirClePlot <- function(DFfromBiomTTest, plotMain = "X"){ ## pvalue = size, logfc = color(orange:case/green:control)
  library(ggraph); library(igraph)
  #DFfromBiomTTest = DFpos
  dfx1 <- cbind(DFfromBiomTTest, colsplit(DFfromBiomTTest$taxa, pattern = "; ", names = c('k', 'p', 'c', 'o', 'f', 'g')))
  #colnames(dfx1)
  
  Edge1 <- dfx1[,c('p', 'o', 'g')]
  Edge1 <- Edge1[which(!duplicated(Edge1$o)),]
  colnames(Edge1) <- c("from", 'to', 'genus')
  
  Edge2 <- dfx1[,c('o', 'taxa', 'g')]
  colnames(Edge2) <- c("from", 'to', 'genus')
  
  EdgeRB <- rbind(Edge1, Edge2)
  
  Vetis1 <- dfx1[,c('taxa', 'g', 'bonferroni', 'UpDown')]
  colnames(Vetis1) <- c("ID", 'name', 'pval', 'fc')
  Vetis1$txtCol = 1
  
  Vetis2 <- dfx1[,c('p', 'g', 'bonferroni', 'UpDown')]
  Vetis2$g <- Vetis2$p; Vetis2$bonferroni <- 7; Vetis2$UpDown <- NA
  colnames(Vetis2) <- c("ID", 'name', 'pval', 'fc')
  Vetis2 <- Vetis2[which(!duplicated(Vetis2$ID)),]
  Vetis2$txtCol = NA
  
  Vetis3 <- dfx1[,c('o', 'g', 'bonferroni', 'UpDown')]
  Vetis3$g <- Vetis3$o; Vetis3$bonferroni <- 3; Vetis3$UpDown <- 0
  colnames(Vetis3) <- c("ID", 'name', 'pval', 'fc')
  Vetis3 <- Vetis3[which(!duplicated(Vetis3$ID)),]
  Vetis3$txtCol = 1
  
  VetisRB <- rbind(Vetis1, Vetis2, Vetis3)
  
  mygraph <- graph_from_data_frame( EdgeRB, vertices=VetisRB )
  
  
  # ggraph(mygraph, layout = 'circlepack', weight=order(pval)**2) + 
  #   geom_node_circle() +
  #   geom_node_label(aes(label=name, size = order(pval)), fill = 'white', position = "identity") +
  #   scale_size_continuous(range = c(1,5)) +
  #   theme_void()
  
  #png("PAPER_EVblood_200526/Ashma_Control_ttest.png", width = 10000, height = 10000, res = 350)
  ggraph(mygraph, layout = "nicely") +
    geom_edge_link() + 
    geom_node_point(size = 1) +
    geom_node_label(aes(label=name, size = pval, fill = fc, color = txtCol), position = "identity") +
    scale_size_continuous(range = c(1,6)) +
    scale_fill_gradient2(low = 'blue', high = 'red', midpoint = 0, na.value = '#455421') +
    scale_color_gradient2(low = 'black', mid = 'black', high = 'black', na.value = 'white') +
    ggtitle("Asthma vs Control") + 
    theme_void() -> cirPlot 
  
  return(list(CirPlot = cirPlot, TTestDF = ttDF))
}


### alpha Diversity (Richness) ###
#Biom <- phyloseq::import_biom(dir(file, recursive = T, full.names = T) %>% grep(pattern = 'otu_table_mc2_w.biom', value = T))
alpDiversity <- function(BiomMrg, SampleSize = 1000){
  
  library(phyloseq)
  #BiomFileList = list("Profiling/SR_Asthma/output_closed/otu_table_w.biom", "Profiling/SR_Control/output_closed/otu_table_w.biom")
  # BiomFileList = list("Profiling/SR_Asthma/output_closed/otu_table_w.biom",
  #                     "Profiling/SR_Control/output_closed/otu_table_w.biom",
  #                     "Profiling/SR_BileDuctCancer/output_closed/otu_table_w.biom",
  #                     "Profiling/SR_Braintumor/output_closed/otu_table_w.biom")
  # 
  # BiomList <- lapply(BiomFileList, function(Ls){phyloseq::import_biom(Ls)})
  
  # if(length(BiomList) == 1){
  #   print('Merging for 1 Biom file')
  #   BiomMrg <- merge_phyloseq(BiomList[[1]])
  # }else{
  #   print(paste('Merging for', length(BiomList), 'Biom files'))
  #   for(i in 1:length(BiomList)){
  #     if(i == 1){
  #       BiomMrg <- merge_phyloseq(BiomList[[i]], BiomList[[i+1]])
  #       print(paste('Done file:', i))
  #       print(paste('Done file:', i+1))
  #     }else{
  #       if(i != length(BiomList)){
  #         BiomMrg <- merge_phyloseq(BiomMrg, BiomList[[i+1]])
  #         print(paste('Done file:', i+1))
  #       }else{
  #         print('finish')
  #         break
  #       }
  #     }
  #   }
  # }
  
  print(paste("Sample No: ", ncol(BiomMrg@otu_table), sep = ""))
  
  rf <- phyloseq::rarefy_even_depth(BiomMrg, verbose = F, sample.size = SampleSize)
  
  er <- estimate_richness(rf)
  er$ID <- rownames(er)
  return(er)
  #er <- gather(er[,c(1,2,4,6,7,10)], "Method", "Value", -"ID")
  #er$Method <- factor(er$Method, levels = c("Observed", "Chao1", "Shannon", "Simpson", "ACE"))
  # ggplot(er[er$Method %in% c("Observed", "Chao1"),], aes(ID, Value, fill =Method))+
  #   geom_bar(stat="identity", position = 'dodge') +
  #   theme_bw() +
  #   ggtitle("Alpha-Diversity")
}
