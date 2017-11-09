 
# This application Is to perform Meta-Analysis of gene expression data.
# In the initial phase we deal with Affymetrix Codelink and Illumina data to get 
# the genes strongly associated with the disease condition. The genes selected 
# here are going to analysed according to their functions and pathway analysis in order to 
# find the potential cancer biomarkers. 
 
 
cat("\014")  

# Packages requried
library(shiny)
library(markdown)
library(plyr)
library(data.table)
library(affy)
library(dplyr)
#library(codelink)
library(simpleaffy)
library(MetaDE)
library(lumi)
#library(h20kcod.db)
#library(h10kcod.db)
##library(hwgcod.db)
#library(limma)
library(MetaDE)
 

## Defining the size of file to be accepted. Here it can accept any size.
options(shiny.maxRequestSize= -1) 


## Function for different meta-analysis methods
Stouffer.test <- function(p, w) { # p is a vector of p-values
  if (missing(w)) {
    w <- rep(1, length(p))/length(p)
    print(w)
  } else {
    if (length(w) != length(p))
      stop("Length of p and w must equal!")
  }
  Zi <- qnorm(1-p) 
  Z  <- sum(w*Zi)/sqrt(sum(w^2))
  p.val <- 1-pnorm(Z)
  return(c( p.value = p.val))
}

Fisher.test <- function(p) {
  Xsq <- -2*sum(log(p))
  p.val <- pchisq(Xsq, df = 2*length(p), lower.tail = FALSE)
  return(c( p.value = p.val)) #
}


sim.AW.Fisher <- function(dff) {
  #print(df)
  p.sort <- sort(dff)
  fisher.stat <- -2 * cumsum(log(p.sort))
  cat("Calculating weights and fisher's statistics\n")
  fisher.p <- pchisq(fisher.stat, df = 2*length(dff), lower.tail = FALSE)
  cat("Final Fishers value\n")
  #stat <- apply(fisher.p, 1, min)
  stat <- min(fisher.p)
  #sum.weight <- apply(fisher.p, 1, which.min)
  cat("Final Fishers value\n")
  #stat <- apply(fisher.p, 1, min)
  return(stat)
}

 



## Function for affymetric processed data
affydata <- list()
affymetrix <- function(df,files,df1){
   
  if(length(df)<2)
    stop("Less than two files provided, please select 2 or more files to perform meta-analysis")
  
  for(i in 1:files)
  { 
    incProgress(0.1, detail = paste("File",i)) 
    file <- df[[i]]
     
    dataAnot <- merge(df1,file, by= c("ID"),all.y=TRUE)
    cat("Annotation done...",i,"\n")
    
    dataAnot <- na.omit(dataAnot)
    
    dataAnot <- dataAnot[,-c(1,2)]
    
    pvalue <- grepl("PVALUE",colnames(dataAnot))
    
    if(any(pvalue)){
       
      dt <- data.table(dataAnot)
      setkey(dt, PVALUE) ; indx <- dt[,.I[1L], by = SYMBOL]$V1 
      dataWoDuplic <- dt[indx]
    }else{
      
      cat("Trying to remove the duplicates\n")
      dataWoDuplic <- setDT(dataAnot)[order(factor(CALL, levels=c('P', 'A', 'M'))),
                .SD[1L], by = SYMBOL]
      
      print(head(dataWoDuplic,20))
    }
    
    affydata[[length(affydata)+1]] <- dataWoDuplic
  }
  
  
  affyfinal <- lapply(seq_along(affydata), function(i) setnames(affydata[[i]],
                                                                2:ncol(affydata[[i]]), paste0(names(affydata[[i]])[-1],i)))
  
  
  affyfinal <- Reduce(function(x, y) merge(x, y, by =c("SYMBOL"), all=TRUE), affyfinal)
  print(affyfinal)
  return(affyfinal)
}

###############################  End  ############################################## 


## Function to process Affymetrix raw  data.
processed <- list()
processed_data <- function(affy.data){
  if(length(affy.data)<2)
    stop("Less than two files provided, please select 2 or more files to perform meta-analysis")
  
  for(i in 1: length(affy.data))
  {
    incProgress(0.1, detail = paste("File",i))  
    eset.mas5 = mas5(affy.data[,i])
    
    ## getting the expression matrix (probesets/genes in rows, chips in columns).
    exprSet.nologs = exprs(eset.mas5)
    
    #print(head(exprSet.nologs))
    
    ## At this time let's log-transform the expression values to get a more normal distribution. 
    ## We have to remember we've done this when we calculate ratios. Logarithms can use any
    ## base, but base 2 is easiest when transforming ratios, since transformed 2-fold
    ## ratios up or down will be +1 or -1. As a result, we'll do all logs with base
    ## 2 to keep thing simplest.
    
    
    ## While we're doing Affymetrix-specific preprocessing, let's calculate an Absent/Present call for each probeset.
    # Run the Affy A/P call algorithm on the CEL files we processed above
    data.mas5calls = mas5calls(affy.data[,i])
    
    
    data.mas5calls.calls = exprs(data.mas5calls)
     
    pvalue <- assayData(data.mas5calls)[["se.exprs"]]
    
    data.full <- cbind(exprSet.nologs,data.mas5calls.calls,pvalue)
    
    data.processed <- data.frame(data.full)
    
    data.processed <- add_rownames(data.processed, "VALUE")
    
    setnames(data.processed,c("ID","INTENSITY","CALL","PVALUE")) 
    
    processed[[length(processed)+1]] <- data.processed
    
    cat("Affymetrix processing done...\n")
  }
  return(processed)
}

####################################  End  ###############################################


## Function for raw codelink data
codelink_processed <- list()

codelink_data <- function(codelinkpath){
  
  if(length(codelinkpath)<2)
    stop("Less than two files provided, please select 2 or more files to perform meta-analysis")
  
  for(i in 1:length(codelinkpath))
  {
    incProgress(0.1, detail = paste("File",i))  
    
    codset = readCodelinkSet(filename = codelinkpath[[i]])
    
    #print(head(codset))
    
    features <- readCodelink(codelinkpath[[i]])  ## An addition to get the ids
    
    ids <- features$id
    
    if(all(!is.na(ids))){
      
      codset = codCorrect(codset, method = "half", offset = 0)
      
      codset = codNormalize(codset, method = "loess", weights = getWeight(codset), loess.method = "fast")
      
      exprs <- exprs(codset)
      print(head(exprs))
      
      snr  <- getSNR(codset)
      print(head(snr))
      
      data <- data.frame(cbind(ids,exprs,snr),stringsAsFactors=FALSE)
      
      setnames(data,c("ID","SIGNALINTENSITY","SNR"))
      
      codelink_processed[[length(codelink_processed)+1]] <- data
      
      cat("Codelink processing done...\n") 
    }
  }
  if(length(codelink_processed)<2)
    stop("Either none or less than two studies has id column present.
          Please provide data with ids present and try again.")
  
  return(codelink_processed)
}


## Function for codelink processed data.
codelinkdata <- list()

codelink <- function(df,files,df1){
  
  if(length(df)<2)
    stop("Less than two files provided, please select 2 or more files to perform meta-analysis")
  
  for(i in 1:files)
  {
    incProgress(0.1, detail = paste("File",i))
    
    file <- df[[i]]
    
    #setnames(file,c("ID","SIGNALINTENSITY","SNR"))
     
    dataAnot <- merge(df1,file, by= c("ID"),all.y=TRUE)
    
    dataAnot <- na.omit(dataAnot)
    
    dataAnot <- dataAnot[,-c(1,2)]    
    #print(head(dataAnot),10)
    
    dataAnot[with(dataAnot, order(SNR,decreasing=TRUE)),]
    
    dt <- data.table(dataAnot)
    
    #print(head(dt),10)
    
    ## Removing duplicate genes, set SNR as key column to remove the duplicates
    ## here max need to be kept, so use order() to sort the data frame and then removing the duplicates.
     
    indx <- dt[,.I[1L], by = SYMBOL]$V1 
    
    dataWoDuplic <- dt[indx]    
    #print(head(dataWoDuplic))
    
    codelinkdata[[length(codelinkdata)+1]] <- dataWoDuplic
  }
  ##Naming the columns uniquely and merging the list of data frames
  cdlnkdata <- lapply(seq_along(codelinkdata), function(i) setnames(codelinkdata[[i]],
                                                                    2:ncol(codelinkdata[[i]]), paste0(names(codelinkdata[[i]])[-1],i)))
  
  cdlnkdata <- Reduce(function(x, y) merge(x, y, by =c("SYMBOL"), all=TRUE), cdlnkdata)
  
  return(cdlnkdata)
}



######################################3  End  ######################################

## Function for Differential Expression of affymetrix data (Fold Change)
foldchange <- list()
fold.change <- function(lst,newpath,df1,fold,pval,newnames,conditions){
  
  if(length(lst)<2)
    stop("Less than two files provided, please select 2 or more files to perform differential expression")
  
  ext <- grepl("*.txt$", lst)
  for(i in 1:length(lst))
  { 
    if(ext[i])
    {
      incProgress(0.1, detail = paste("File",i))  
      cond <- conditions
      cat("................................\n")
      
      ## Reading the phenotype files to construct the model matrix 
      lines <- readLines(newnames[i])
      
      pdata <- read.table(text=sub('.*(#.*)', '\\1', lines),
                         check.names=FALSE, stringsAsFactors=FALSE, 
                         comment.char='#',header=TRUE)
      #print(pdata)
      
      raw.data <-  read.affy(lst[i], path=newpath)
      
      x.rma <- call.exprs(raw.data,"rma")
         
      unique <- unique(pdata[,2])
      f <- factor(pdata[,2],levels=unique) 
      
      combination <- grep('(.*)-\\1',unique(as.vector(outer(unique,unique,FUN=paste,sep='-'))), 
                          value=TRUE,invert=TRUE)
      
      cat("possible combination:\n")
      #print(combination)
      
      if(!all(conditions[[1]]%in%combination))
       stop("Choosen contrasts is not a possible combination, please choose valid contrasts and try again")
         
      design <- model.matrix(~0+f)
      colnames(design) <- unique  
      fit <- lmFit(x.rma, design)
      
      contrast.matrix <- makeContrasts(contrasts=conditions[[1]], levels=design)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      fit2
      
      for(i in 1:ncol(contrast.matrix))
      {
        top <- topTable(fit2, coef=i, n=Inf, p.value=pval,lfc=fold, adjust.method="none")
         
        if(nrow(top)>0)
        {
          foldPval <- top[,c(1,4)]
          foldPval <- add_rownames(foldPval)
          
          setnames(foldPval,c("ID","AlogFC","PVALUE"))
          foldchange[[length(foldchange)+1]] <- foldPval
        }
      }
    }
  }
  if (length(foldchange)<2)
    stop("No gene found significant at the threshold p-value= ", pval , " and foldchange = ",fc, 
         " please change your pvalue and foldchange and try again")
  
  return(foldchange)
}

## Function to annotate for fold change, linked to above code
foldchange1 <- list()
fold.change1 <- function(df,df1){
  for(i in 1:length(df))
  {
    incProgress(0.1, detail = paste("File",i))  
    file <- df[i]
    afterAnot <- merge(df1,file, by="ID", all.y=TRUE)
         
    afterAnot <- afterAnot[,-c(1,2)]
    afterAnot <- na.omit(afterAnot)
    #print(head(afterAnot))
    
    dt <- data.table(afterAnot)
    #print(dim(dt))
    
    setkey(dt, PVALUE) ; indx <- dt[,.I[1L], by = SYMBOL]$V1 
    dataWoDuplic <- dt[indx]
    #print(head(dataWoDuplic))
    
    foldchange1[[length(foldchange1)+1]] <- dataWoDuplic
  }
  
  foldchange <- lapply(seq_along(foldchange1), function(i) setnames(foldchange1[[i]],
                                                               2:ncol(foldchange1[[i]]), paste0(names(foldchange1[[i]])[-1],i)))
  
  foldchange<- Reduce(function(x, y) merge(x, y, by =c("SYMBOL"), all=TRUE), foldchange)
  #print(head(foldchange))
  
  return(foldchange)       
}

###########################  End ######################################


## Function for codlink differential expression
coddiff <- list()

codelinkdiff <- function(lst, newpath,cntrsts,pval,fc,annotation,newnames){
  
  ext <- grepl("*.txt$", lst)
  
  if(length(lst)<2)
    stop("Less than two files provided, please select 2 or more files with at least two 
         conditions to perform differential expression")
  
  if(!any(ext))
    stop("No target files found, please check target files and try again")
  
  if(all(ext))
    stop("No data files found, please check input data files and try again")
    
    for(i in 1:length(lst))
    { 
     if(ext[i])
      {
        cat("Differential expression function called \n")
        incProgress(0.1, detail = paste("File",i))  
        
       ## Reading phenotype data to construct model matrix        
        lines <- readLines(newnames[i])
        pdata <- read.table(text=sub('.*(#.*)', '\\1', lines),
                            check.names=FALSE, stringsAsFactors=FALSE, 
                            comment.char='#',header=TRUE)
        #print(pdata)
       
        unique <- unique(pdata[,2])
        f <- factor(pdata[,2],levels=unique)
       
        cat("Treatment conditions:",unique,"\n")
        
       combination <- grep('(.*)-\\1',unique(as.vector(outer(unique,unique,FUN=paste,sep='-'))), 
                           value=TRUE,invert=TRUE)
       #print(combination)
       
       if(!all(cntrsts[[1]]%in%combination))
           stop("Choosen contrast is not a possible combination, please choose the possible combination and try again")
         
        pdata$FileName <- file.path(newpath, pdata$FileName)
        cat("setting path done..\n")
        #print(pdata)
         
        codset <- readCodelinkSet(filename = pdata$FileName) 
        codset = codCorrect(codset, method = "half", offset = 0)
        codset = codNormalize(codset, method = "loess", weights = getWeight(codset), loess.method = "fast")
    
        design <- model.matrix(~0+f)
        cat("Experimental design:\n",design)
        colnames(design) <- unique ## c("C","TT","D","DT")
        fit <- lmFit(codset, design,weights = getWeight(codset))
        contrast.matrix <- makeContrasts(contrasts=cntrsts[[1]],levels=design)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2)
        
        for(i in 1:ncol(contrast.matrix))
        {
          top <- topTable(fit2, coef=i, n=Inf, p.value=pval,lfc=fc, adjust.method="none")
           
          if(nrow(top)>0)
          {
            foldPval <- top[,c("probeName","logFC","P.Value")]
            cat("Fold change:\n")
        
            setnames(foldPval,c("PROBEID","ClogFC","P.Value"))
            #print(head(foldPval))
            
            ## Extracting information from ids to have a match with accnum 
            ids <- foldPval$PROBEID
            test <- grepl("_",ids[1])
            newids <- gsub("\\..*|_PROBE.*", "", ids)
            test1 <- grepl("GE",newids[1])
        
            if(test)
            {
              foldPval <- cbind(newids,foldPval[,-1])
              colnames(foldPval)[1] <- "PROBEID"
            }
            if(annotation=="hwgcod")
            {
            if(test)
              keys <- select(hwgcod.db, newids, c("SYMBOL"), "ACCNUM")
          
            else if(test1)
              keys <- select(hwgcod.db,newids,c("SYMBOL"),"PROBEID")
            }
        
            else if(annotation=="h20kcod")
            {
            if(test)
              keys <- select(h20kcod.db, newids, c("SYMBOL"), "ACCNUM") 
          
            else if(test1)
              keys <- select(h20kcod.db,newids,c("SYMBOL"),"PROBEID")
            }
        
            else if(annotation=="h10kcod")
            {
            if(test)
              keys <- select(h10kcod.db, newids, c("SYMBOL"), "ACCNUM")  
          
            else if(test1)
              keys <- select(h10kcod.db,newids,c("SYMBOL"),"PROBEID")
            }
        
        setnames(keys,c("PROBEID","SYMBOL"))
        cat("keys form database:")
        
        #print(head(keys))
        
        fold <- merge(keys, foldPval, by= c("PROBEID"), all.x=TRUE)
        
        fold1 <- fold[,-1]
        fold1 <- na.omit(fold1)
        
        cat("fold is:\n")
        #print(head(fold))
        #print(head(fold1))
        
        dt <- data.table(fold1)
        setkey(dt, P.Value) ; indx <- dt[,.I[1L], by = SYMBOL]$V1 
        dataWoDuplic <- dt[indx]
        coddiff[[length(coddiff)+1]] <- dataWoDuplic
      }
    }   
   }
  }
  if (length(coddiff)==0)
  stop("No gene found significant at the threshold p-value= ", pval , " and foldchange = ",fc, 
       " please change your pvalue and foldchange and try again")
  
  if (length(coddiff)==1)
    stop("Files selected belong to only one study, please provide atleast two studies to perform meta-analysis of differential expression")

  return(coddiff)
}

codlinkdiff1 <- function(df){
   
  codfinal <- lapply(seq_along(df), function(i) setnames(df[[i]],
                                                              2:ncol(df[[i]]), paste0(names(df[[i]])[-1],i)))
 
  codfinal <- Reduce(function(x, y) merge(x, y, by =c("SYMBOL"), all=TRUE), codfinal)
  #print(codfinal)
  codfinal
}

############################## End ###################################

## Function for illumina processed data to perform meta-analysis
illuminadata <- list()
illumina <- function(df,files,df1,n){
  
  if(length(df)<2)
    stop("Less than two files provided, please select 2 or more files to perform meta-analysis")
  
  for(i in 1:files)
  {
    cat("Illumina in annotation section:\n")
    incProgress(0.1, detail = paste("File",i))  
    
    file <- df[[i]]
    #print(head(file))
    
    #setnames(file,c("ID","INTENSITY","PVALUE"))
    cat("setting the column names done...\n")
     
    dataAnot <- merge(df1,file, by= c("ID"),all.y=TRUE)
    dataAnot <- na.omit(dataAnot)
    #print(head(dataAnot))
    
    dataAnot <- dataAnot[,-c(1,2)]
    dt <- data.table(dataAnot)
    #print(head(dt))
    
    setkey(dt, PVALUE) ; indx <- dt[,.I[1L], by = SYMBOL]$V1 
    dataWoDuplic <- dt[indx]
    illuminadata[[length(illuminadata)+1]] <- dataWoDuplic
  }
  
  ilmnfinal <- lapply(seq_along(illuminadata), function(i) setnames(illuminadata[[i]],
                                                                    2:ncol(illuminadata[[i]]), paste0(names(illuminadata[[i]])[-1],i)))
  
  ilmnfinal <- Reduce(function(x, y) merge(x, y, by =c("SYMBOL"), all=TRUE), ilmnfinal)
  #print(head(ilmnfinal))
  return(ilmnfinal)
}
  
########################### End ################################################


## Function to process Illumina raw data
ilm <- list()
illuminaRaw <- function(illumina){
  for(i in 1:length(illumina))
  {
    incProgress(0.1, detail = paste("File",i))  
    x.lumi <- lumiR(illumina[i])
    lumiExpr <- lumiExpresso(x.lumi, bg.correct = TRUE,normalize = TRUE,verbose = TRUE)
    
    exprs <- exprs(lumiExpr)
    pvalue <- detection(lumiExpr)
    #print(head(exprs))
    #print(head(pvalue))
    
    if(length(exprs)<3)
      {
      lumidata <- cbind(exprs,pvalue)
      data.processed <- data.frame(lumidata)
      data.processed <- add_rownames(data.processed, "VALUE")
      
      colnames(data.processed) <- c("ID","INTENSITY","PVALUE")
      ilm[[length(ilm)+1]] <- data.processed
      }else{
      exprs <- data.frame(exprs)
      #print(head(exprs))
      
      pvalue <- data.frame(pvalue)
      exprs <- add_rownames(exprs,"ID")
      pvalue <- add_rownames(pvalue,"ID") 
      #print(head(pvalue))
      
      names(exprs)[2:ncol(exprs)] <- paste("INTENSITY", 1:ncol(exprs), sep="")
      names(pvalue)[2:ncol(pvalue)] <- paste("PVALUE", 1:ncol(pvalue), sep="")
      length <- ncol(exprs)-1
      
      for(i in 1:length)
      {
        exprs1 <- exprs[,c(1,i+1)]
        pvalue1 <- pvalue[,c(1,i+1)]
        exp.pval <- merge(exprs1,pvalue1, by="ID")
        setnames(exp.pval,c("ID","INTENSITY","PVALUE"))
        ilm[[length(ilm)+1]] <- exp.pval
      }
    }
  }
  cat("Illumina raw data processing done ... ")
  return(ilm)
}

################################# End ######################################


## Function for differential expression of illumina data
ilmndiff <- list()
illuminaDe <- function(ilmNames,ilmTargets,x,fc,pval,newnames){
  
  if(length(ilmNames) && length(ilmTargets)<2)
    stop("Less than two files provided, please select 2 or more files with at least two 
         conditions to perform differential expression")
  
  if(length(ilmTargets)==0)
    stop("No target files found, please check target files and try again")
  
  if(length(ilmNames)==0)
    stop("No data files found, please check input data files and try again")
  
  for(i in 1:length(ilmNames))
  {
    incProgress(0.1, detail = paste("File",i))  
    
    lines <- readLines(ilmTargets[[i]])
    pdata <- read.table(text=sub('.*(#.*)', '\\1', lines),check.names=FALSE, 
                        stringsAsFactors=FALSE,comment.char='#',header=TRUE)
    
    unique <- unique(pdata[,2])
    f <- factor(pdata[,2],levels=unique)
    target <- readTargets(ilmTargets[[i]])
    #cat("Targets:\n")
    ##print(target)
     
    combination <- grep('(.*)-\\1',unique(as.vector(outer(unique,unique,FUN=paste,sep='-'))), 
                        value=TRUE,invert=TRUE)
    #print(combination)
    
    if(!all(x[[1]]%in%combination))
      stop("Choosen contrast is not a possible combination in target file", i , ",please choose the possible combination and try again")
    
    x.lumi <- lumiR(ilmNames[[i]])
    lumiExpr <- lumiExpresso(x.lumi, bg.correct = TRUE,normalize = TRUE,verbose = TRUE)
    
    dataMatrix <- exprs(lumiExpr)
         
    presentCount <- detectionCall(lumiExpr)
        
    selDataMatrix <- dataMatrix[presentCount > 0,]
        
    probeList <- rownames(selDataMatrix)
         
    design <- model.matrix(~0+f)
    colnames(design) <- unique ## c("C","TT","D","DT")
    
    fit <- lmFit(selDataMatrix, design)
    cat("fitting the model done...\n")
    #print(fit)
    contrast.matrix <- makeContrasts(contrasts=x[[1]],levels=design)
    
    cat("contrast matrix construction done:...\n")
    #print(contrast.matrix)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    #print(fit2)
    
    cat("Design contrast done...\n")
    #print(ncol(contrast.matrix))
    for(i in 1:ncol(contrast.matrix))
    {
      top <- topTable(fit2, coef=i, n=Inf, p.value=pval,lfc=fc, adjust.method="none")
      #print(top)
      #top <- topTable(fit2, coef=i,number=geneNum, adjust.method="BH")
      if(nrow(top)>0)
      {
        foldPval <- top[,c(1,4)]
        foldPval <- add_rownames(foldPval)
        setnames(foldPval,c("ID","IlogFC","PVALUE"))
        ilmndiff[[length(ilmndiff)+1]] <- foldPval
        #print(ilmndiff)
      }
    } 
  }
  if (length(ilmndiff)==0)
    stop("No gene found significant at the threshold p-value= ", pval , " and foldchange = ",fc, 
         " please change your pvalue and foldchange and try again")
  
  if (length(ilmndiff)==1)
    stop("Files selected belong to only one study, please provide atleast two studies to perform meta-analysis of differential expression")
  
  cat("ilmdiff is:\n")
  print(ilmndiff)
  return(ilmndiff)
}


## Function to annotate for fold change illumina data, linked to above code
foldchange2 <- list()
fold.change2 <- function(df,df1){
  for(i in 1:length(df))
  {
    incProgress(0.1, detail = paste("File",i))  
    file <- df[i]
    afterAnot <- merge(df1,file, by="ID", all.y=TRUE)
    
    afterAnot <- afterAnot[,-c(1,2)]
    afterAnot <- na.omit(afterAnot)
    #print(head(afterAnot))
    
    dt <- data.table(afterAnot)
    #print(dim(dt))
    
    setkey(dt, PVALUE) ; indx <- dt[,.I[1L], by = SYMBOL]$V1 
    dataWoDuplic <- dt[indx]
    #print(head(dataWoDuplic))
    
    foldchange1[[length(foldchange1)+1]] <- dataWoDuplic
  }
  
  foldchange <- lapply(seq_along(foldchange1), function(i) setnames(foldchange1[[i]],
                                                                    2:ncol(foldchange1[[i]]), paste0(names(foldchange1[[i]])[-1],i)))
  
  foldchange<- Reduce(function(x, y) merge(x, y, by =c("SYMBOL"), all=TRUE), foldchange)
  #print(head(foldchange))
  return(foldchange)       
}


########################################## End ########################################

## Illumina function to split the list of files as Names and Targets
names <- list()
targets <- list()
ilmTest <- function(lst){
  for(i in 1:length(lst))
  {
  
    lines <- readLines(lst[i],n=20)
    file <- read.table(text=sub('.*(#.*)', '\\1', lines),
                       check.names=FALSE, stringsAsFactors=FALSE, comment.char='#',fill=TRUE,header=TRUE)
    
    #print(head(file))
    cols1 <- colnames(file)
    fileTst <- grepl("bead", cols1,ignore.case=TRUE)
    fileTst1 <- grepl("AVG_Signal",cols1,ignore.case=TRUE)
     
    cat("Testing for target file :\n")
    
    if(any(fileTst,fileTst1))
      names[[length(names)+1]] <- lst[i]
    else
      targets[[length(targets)+1]] <- lst[i]
  }
  newList <- list("names"=names,"targets"=targets)
  return(newList)
}

############################ End #####################################

## Function for codelink data to reorder the dataframes(Automatic 
## Detection of data platform)
fun1 <- function(df){
  #print(df)
  if(ncol(df)>=2){
    if(ncol(df)==4)
    {
      df <- df[,-1]
      cat("after removing columns\n")
      #print(head(df))
    }
    indx <- sapply(df, is.character)
    df[indx] <- lapply(df[indx], type.convert)
    means <- colMeans(df)
    df <- df[order(means,decreasing=TRUE)]
    cat("data frame ordering done ...\n")
     
    if(ncol(df)==3)
      setnames(df,c("ID","SIGNALINTENSITY","SNR"))
    
    if(ncol(df)==2)
      setnames(df,c("ID","SNR"))
    
    cat("names setting done ... \n")
    df
    #print(head(df))
  }
}


## Function for Affymetrix data to reorder(Automatic Detection of data platform)
fun2 <- function(data){
  
  if(ncol(data)>4)
    stop("Data is not of correct formate, please check the files and try again")
  
  dim(data)
  indx <- sapply(data, is.character)
  data[indx] <- lapply(data[indx], type.convert)
  d1 <- data[1, , drop = FALSE]
  
  nums <- d1[, nn <- sapply(d1, is.numeric),drop=FALSE]
  #print(nums)
  
  if(length(nums)==2){
    p <- names(nums[, nums < 1, drop = FALSE])
    val <- setdiff(names(nums), p)
  }
  
  ch <- d1[, !nn, drop = FALSE]
  id <- names(ch[, grepl('_at$', as.character(unlist(ch))), drop = FALSE])
  
  if(ncol(data)==4){
    abs <- setdiff(names(ch), id)
    
    d <- data[, c(id, val, abs, p)]  
    setnames(d,c("ID","INTENSITY","CALL","PVALUE"))
  }
  else if(ncol(data)==3)
    {
    if(length(nums)==2)
      {
      p <- names(nums[, nums < 1, drop = FALSE])
      d <- data[, c(id, val, p)]
      setnames(d,c("ID","INTENSITY","PVALUE"))
    }else{
      abs <- setdiff(names(ch), id)
      d <- data[,c(id,names(nums),abs)]
      setnames(d,c("ID","INTENSITY","CALL"))
    }
  }
  else if(ncol(data)==2)
    {
    p <- names(nums)
    d <- data[,c(id,p)]
    setnames(d,c("ID","PVALUE"))
  }
   
  return(d)
}


## Function to order illumina data
ilmnOrder <- function(data){
  dim(data)
  indx <- sapply(data, is.character)
  data[indx] <- lapply(data[indx], type.convert)
  d1 <- data[1, , drop = FALSE]
  
  nums <- d1[, nn <- sapply(d1, is.numeric),drop=FALSE]
  #print(nums)
  if(length(nums)==2){
    p <- names(nums[, nums <= 1, drop = FALSE])
    
    if(length(p)==0)
      stop("Data is not of illumina formate, please check the data files and try again")
    
    val <- setdiff(names(nums), p)
    ch <- d1[, !nn, drop = FALSE]
    id <- names(ch)
    
    d <- data[, c(id, val, p)]
    setnames(d,c("ID","INTENSITY","PVALUE"))
  }else{
    ch <- d1[, !nn, drop = FALSE]
    id <- names(ch)
    
    if(nums>1)
      stop("Data is not of correct formate, please check the data files and try again")
    
    p <- names(nums)
    d <- data[,c(id,p)]
    setnames(d,c("ID","PVALUE"))
  }
  return(d)
}

#################################### End ##############################

## Function to read the processed data to store into a list
ready_data <- list()
proc_data <- function(path, leng){
  for (i in 1:leng)
  {  
    incProgress(0.1, detail = paste("File",i)) 
    process <- fread(path[[i]], data.table = FALSE, stringsAsFactors = FALSE)
    ready_data[[length(ready_data)+1]] = process
  }
  return(ready_data)
}

############################# End ######################################

 
shinyServer(function(input, output,session) { 
  
  #wd=getwd()
  #print(wd)
  #testdata=c("data.csv","data1.csv")
  #print(testdata)
  #testdata_full_path=path.expand(paste0(wd,"\\www\\",testdata))
  ## Reactive to accept input files from source directory
  infiles <- eventReactive(input$Submit, {        
    if (is.null(input$files))
    {
      return(NULL)
    }
    input$files
  })
  
  output$file <- renderDataTable({
    infiles()[,c(1, 2, 3)]
     
  })
  
  list_files <- eventReactive(input$Submit,{
    lst <-  input$files[['name']]  
  
    path1 <- input$files[['datapath']][1]
    
    lines <- readLines(path1,n=10)
   
    forIlm <- read.table(text=sub('.*(#.*)', '\\1', lines),
                        check.names=FALSE, stringsAsFactors=FALSE, 
                        comment.char='#',header=TRUE,fill=TRUE,row.names=NULL)
    #print(forIlm)
     
    na.omit(forIlm)
    names <- names(forIlm)
    cols <- ncol(forIlm)
    
    n <- grepl("bead", names,ignore.case=TRUE)
    n1 <- grepl("AVG_Signal", names,ignore.case=TRUE)
    n2 <- grep('^(detection|Pval)', names,ignore.case=TRUE)
    n3 <- all(length(names) >4 ,length(n2) >=2)
    
    cat("Is this a bead array:",any(n,n1,n3), "\n")
    
    an <- any(n,n1,n3)
     
   
  ext1 <- grepl("*.cel$", lst,ignore.case=TRUE)
  
  ext2 <- grepl("*.TXT$", lst)
   
  ################################################################
  #allext3 <- all(sapply(c('\\.CEL$', '\\.txt$'), function(x) any(grepl(x,lst, ignore.case=TRUE))))
    
  ## check for if the files belong to DE for codelink
  #ext4 <- c('TXT','txt')
  #allext4 <- all(ext4 %in% sub('.*\\.', '', lst))
  #print(allext4)
  
  ## checking if the files belong to DE of illumina bead array
  #ext5 <- grepl("*.txt$", lst)
  #chkbox <- input$checkbox
  #print(lst)
  
  ## To check if the checkbox is needed for differential expression
  #if(any(allext3,allext4)){
    #if(!(chkbox))
      #stop("You seem to be interested in differential expression, please make use of 
           #differential expression checkbox and try again")
  #}
  
  ## Test to check if check box is needed for illumina data. If the user loads differential
  ## expression data but forgot to check the differential expression box, the program will 
  ## try to identify the differential data here by dividing the list into names and targets
  
  #if(all(ext5)){
    #ilmpath = input$files[['datapath']]
    #nameLst <- ilmTest(ilmpath)
    #print(nameLst)
    
    ## creating seperate list for data files and phenotype files
    #ilmNames <- nameLst$names
    #ilmTargets <- nameLst$targets
    #if(length(ilmNames)>0 && length(ilmNames)<length(ilmpath)){
      #if(!chkbox)
        #stop("You seem to be interested in differential expression, please make use of 
             #differential expression checkbox and try again") 
    #}  
  #}
  
  #cat("Test for affymetrix DE:",allext3,"\n")
  #cat("Test for codelink DE:",allext4,"\n")
  #cat("Test for illumina DE:",all(ext5,chkbox),"\n")
    
  withProgress(message = 'Processing Data . . .',{
  if(all(ext1)) ## && !chkbox)
  {
    if(input$radio ==2)
      stop("Data files uploaded are not of codelink format, please check the data files and try again")
    
    if(input$radio ==3)
      stop("Data files uploaded are not of illumina format, please check the data files and try again")
    
    affy.data = ReadAffy(filenames = input$files[['datapath']])
      cat("Affymetrix Raw data selected\n")
      df <- processed_data(affy.data)
    }
  else if(all(ext2)) ## && !chkbox)
  {
    if(input$radio ==1)
      stop("Data uploaded is not of Affymetrix format, please check the data files and try again")
    
    if(input$radio ==3)
      stop("Data files uploaded is not of Illumina format, please check the data files and try again")
    
    cat("Codelink raw data selected\n")
    codelinkpath = input$files[['datapath']]
    code <- codelink_data(codelinkpath)
    code <- lapply(code,na.omit)
     
  }
  #else if(allext3 && chkbox)
    #{
    #if(input$radio ==2)
      #stop("Data uploaded is not of Codelink format, please check the data files and try again")
    
    #if(input$radio ==3)
      #stop("Data files uploaded is not of illumina format, please check the data files and try again")
    
    #cat("Checking for making contrasts of affymetrix data:")
    
    #fold <- input$num
    #if(is.na(fold))
      #stop("Foldchange is missing, please choose and try again")
    
    #pval <- input$num1
    #if(is.na(pval))
     # stop("p-value is missing, please choose and try again")
     
    #conditions <- input$text
    #if(is.na(conditions))
      #stop("Comparision information is missing, please provide and try again")
    
    #conditions <- strsplit(conditions,split=",")
    #print(conditions)
    
    #cat("Disease conditions choosen to compare:\n")
    
    #cat("Differential expression choosen\n")
    
    ## Renaming the files, as it gets changed on uploading the files
    #from <- input$files[['datapath']]
    #newnames <- file.path(dirname(from), basename(input$files[['name']]))
    #newone <- newnames[1]
    
    ##Removing the file name to set the path
    #newpath <- sub("/[^/]*$", "", newone) 
    #n <- file.rename(from, newnames)
    
    ## Calling affymetrix differential expression function
    #fc1 <- fold.change(lst,newpath,df1,fold,pval,newnames,conditions)
    #print(fc1)
    #}
  #else if(allext4 && chkbox)
  #{
    #if(input$radio ==1)
      #stop("Data uploaded is not of Affymetrix format, please check the data files and try again")
    
    #if(input$radio ==3)
      #stop("Data uploaded is not of illumina format, please check the data files and try again")
    
    #cat("Checking for making contrasts of codelink data:")
    #cntrst <- input$text1
    #print(cntrst)
    
    #if(is.na(cntrst))
      #stop("Comparison information is missing, please provide and try again")
    
    #fc <- input$num3
    #if(is.na(fc))
      #stop("Foldchange is missing, please provide and try again")
    
    #pval <- input$num4
    #if(is.na(pval))
      #stop("p-value is missing, please choose and try again")
    
    #cntrsts <- strsplit(cntrst,split=",")
    #print(cntrsts)
    
    #cat("checking for annotation DB:")
    #annotation <- input$text2
    #annotation <- input$select
    
    #if(is.na(annotation))
      #stop("Annotation data base is missing, please choose and try again")
    
    #print(annotation)
    
    #from <- input$files[['datapath']]
    #cat("From datapath:",from,"\n")
    #newnames <- file.path(dirname(from), basename(input$files[['name']]))
    #cat("New names:",newnames,"\n")
    #newone <- newnames[1]
    
    ##Removing the file name to set the path
    #newpath <- sub("/[^/]*$", "", newone) 
    #cat("New path:",newpath,"\n")
    #n <- file.rename(from, newnames)
    
    #codelink <- codelinkdiff(lst,newpath,cntrsts,pval,fc,annotation,newnames)
    #print(codelink)
  #}
  #else if(all(ext5) && chkbox)
  #{
    #cat("you are here in illumina")
    
    #if(length(ilmNames)>0 && length(ilmNames)<length(ilmpath)){
      
      #if(input$radio ==2)
        #stop("Data uploaded is not of codelink format, please check the data files and try again")
      
      #if(input$radio ==1)
       # stop("Data files uploaded is not of Affymetrix format, please check the data files and try again")
      
      ## To read the phenotype file to make sure the model matrix 
      ##construction is done correctly
      #from <- input$files[['datapath']]
      #newnames <- file.path(dirname(from), basename(input$files[['name']]))
      
      #x <- input$text3
      #cat("contrasts choosen:",x,"\n")
      
      #if(is.na(x))
        #stop("Make contrasts information is missing, please provide and try again")
      
      #cat("Spliting names of list:\n")
      #x <- strsplit(x,split=",")
      #print(x)
      
      #fc <- input$num5
      #if(is.na(fc))
       # stop("Foldchange information is missing, please choose and try again")
      
      #pval <- input$num6
      #if(is.na(pval))
        #stop("pvalue information is missing, please provide and try again")
      
      ## Both the phenotype files and raw data files are of .txt format, so we divide the
      ## list of files into phenotype and raw data files.
      #ilmpath = input$files[['datapath']]
      #nameLst <- ilmTest(ilmpath)
      #print(nameLst)
      
      ## creating seperate list for data files and phenotype files
      #ilmNames <- nameLst$names
      #ilmTargets <- nameLst$targets
      
      #cat("Name list done:\n")
      #print(ilmNames)
      
      #cat("Target list done:\n")
      #print(ilmTargets)
      
      #illumina <- illuminaDe(ilmNames,ilmTargets,x,fc,pval,newnames)
      
      #cat("Illumina df's for DE:\n")
      #print(head(illumina))
    #}
    #else if(length(ilmNames)==length(ilmpath)){
      #stop("No target files found, please check the files and try again")
    #}
    #else if(length(ilmNames)==0){
      #stop("No file with bead information found, please check your data files and try again")
    #}
  #}
  else if(any(n,n1,n3)) 
  {
    #if(chkbox)
      #stop("No file compartible for differential expression, please check the data file and try again")
    
    if(input$radio ==1)
      stop("Data uploaded is not of Affymetrix format, please check the data files and try again")
    
    if(input$radio ==2)
      stop("Data uploaded is not of Codelink format, please check the data files and try again")
      
    ## Illumina raw data processing
      cat("Illumina raw data processing:\n")
      illumina <- input$files[['datapath']]
      illuminaDat <- illuminaRaw(illumina)
      #print(illuminaDat)
  }
  #else if(all(ext1)|all(ext2) && chkbox){
    #stop("No target files found, please check target files and try again")
  #}
  #else if(all(ext5) && chkbox){
    #stop("No data files found, please check data files and try again")
  #}
    else{
    #cat("you are here\n")
    #if(chkbox)
      #stop("Are you sure of perfoming differential expression on the loaded data files,as the files 
           #are not compatible for differential expression, please check your data and try again")
    
      cat("Processed data  choosen\n")
      path <- input$files[['datapath']]
      leng <- length(path)
      
      ## Reading all the processed files and storing in a list
      procs <- proc_data(path,leng)
      procs <- lapply(procs, na.omit)
      #print(head(procs[[1]]))
      vec <- procs[[1]][1, ]
      
      aff <- sapply(vec,function(x) grepl("_at$",x))
    
      ilm <- sapply(vec,function(x) grepl("ILMN_",x))
      print(ilm)
    
      if('TRUE'%in%aff)
      {
        cat("This part belongs to affymetrix data\n")
        if(input$radio ==2)
          stop("Data uploaded is not of codelink format, please check the data files and try again")
        
        if(input$radio ==3)
          stop("Data files uploaded is not of Illumina format, please check the data files and try again")
        
        cat("Affymetrix data choosen: \n")
        procss <- lapply(procs, fun2)
        print(head(procss))
      }
      else if(!any(aff,ilm,an))##input$radio==2)
        {
        cat("Here only for codelink data\n")
        if(input$radio == 1)
          stop("Data uploaded is not of Affymetrix format, please check the data files and try again")
        
        if(input$radio == 3)
          stop("Data uploaded is not of Illumina format, please check the data files and try again")
          
        cat("Codelink data choosen\n")
        
        procss <- lapply(procs,fun1)
        #print(head(procss))
        
        if(length(procss)==0)
          stop("Data files seem not of correct formate, please check the data files and try again")
        
        if(length(procss)<length(procs))
          warning("Some file are skipped")
        
        procss <- Filter(function(x) !is.null(x), procss)
        procss <- lapply(procss, na.omit)
        #print(procs)
        } 
      else if('TRUE'%in%ilm)##input$radio==3)
      {
        cat("U should be here for ilumina\n")
        if(input$radio == 2)
          stop("Data uploaded is not of codelink format, please check the data files and try again")
        
        if(input$radio ==1)
          stop("Data uploaded is not of Affymetrix format, please check the data files and try again")
        
        cat("Illumina data choosen\n")
        
        procs <- lapply(procs,ilmnOrder)
        
        procs <- lapply(procs, na.omit)
      }
    }
  })
  })
  
  ## Displaying content of one file
  output$sourced <- renderDataTable({ 
    if (is.null(input$files))
    {
      return(NULL)
    }
      first_file <- list_files()[[1]]  
  })
  
  
  ## Reactive to select a file from the data-base (Annotation File)
  dataInput <- eventReactive(input$Submit, {
    
    if (grepl("[/\\\\]", input$dataset)) {
      stop("Invalid dataset")
    }
    annotation <- input$select
    
    #if(annotation=="h10kcod"|annotation=="h20kcod"|annotation=="hwgcod"){
      #cat("No annotation file to display,bioconductor annotation choosen\n")
    #}else{
      gpl <- read.csv(file.path("data", input$dataset),sep = '\t', comment.char = "#", check.names=FALSE,stringsAsFactors=FALSE)
      colnames(gpl) <- c("ID","NAME","SYMBOL")
      #print(head(gpl))
      gpl
    #}
  })
  
  ## To clean the annotation file
  inputdata <- eventReactive(input$Submit, {
    
    df1 <- dataInput()
    
    if (is.null(df1))
      return(NULL)
    
    df1[] <- lapply(df1, as.character)
    df1 <- df1[df1[,3]!='',]
    df1 <- df1[!grepl('[/]', df1$SYMBOL),]
    df1 <- df1[!grepl('NM_', df1$SYMBOL),]
    df1 <- df1[,c(1:3)]
    df1 <- na.omit(df1)
    
    cat("Testing GPL file done ...\n")
    df1
  }) 
  
  
  output$annotation <- renderDataTable({  ## To display the annotaion file
    withProgress(message = 'Loading Data . . .', {
      
      if (is.null(input$files))
        return(NULL)
      
      inputdata() ## dataInput() 
    })
  })
  
     
  ## Main Program ##
  data <- eventReactive(input$Submit,{
     
    if (is.null(input$files))
      return(NULL)
    
    if(is.null(input$dataset))
     return(NULL)
    
    df1 <-  inputdata()
    withProgress(  message = 'Annotating . . .',value= 0.01,{
      
      rout <- input$radio
      df <- list_files()
      files <- length(df)
       
      cat("You are in the main program... \n")
      cat("Sample frame generation done...\n")
      print(head(df[[1]]))
      
      #forFold <- grepl('AlogFC',colnames(df[[1]]))[2]
      #forFoldIlmn <- grepl('IlogFC',colnames(df[[1]]))[2]
      #print(forFold)
      
      #if(forFold)
        #{
          #foldchange <- fold.change1(df,df1)
          
          #if(nrow(foldchange)==0)
            #stop("Data files uploaded and annotation does not matching, please check the files and try again")
          
          #foldchange
        #}
      #else if(forFoldIlmn)
        #{
          #cat(" Processing Illumina differential data: ")
          #foldchange <- fold.change2(df,df1)
          
          #if(nrow(foldchange)==0)
            #stop("Data files uploaded and annotation does not matching, please check the files and try again")
          
          #foldchange
        #}
      
        if(rout==1 && any(length(df[[1]])==4,length(df[[1]])==3,length(df[[1]]==2)))
        {
          cat("Affymetrix  data selected \n")
          affyFin <- affymetrix(df,files,df1)
          affyFin <- data.frame(affyFin,stringsAsFactors=FALSE)
          print(affyFin)
        
          if(nrow(affyFin)==0)
            stop("Data files uploaded and annotation does not matching, please check the files and try again")
          
          #print(affyFin)
          affyFin 
        }
      else if(any(length(df[[1]])==3,length(df[[1]]==2)) && rout==2 )
        {
          cat("Codelink data selected \n")
          symbols <-  grepl("SYMBOL",colnames(df[[1]]))
          
          if(any(symbols)){
            cdlnkFin <- codlinkdiff1(df)
          }else{
            cdlnkFin <- codelink(df,files,df1)
          } 
          cdlnkFin <- data.frame(cdlnkFin,stringsAsFactors=FALSE)
          
          if(nrow(cdlnkFin)==0)
            stop("Data files uploaded and annotation does not matching, please check the files and try again")
          
          cdlnkFin 
        }
      else if(any(length(df[[1]])==3,length(df[[1]]==2)) && rout==3)
        {
        cat("Illumina data selected \n")
        ilmdata <- illumina(df,files,df1)
        #print(ilmdata)
        
        ilmdata <- data.frame(ilmdata,stringsAsFactors=FALSE)
        
        if(nrow(ilmdata)==0)
          stop("Data files uploaded and annotation does not matching, please check the files and try again") 
        
        ilmdata 
        }
    })
  })
  
  
   ## Summarizing data
  fulldata <- eventReactive(input$Submit,{ 
    withProgress(  message = 'Generating Final Summary . . .',value= 0.01,{
      
      cat ("Data Summary done ...\n")
      affy <- data()
    })
  })
    
  
  ## To  summarize the results and assign score, and then cummulative score.
  summary <- eventReactive(input$Submit,{
    final <- fulldata()
    
    if (is.null(final))
      return(NULL)
     
    if(ncol(final)<3)
      stop("Less than two data files are selected, please provide two or more file and try again")
    
    
    logical <- grepl('PVALUE',colnames(final))
    cols <- ncol(final)-1
    sm <- sum(logical, na.rm=TRUE)
    
    logical3 <- grepl('INTENSITY',colnames(final))
    logical4 <- grepl('CALL',colnames(final))
    logical1 <- grepl('SNR',colnames(final))[3]
    logical2 <- grepl('AlogFC|ClogFC|IlogFC',colnames(final))[2]
    chkbox <- input$checkbox
    radio1 <- input$radio1
    
    cat("Select the meta-analysis method\n")
    print(radio1)
  
    if(any(logical) && any(logical3) | sm == cols)
    {
      cat("Final data generation reached:\n")
      #print(head(final))
      
      final <- data.frame(final)
      finaldata <- final[grep('^(SYMBOL|PVALUE)', names(final))]
      finaldata[is.na(finaldata)] <- 0
      cat("i AM IN AFFYMETRIX PART\n")
      #print(head(finaldata))
      symbolpvalue <- finaldata
      
      ## converting the column factors
      indx <- sapply(finaldata[,-1], is.factor)[1]
      
      if(indx){
        finaldata[, 2:ncol(finaldata)] = lapply(finaldata[, -1], function(x) as.numeric(levels(x))[x])
      }
      intensity.cols <- grep("INTENSITY", names(finaldata))
      print(intensity.cols)
      pvalue.cols <- grep("PVALUE", names(finaldata))
      if(radio1 == 1)
      {
        cat("Fisher method choosen\n")
        ## Conducting Fishers test
        #finaldata$Fishers <- apply(finaldata[, names(finaldata)[pvalue.cols]], 1, 
                                   #function(x) Fisher.test(x[x!=0]))
        #finaldata$Fishers <- as.numeric(finaldata$Fishers)
        
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
        #print(finaldata)
        finaldata <- MetaDE.pvalue(finaldata, meta.method = 'Fisher', miss.tol = 0.9)
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        print(finaldata)
        ## Sorting the results in incresing order
        finaldata = finaldata[ order(finaldata$Fisher),]
        finaldata$FinalCall <- ifelse(finaldata$Fisher ==0, 0,
                                      ifelse(finaldata$Fisher < 0.05, 'P',
                                             ifelse(finaldata$Fisher >= 0.065,'A',0)))
      }
      else if(radio1 == 2)
      {
        ## Conducting Souffers test
        #finaldata$Stouffers <- apply(finaldata[, names(finaldata)[pvalue.cols]], 1, 
                                     #function(x) Stouffer.test(x[x!=0]))
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
        
        #finaldata$Stouffers <- as.numeric(finaldata$Stouffers)
        
        finaldata <- MetaDE.pvalue(finaldata, meta.method = 'Stouffer')
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        print(finaldata)
        
        ## Sorting the resutls in increasing order
        finaldata = finaldata[ order(finaldata$Stouffer),]
        finaldata$FinalCall <- ifelse(finaldata$Stouffer == 0, 0,
                                      ifelse(finaldata$Stouffer < 0.05, 'P',
                                             ifelse(finaldata$Stouffer >= 0.065, 'A', 0)))
      }
      else if(radio1==3)
      {
        ## Conducting AW.Fishers test
        finaldata$AW.Fishers <- apply(finaldata[, names(finaldata)[pvalue.cols]], 1, 
                                   function(x) sim.AW.Fisher(x[x != 0]))
        finaldata$AW.Fishers <- as.numeric(finaldata$AW.Fishers)
        
        ## Sorting the results in incresing order
        finaldata = finaldata[order(finaldata$AW.Fishers), ]
        finaldata$FinalCall <- ifelse(finaldata$AW.Fishers == 0, 0,
                                      ifelse(finaldata$AW.Fishers < 0.05, 'P',
                                             ifelse(finaldata$AW.Fishers >= 0.065,'A',0)))
      }
      else if(radio1 == 4)
        {
        ##minP test
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
        finaldata <- MetaDE.pvalue(finaldata, meta.method = 'minP')
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        finaldata <- finaldata[order(finaldata$minP), ]
        finaldata$FinalCall <- ifelse(finaldata$minP == 0, 0,
                                      ifelse(finaldata$minP < 0.05, 'P',
                                             ifelse(finaldata$minP >= 0.065,'A',0)))
        cat("Checking minP method\n")
        #print(finaldata)
      }
      
      else if(radio1 == 5)
      {
        ##minP test
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
        #print(finaldata)
        finaldata <- MetaDE.pvalue(finaldata, meta.method = 'maxP')
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        finaldata <- finaldata[order(finaldata$maxP), ]
        finaldata$FinalCall <- ifelse(finaldata$maxP == 0, 0,
                                      ifelse(finaldata$maxP < 0.05, 'P',
                                             ifelse(finaldata$maxP >= 0.065,'A',0)))
        cat("Checking minP method\n")
        #print(finaldata)
      }
      else if(radio1 == 6)
      {
        ##minP test
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
        print(finaldata)
        finaldata <- MetaDE.pvalue(finaldata, meta.method = 'SR')
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        finaldata <- finaldata[order(finaldata$SR), ]
        finaldata$FinalCall <- ifelse(finaldata$SR == 0, 0,
                                      ifelse(finaldata$SR < 0.05, 'P',
                                             ifelse(finaldata$SR >= 0.065,'A',0)))
        cat("Checking minP method\n")
        #print(finaldata)
      }
      else if(radio1 == 7)
      {
        ##minP test
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
        #print(finaldata)
        finaldata <- MetaDE.pvalue(finaldata, meta.method = 'PR')
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        finaldata <- finaldata[order(finaldata$PR), ]
        finaldata$FinalCall <- ifelse(finaldata$PR == 0, 0,
                                      ifelse(finaldata$PR < 0.05, 'P',
                                             ifelse(finaldata$PR >= 0.065,'A',0)))
        cat("Checking minP method\n")
        #finaldata
      }
      else if(radio1 == 8)
      {
        ##minP test
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
        #print(finaldata)
        finaldata <- MetaDE.pvalue(finaldata, meta.method = 'minP.OC')
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        finaldata <- finaldata[order(finaldata$minP.OC), ]
        finaldata$FinalCall <- ifelse(finaldata$minP.OC == 0, 0,
                                      ifelse(finaldata$minP.OC < 0.05, 'P',
                                             ifelse(finaldata$minP.OC >= 0.065,'A',0)))
        cat("Checking minP method\n")
        #print(finaldata)
      }
      else if(radio1 == 9)
      {
        ##minP test
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
        #print(finaldata)
        finaldata <- MetaDE.pvalue(finaldata, meta.method = 'maxP.OC')
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        finaldata <- finaldata[order(finaldata$maxP.OC), ]
        finaldata$FinalCall <- ifelse(finaldata$maxP.OC == 0, 0,
                                      ifelse(finaldata$maxP.OC < 0.05, 'P',
                                             ifelse(finaldata$maxP.OC >= 0.065,'A',0)))
        cat("Checking minP method\n")
        #print(finaldata)
      }
      else if(radio1 == 10)
      {
        cat("Fiser one sided correction\n")
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
        #print(finaldata)
        finaldata <- MetaDE.pvalue(finaldata, meta.method = 'Fisher.OC')
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        finaldata <- finaldata[order(finaldata$Fisher.OC), ]
        finaldata$FinalCall <- ifelse(finaldata$Fisher.OC == 0, 0,
                                      ifelse(finaldata$Fisher.OC < 0.05, 'P',
                                             ifelse(finaldata$Fisher.OC >= 0.065,'A',0)))
        cat("Checking minP method\n")
        #print(finaldata)
      }
      else if(radio1 == 11)
      {
        ##minP test
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
        #print(finaldata)
        finaldata <- MetaDE.pvalue(finaldata, meta.method = 'Stouffer.OC')
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        finaldata <- finaldata[order(finaldata$Stouffer.OC), ]
        finaldata$FinalCall <- ifelse(finaldata$Stouffer.OC == 0, 0,
                                      ifelse(finaldata$Stouffer.OC < 0.05, 'P',
                                             ifelse(finaldata$Stouffer.OC >= 0.065,'A',0)))
        cat("Checking minP method\n")
        #print(finaldata)
      }
        
      #finaldata$CummuScore <- rowSums(ifelse(finaldata[,-1]==0, 0, 
                                             #ifelse(finaldata[, -1]<= 0.05, 2,
                                                    #ifelse(finaldata[,-1]>= 0.065,-2,0))))
      #finaldata$FinalCall <- ifelse(finaldata$Fishers >0, 'P',
                                    #ifelse(finaldata$CummuScore < 0,'A',0))
      
      cat("Affymetrix final data score assignment done ...\n")
      finaldata
    }
    else if(!any(logical) && any(logical4))
    {
      cat("Final data has CALL and INTENSITY present\n")
      #print(head(final,20))
      
      final <- data.frame(final)
      finaldata <- final[grep('^(SYMBOL|CALL)', names(final))]
       
      finaldata[] = lapply(finaldata, as.character)
      finaldata[is.na(finaldata)] <- 0
      
      finaldata$CummuScore <- rowSums(ifelse(finaldata[,-1]=='P', 2, 
                                             ifelse(finaldata[, -1]=='A', -2,0)))
      
      finaldata$CummuScore <- as.numeric(finaldata$CummuScore)
      #print(head(finaldata,100))
      
      finaldata = finaldata[ order(finaldata$CummuScore,decreasing=TRUE),]
      
      finaldata$FinalCall <- ifelse(finaldata$CummuScore >0, 'P',
                                    ifelse(finaldata$CummuScore < 0,'A',0))
      
      cat("Affymetrix final data score assignment done ...\n")
      
      finaldata
    }
    else if(logical1)
    {
      finaldata <- final[grep('^(SYMBOL|SNR)', names(final))]
      finaldata[is.na(finaldata)] <- 0
      
      indx <- sapply(finaldata[,-1], is.factor)[1]
      
      if(indx)
        finaldata[, 2:ncol(finaldata)] = lapply(finaldata[, -1], function(x) as.numeric(levels(x))[x])
  
      finaldata$CummuScore <- rowSums(ifelse(finaldata[,-1]==0, 0, 
                                             ifelse(finaldata[,-1]>= 1,2,-2)))
      
      finaldata$CummuScore <- as.numeric(finaldata$CummuScore)
      
      finaldata = finaldata[ order(finaldata$CummuScore,decreasing=TRUE),] 
      
      finaldata$FinalCall <- ifelse(finaldata$CummuScore >0, 'P',
                                    ifelse(finaldata$CummuScore < 0,'A',0))
      
      cat("Codelink final data score assignment done ... \n")
      
      finaldata
    }
    else if(logical2 && chkbox)
    {
      cat("Performing Differential Expression\n")
      chk <- input$checkbox
       
      fold <- data.frame(final)
      finaldata <- fold[, c(1,2,seq(4,ncol(fold),2))]
      finaldata[is.na(finaldata)] <- 0
      
      finaldata$CummuScore <- rowSums(ifelse(finaldata[,-1]==0, 0, 
                                             ifelse(finaldata[, -1]>0, 2, -2)))
      
      finaldata$CummuScore <- as.numeric(finaldata$CummuScore)
      
      finaldata = finaldata[ order(finaldata$CummuScore,decreasing=TRUE),] 
      
      finaldata$FinalCall <- ifelse(finaldata$CummuScore >0, 'P',
                                    ifelse(finaldata$CummuScore < 0,'A',0))
      
      finaldata
    }
    else if(input$radio == 3)
    {
      intensity.cols <- grep("INTENSITY", names(finaldata))
      pvalue.cols <- grep("PVALUE", names(finaldata))
      
      # Extracting pvalue columns only
      finaldata <- final[, c(1,seq(3,ncol(final),2))]
      #print(head(finaldata))
      #finaldata[is.na(finaldata)] <- 0
    
      ## converting the column factors
      indx <- sapply(finaldata[,-1], is.factor)[1]
      if(indx){
      finaldata[, 2:ncol(finaldata)] = lapply(finaldata[, -1], function(x) as.numeric(levels(x))[x])
      }
      if(radio1 == 1)
      {
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
         
        finaldata <- MetaDE.pvalue(finaldata[finaldata != 0], meta.method = 'Fisher')
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        finaldata <- finaldata[order(finaldata$Fisher), ]
        finaldata$FinalCall <- ifelse(finaldata$Fisher == 0, 0,
                                      ifelse(finaldata$Fisher < 0.05, 'P',
                                             ifelse(finaldata$Fisher >= 0.065,'A',0)))
        cat("Checking minP method\n")
        #print(finaldata)
      }
      else if(radio1 == 2)
      {
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
         
        finaldata <- MetaDE.pvalue(finaldata, meta.method = 'Stouffer')
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        finaldata <- finaldata[order(finaldata$Stouffer), ]
        finaldata$FinalCall <- ifelse(finaldata$Stouffer == 0, 0,
                                      ifelse(finaldata$Stouffer < 0.05, 'P',
                                             ifelse(finaldata$Stouffer >= 0.065,'A',0)))
        cat("Checking minP method\n")
        #print(finaldata)
      }
      else if(radio1 == 3)
      {
        ## Conducting AW.Fishers test
        finaldata$AW.Fishers <- apply(finaldata[, names(finaldata)[pvalue.cols]], 1, 
                                      function(x) sim.AW.Fisher(x[x != 0]))
        finaldata$AW.Fishers <- as.numeric(finaldata$AW.Fishers)
        
        ## Sorting the results in incresing order
        finaldata = finaldata[order(finaldata$AW.Fishers), ]
        finaldata$FinalCall <- ifelse(finaldata$AW.Fishers == 0, 0,
                                      ifelse(finaldata$AW.Fishers < 0.05, 'P',
                                             ifelse(finaldata$AW.Fishers >= 0.065,'A',0)))
      }
      else if(radio1 == 4)
      {
        ##minP test
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
         
        finaldata <- MetaDE.pvalue(finaldata, meta.method = 'minP')
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        finaldata <- finaldata[order(finaldata$minP), ]
        finaldata$FinalCall <- ifelse(finaldata$minP == 0, 0,
                                      ifelse(finaldata$minP < 0.05, 'P',
                                             ifelse(finaldata$minP >= 0.065,'A',0)))
        cat("Checking minP method\n")
        #print(finaldata)
      }
      
      else if(radio1 == 5)
      {
        ##maxP test
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
         
        finaldata <- MetaDE.pvalue(finaldata, meta.method = 'maxP')
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        finaldata <- finaldata[order(finaldata$maxP), ]
        finaldata$FinalCall <- ifelse(finaldata$maxP == 0, 0,
                                      ifelse(finaldata$maxP < 0.05, 'P',
                                             ifelse(finaldata$maxP >= 0.065,'A',0)))
        cat("Checking minP method\n")
        #print(finaldata)
      }
      else if(radio1 == 6)
      {
        ##Sum of Ranks (SR) method
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
        #print(finaldata)
        finaldata <- MetaDE.pvalue(finaldata, meta.method = 'SR')
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        finaldata <- finaldata[order(finaldata$SR), ]
        #finaldata$FinalCall <- ifelse(finaldata$SR == 0, 0,
                                      #ifelse(finaldata$SR < 0.05, 'P',
                                             #ifelse(finaldata$SR >= 0.065,'A',0)))
        cat("Checking minP method\n")
        #print(finaldata)
      }
      else if(radio1 == 7)
      {
        ##Product of Ranks (PR) method
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
        #print(finaldata)
        finaldata <- MetaDE.pvalue(finaldata, meta.method = 'PR')
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        finaldata <- finaldata[order(finaldata$PR), ]
        finaldata$FinalCall <- ifelse(finaldata$PR == 0, 0,
                                      ifelse(finaldata$PR < 0.05, 'P',
                                             ifelse(finaldata$PR >= 0.065,'A',0)))
        cat("Checking minP method\n")
        #print(finaldata)
      }
      else if(radio1 == 8)
      {
        ##minP one sided correction
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
         
        finaldata <- MetaDE.pvalue(finaldata, meta.method = 'minP.OC')
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        finaldata <- finaldata[order(finaldata$minP.OC), ]
        finaldata$FinalCall <- ifelse(finaldata$minP.OC == 0, 0,
                                      ifelse(finaldata$minP.OC < 0.05, 'P',
                                             ifelse(finaldata$minP.OC >= 0.065,'A',0)))
        cat("Checking minP method\n")
        #print(finaldata)
      }
      else if(radio1 == 9)
      {
        ##maxP one sided correction
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
         
        finaldata <- MetaDE.pvalue(finaldata, meta.method = 'maxP.OC')
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        finaldata <- finaldata[order(finaldata$maxP.OC), ]
        finaldata$FinalCall <- ifelse(finaldata$maxP.OC == 0, 0,
                                      ifelse(finaldata$maxP.OC < 0.05, 'P',
                                             ifelse(finaldata$maxP.OC >= 0.065,'A',0)))
        cat("Checking minP method\n")
        #print(finaldata)
      }
      else if(radio1 == 10)
      {
        ##minP test
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
        #print(finaldata)
        finaldata <- MetaDE.pvalue(finaldata, meta.method = 'Fisher.OC')
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        finaldata <- finaldata[order(finaldata$Fisher.OC), ]
        finaldata$FinalCall <- ifelse(finaldata$Fisher.OC == 0, 0,
                                      ifelse(finaldata$Fisher.OC < 0.05, 'P',
                                             ifelse(finaldata$Fisher.OC >= 0.065,'A',0)))
        cat("Checking minP method\n")
        #print(finaldata)
      }
      else if(radio1 == 11)
      {
        ##Stouffer one sided correction 
        finaldata <- list(intensity = finaldata[, names(finaldata)[intensity.cols]],
                          pval = finaldata[, names(finaldata)[pvalue.cols]])
         
        finaldata <- MetaDE.pvalue(finaldata, meta.method = 'Stouffer.OC')
        FDR <- finaldata$meta.analysis$FDR
        finaldata <- finaldata$meta.analysis$pval
        finaldata <- cbind(symbolpvalue, finaldata)
        finaldata$FDR <- FDR
        finaldata <- finaldata[order(finaldata$Stouffer.OC), ]
        finaldata$FinalCall <- ifelse(finaldata$Stouffer.OC == 0, 0,
                                      ifelse(finaldata$Stouffer.OC < 0.05, 'P',
                                             ifelse(finaldata$Stouffer.OC >= 0.065,'A',0)))
        cat("Checking minP method\n")
        #print(finaldata)
      }
    }
     
    
    ## Assigning values and summing up.
    #finaldata$CummuScore <- rowSums(ifelse(finaldata[,-1]==0, 0, 
                                           #ifelse(finaldata[, -1]<= 0.05, 2, ifelse(finaldata[,-1]>= 0.065,-2,0))))
    
    ## Sorting Cummulative score in decreasing order
    #finaldata$CummuScore <- as.numeric(finaldata$CummuScore)
    
    #finaldata = finaldata[order(finaldata$CummuScore,decreasing=TRUE),] 
    
    #finaldata$FinalCall <- ifelse(finaldata$CummuScore >0, 'P',
                                  #ifelse(finaldata$CummuScore < 0,'A',0))
    #}
     finaldata
  })
    
  ## To display in table for final output
  output$final <- renderDataTable({
    
    if (is.null(data()))
      return(NULL)
  
    summary()  
  })
  
  ## To display in table with full information including p-value and signal intensity
  output$full <- renderDataTable({ 
    
    if (is.null(fulldata()))
    {
      return(NULL)
    }
    
    alldata <- fulldata()
    #print(head(alldata,20))
    
    logical <- grepl('INTENSITY',colnames(alldata))
    logical1 <- grepl("CALL",colnames(alldata))
    logical2 <- grepl("PVALUE",colnames(alldata))
    
    if(any(logical1) && !any(logical2))
      {
       
      finaldata <- alldata[grep('^(SYMBOL|INTENSITY|CALL)',names(alldata))]
      cat("print head of final data")
      #print(head(finaldata))
    }
    else if(any(logical))
    {
      cat("U R almost done...\n")
      finaldata <- alldata[grep('^(SYMBOL|INTENSITY|SIGNALINTENSITY|SNR|PVALUE)', names(alldata))]
      #print(head(finaldata))
    }else{
      alldata 
    }
  })
  
  ## To display in table for fold change
  output$fold <- renderDataTable({
    
    if (is.null(summary1()))
      return(NULL)
    
    summary1()  
  })
  
  
  ## To download summary data
  output$downloadData1 <- downloadHandler(
    filename <- function() { paste('Summary', '.csv', sep='')},
    content <- function(file) {
    write.csv(summary(), file)
    })
  
  
  ## Download Supplementary file
  output$downloadData2 <- downloadHandler(
    filename <- function() { paste('Supplementary', '.csv', sep='')},
    content <- function(file){
    write.csv(fulldata(), file)
    })
})

##################### End ###################################