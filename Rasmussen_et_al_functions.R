###### produces abundance matrix from entry matrix with accepted_names | youngest | oldest
ab.mat<-function(range,tsx) 
{
  bin.name <- array()
  r1x <- range[which(as.character(range$oldest)==as.character(range$youngest)),]
  r1 <- r1x[c("accepted_name", "oldest")]
  r2x <- range[which(!as.character(range$oldest)==as.character(range$youngest)),]
  r2o <- r2x[c("accepted_name", "oldest")]
  r2y <- r2x[c("accepted_name", "youngest")]
  colnames(r2y) <- c("accepted_name", "oldest")
  rc <- rbind(r1, r2y, r2o)
  rc <- drop.levels(rc)
  rx <- table(rc$accepted_name, rc$oldest)
  rxh <- matrix(,nrow(rx), nrow(tsx))
  
  for (k in 1:(nrow(tsx)))
  {
    cns <- as.character(tsx[k,1])
    bin.name[k] <- cns
    if (cns %in% colnames(rx) == TRUE) 
    {
      rxh[,k] <- rx[,cns]
    }
  }
  colnames(rxh) <- bin.name
  rxh[is.na(rxh)] <- 0
  return(rxh)
}


#### produces in mat for marked at Rasmussen slices level

in.mat<-function(range,tsx) 
{
  bin.name <- array()
  r1x <- range[which(as.character(range$oldest)==as.character(range$youngest)),]
  r1 <- r1x[c("accepted_name", "oldest")]
  r2x <- range[which(!as.character(range$oldest)==as.character(range$youngest)),]
  r2o <- r2x[c("accepted_name", "oldest")]
  r2y <- r2x[c("accepted_name", "youngest")]
  colnames(r2y) <- c("accepted_name", "oldest")
  rc <- rbind(r1, r2y, r2o)
  rc <- drop.levels(rc)
  rx <- table(rc$accepted_name, rc$oldest)
  rxh <- matrix(,nrow(rx), nrow(tsx))
  
  for (k in 1:(nrow(tsx)))
  {
    cns <- as.character(tsx[k,1])
    bin.name[k] <- cns
    if (cns %in% colnames(rx) == TRUE) 
    {
      rxh[,k] <- rx[,cns]
    }
  }
  colnames(rxh) <- bin.name
  rxh[is.na(rxh)] <- 0
  
  rxh[rxh > 0] <- 1
  rxh[is.na(rxh)] <- 0
  rxh <- subset(rxh, rowSums(rxh)>0 )
  cmrarr <- as.data.frame(matrix(,nrow(rxh), 1))
  for(kk in 1:nrow(rxh)) 
  {
    cmrarr[kk,] <- paste(rxh[kk,1], rxh[kk,2], rxh[kk,3], rxh[kk,4], rxh[kk,5], rxh[kk,6], rxh[kk,7], rxh[kk,8], rxh[kk,9], 
                         rxh[kk,10], rxh[kk,11], rxh[kk,12], rxh[kk,13], rxh[kk,14], rxh[kk,15], rxh[kk,16], rxh[kk,17], rxh[kk,18],
                         rxh[kk,19], rxh[kk,20], rxh[kk,21], rxh[kk,22], rxh[kk,23], rxh[kk,24], rxh[kk,25], rxh[kk,26], rxh[kk,27],
                         rxh[kk,28], rxh[kk,29], rxh[kk,30], rxh[kk,31], rxh[kk,32], rxh[kk,33], rxh[kk,34], rxh[kk,35], rxh[kk,36],
                         rxh[kk,37], rxh[kk,38], rxh[kk,39], rxh[kk,40], rxh[kk,41], rxh[kk,42], rxh[kk,43], rxh[kk,44], rxh[kk,45],
                         rxh[kk,46], rxh[kk,47], rxh[kk,48], rxh[kk,49], rxh[kk,50], rxh[kk,51], rxh[kk,52], rxh[kk,53], rxh[kk,54], sep="") # with Sil
    cmrarr[kk,] <- as.character(cmrarr[kk,])
  }
  #cmrarr <- as.data.frame(cmrarr)
  colnames(cmrarr) <- "ch"
  return(cmrarr)
}

hill.div<-function(rxh,tsx)
{
  Cm <- array()
  bin.name <- array()
  md2a <- matrix(,(nrow(tsx)),5)
  rf.s <- array()
  
  for (k in 1:nrow(tsx)) 
  {
    rf.k <- rxh[,k]
    
    if (sum(rf.k)>99)
    {
      Cm[k] <-Coverage(rxh[,k]) # using Package entropart
    }
  }
  Cm <- Cm[!is.na(Cm)]
  
  for (i in 1:(nrow(tsx))) {
    cns <- as.character(tsx[i,1])
    bin.name.a <- cns
    bin.name <- rbind(bin.name, bin.name.a)
    if (cns %in% colnames(rxh) == TRUE) {
      rf.s <- rxh[,cns]
      
      md2a[i,1] <- sum(rf.s)
      md2a[i,2] <- length(rf.s[rf.s>0])
      md2a[i,3] <- length(rf.s[rf.s>0])
      md2a[i,4] <- 0
      md2a[i,5] <- 0
      
      if (sum(rf.s)>99)
      {
        rarx <- estimateD(rf.s, datatype = "abundance", base = "coverage", level = min(Cm), conf = 0.95) # q=1 used Shannon Entropy
        md2a[i,3] <- rarx[1,5]
        md2a[i,4] <- rarx[1,6]
        md2a[i,5] <- rarx[1,7]
      }
    }
  }
  row.names(md2a) <- tsx[,1]
  colnames(md2a) <- c("collections", "Sobs", "S_Hill", "+ci", "-ci")
  md2a[is.na(md2a)] <- 0
  md2bins <- as.data.frame(md2a)
  return(md2bins)
}

###### produces SQS
sqs.div<-function(rxh,tsx)
{
  Cm <- array()
  bin.name <- array()
  md2a <- matrix(,(nrow(tsx)),3)
  rf.s <- array()
  
  for (k in 1:nrow(tsx)) 
  {
    rf.k <- rxh[,k]
    if (sum(rf.k)>99)
    {
      Cm[k] <-Coverage(rxh[,k]) # using Package entropart
    }
    Cm <- Cm[!is.na(Cm)]
    
  }
  
  for (i in 1:(nrow(tsx))) 
  {
    cns <- as.character(tsx[i,1])
    bin.name.a <- cns
    bin.name <- rbind(bin.name, bin.name.a)
    if (cns %in% colnames(rxh) == TRUE) 
    {
      rf.s <- rxh[,cns]
      md2a[i,1] <- sum(rf.s)
      md2a[i,2] <- length(rf.s[rf.s>0])
      md2a[i,3] <- length(rf.s[rf.s>0])
      
      if (sum(rf.s)>99)
      {
        rarx <- sqs(rf.s[rf.s>0], min(Cm)) # q=1 used Shannon Entropy
        md2a[i,3] <- rarx[3]
      }
    }
  }
  row.names(md2a) <- tsx[,1]
  colnames(md2a) <- c("collections", "Sobs", "S_SQS")
  md2a[is.na(md2a)] <- 0
  md2bins <- as.data.frame(md2a)
  return(md2bins)
}