###############################################################################################################################
###### This file calculates diversities based on data downloaded from PaleoDb & matched with ##################################
###### manually and automatically binned stratigraphic data, see Rasmussen et al. (201xx)    ##################################
###############################################################################################################################
###### before running specify work directory  and select time resolution ######################################################

library(RMark)
library(gdata)
library(iNEXT)
library(entropart)
library(stats4) # used for drop.levels()
library(vegan)
setwd("xx")
source("SQS3.r")
source("Rasmussen_et_al_functions.r")
load("Rasmussen_et_al.Rdata")


###################################
# variables loaded from Rassdata.Rdata:
# devi_tot = PBDB occurrences Lochkovian only, downloaded  14.09.2017
# ordo_tot = PBDB occurrences occurrences Cambrian - Lochkovian, , downloaded  14.09.2017
# ts.r = stage slices
# ts.lr = stage slice length in ma
# xoc = Rassmussen et al. binnings
# ordoRn = automatically binned units from RNames (Ordovician only)


###################################
################### clean downloads

ordox_tot <- subset(ordox_tot, ordox_tot$identified_rank == "genus" | ordox_tot$identified_rank == "species" )
ordox_tot <- subset(ordox_tot, ordox_tot$identified_no>0)
ordx <- ordox_tot[, c("accepted_name", "collection_no", "formation")]

devi_tot <- subset(devi_tot, devi_tot$identified_rank == "genus" | devi_tot$identified_rank == "species" )
devi_tot <- subset(devi_tot, devi_tot$identified_no>0)
devi <- devi_tot[, c("accepted_name", "collection_no", "formation")]


###################################
# merge binned units with PBDB data

#####merge all manually binned units
mo <- merge(xoc, ts.r, by.x= "oldest", by.y ="V1", all.y = TRUE) # merge time slices with binnings
my <- merge(mo, ts.r, by.x= "youngest", by.y ="V1", all.y = TRUE)
my[,7] <- my[,6]-my[,5]
mf <- my[,c("formation", "youngest", "oldest", "V7")]
ordxm <- merge (ordx, mf, by="formation") # merge binned formations with occurrences (this ranges only until Pridoli)
devixm <- cbind(as.data.frame(as.character(devi$formation)), as.data.frame(as.character(devi$accepted_name)), 
                as.data.frame(as.numeric(as.character(devi$collection_no))), rep("Lk", nrow(devi)), rep("Lk", nrow(devi)), 
                rep(1, nrow(devi)))
colnames(devixm) <- colnames(ordxm)
ordxm <- rbind(ordxm, devixm) # now the Devonian occurrences are added
ordxm <- subset(ordxm, ordxm$V7<=1) # subset of accuracy

###### Prepare automatically binned units
ordRn_tot <- subset(ordoRn, ordoRn$identified_rank == "genus" | ordoRn$identified_rank == "species" )
ordRn_tot <- subset(ordRn_tot, ordRn_tot$identified_no>0)
ordRn_tot <- subset(ordRn_tot, ordRn_tot$ts_count<3)
ordRnx <- ordRn_tot[, c("accepted_name", "collection_no", "oldest", "youngest")]
ordRnx_reduct <- subset(ordRnx, !ordRnx$oldest %in% c("Camb", "Sil") | !ordRnx$youngest %in% c("Camb", "Sil"))
ordRnx_reduct <- drop.levels(ordRnx_reduct)
ordxm_reduct <- subset(ordxm, !ordxm$oldest %in% ts.r[25:42, 1] | !ordxm$youngest %in% ts.r[25:42, 1])
ordxm_reduct <- ordxm_reduct[, c("accepted_name", "collection_no", "oldest", "youngest")]
ordRnx_combined <- rbind(ordRnx_reduct, ordxm_reduct)

###################################
#### further preparation of data

ord_all <- ordxm[,c(2,4,5)]
for (i in 1:nrow(ord_all)) {
  if(grepl(" ", ord_all$accepted_name[i])){ ord_all$accepted_name[i] <- strsplit(as.character(ord_all$accepted_name[i])," ")[[1]][1]}
}
ord_all <- na.omit(ord_all)
ord_all <- drop.levels(ord_all)

ord_all_r <- ordRnx_combined[,c(1,3,4)]
for (i in 1:nrow(ord_all_r)) {
  if(grepl(" ", ord_all_r$accepted_name[i])){ ord_all_r$accepted_name[i] <- strsplit(as.character(ord_all_r$accepted_name[i])," ")[[1]][1]}
}
ord_all_r <- na.omit(ord_all_r)
ord_all_r <- drop.levels(ord_all_r)

###################################
###################################
###################################
############ DIVERSITY CALCULATIONS

################### some raw counts
ox_all <- ab.mat(ord_all, ts.r)
ox_all_r <- ab.mat(ord_all_r, ts.r)
ox.r <-  ox_all
ox.rr <-  ox_all_r
ox.r[ox.r > 0] <- 1
ox.rr[ox.rr > 0] <- 1
colSums(ox.rr)-colSums(ox.r) # difference between number of genus occurrences/bin of automated and manual approach
((colSums(ox.rr)-colSums(ox.r))*100)/(colSums(ox.rr)+colSums(ox.r)) # same as above in %

###################################
############## SQS & Hill diversity

SQSoll <- sqs.div(ox_all, ts.r) # manually
SQSollr <- sqs.div(ox_all_r, ts.r) # automated
Dord <- hill.div(ox_all, ts.r) # manually
Dordr <- hill.div(ox_all_r, ts.r) # automated

###################################
##################### CMR diversity

###### prepare inp files
in_ord <- in.mat(ord_all, ts.r)
in_ord_r <- in.mat(ord_all_r,ts.r)

###### models
Phitime=list(formula=~time)
ptime=list(formula=~time)
Penttime=list(formula=~time)
Gammatime=list(formula=~time)

###### POPAN complete time slices, manually binned
all.processed <- process.data(in_ord, model="POPAN", time.intervals = ts.lr)
all.ddl <- make.design.data(all.processed)
all.time.ch <- mark(all.processed,all.ddl, model.parameters=list(Phi=Phitime, p=ptime, pent=Penttime)) #chat = 22.14
results.popan.all.ch <- all.time.ch$results$derived$N
sampl.T <- all.time.ch$results$real[54:107,]

# check one time results
plot(2:54, results.popan.all.ch[2:54,1], type="o", pch=21, bg= 142, ylab="CMR manually Filtered")
arrows(2:54, all.time.ch$results$derived$N[2:54,2],2:54, all.time.ch$results$derived$N[2:54,3], length=0.05, angle=90, code=3)

# collect (by repeating) a number of calculations and calculate means
estims <- matrix(,25,54)
uci <- matrix(,25,54)
lci <- matrix(,25,54)
sestims <- matrix(,25,54)
suci <- matrix(,25,54)
slci <- matrix(,25,54)

all.processed <- process.data(in_ord, model="POPAN", time.intervals = ts.lr)
all.ddl <- make.design.data(all.processed)
for (k in 1:25) # 25 repetitions
{
  all.time <- mark(all.processed,all.ddl, model.parameters=list(Phi=Phitime, p=ptime, pent=Penttime))
  results.popan.all <- all.time$results$derived$N
  sampl.T <- all.time$results$real[54:107,]
  if (sum(results.popan.all[2:54,1])>0)
  {
    estims[k,] <- results.popan.all[,1]
    lci[k,] <- results.popan.all[,2]
    uci[k,] <- results.popan.all[,3]
    sestims[k,] <- sampl.T[,1]
    slci[k,] <- sampl.T[,3]
    suci[k,] <- sampl.T[,4]
  }
}

mest <- colMeans(estims, na.rm = TRUE) # manually binned
mlci <- colMeans(lci, na.rm = TRUE)
muci <- colMeans(uci, na.rm = TRUE)
smest <- colMeans(sestims, na.rm = TRUE) # manually binned
smlci <- colMeans(slci, na.rm = TRUE)
smuci <- colMeans(suci, na.rm = TRUE)

plot(1:54, mest, type="o", pch=21, bg= 142, , ylab="CMR means manually Filtered")
arrows(1:54, mlci,1:54, muci, length=0.05, angle=90, code=3)
plot(1:54, smest, type="o", pch=21, bg= 142, , ylab="sampling")
arrows(1:54, smlci,1:54, smuci, length=0.05, angle=90, code=3)

###### POPAN complete time slices, automatically binned 
r.processed <- process.data(in_ord_r, model="POPAN", time.intervals = ts.lr)
r.ddl <- make.design.data(r.processed)
r.time.ch <- mark(r.processed, r.ddl, model.parameters=list(Phi=Phitime, p=ptime, pent=Penttime)) #chat = 22.14
results.popan.r.ch <- r.time.ch$results$derived$N

# check one time results
plot(2:54, results.popan.r.ch[2:54,1], type="o", pch=21, bg= 142, ylab="CMR RNames Filtered")
arrows(2:54, r.time.ch$results$derived$N[2:54,2],2:54, r.time.ch$results$derived$N[2:54,3], length=0.05, angle=90, code=3)

# collect (by repeating) a number of calculations and calculate means
estimsR <- matrix(,25,53)
uciR <- matrix(,25,53)
lciR <- matrix(,25,53)

r.processed <- process.data(in_ord_r, model="POPAN", time.intervals = ts.lr)
r.ddl <- make.design.data(r.processed)

for (k in 1:25) # 25 repetitions
{
  r.time <- mark(r.processed,r.ddl, model.parameters=list(Phi=Phitime, p=ptime, pent=Penttime))
  results.popan.r <- r.time$results$derived$N
  if (sum(results.popan.r[2:54,1])>0)
  {
    estimsR[k,] <- results.popan.r[2:54,1]
    lciR[k,] <- r.time$results$derived$N[2:54,2]
    uciR[k,] <- r.time$results$derived$N[2:54,3]
  }
}

mestr <- colMeans(estimsR, na.rm = TRUE)
mlcir <- colMeans(lciR, na.rm = TRUE)
mucir <- colMeans(uciR, na.rm = TRUE)

plot(2:54, mestr, type="o", pch=21, bg= 136, ylab="CMR dots:manual; squares:automated ")
lines(2:54,mest)
points(2:54, mest, type="p", pch=22, bg= 136)
arrows(2:54, mlci,2:54, muci, length=0.05, angle=90, code=3)
arrows(2:54, mlcir,2:54, mucir, length=0.05, angle=90, code=3)


###################################
###################################
###################################
###################################
###### TESTS

###### Number of observations
N_obs <- colSums(ox_all)
N_obsr <- colSums(ox_all_r)

###### Number of formations
ofm <- unique(ordxm[,c(1,4,5)])
bin.name <- array()
r1x <- ofm[which(as.character(ofm$oldest)==as.character(ofm$youngest)),]
r1 <- r1x[c("formation", "oldest")]
r2x <- ofm[which(!as.character(ofm$oldest)==as.character(ofm$youngest)),]
r2o <- r2x[c("formation", "oldest")]
r2y <- r2x[c("formation", "youngest")]
colnames(r2y) <- c("formation", "oldest")
rc <- rbind(r1, r2y, r2o)
rc <- drop.levels(rc)
rx <- table(rc$formation, rc$oldest)
orxh <- matrix(,nrow(rx), nrow(ts.r))

for (k in 1:(nrow(ts.r)))
{
  cns <- as.character(ts.r[k,1])
  bin.name[k] <- cns
  if (cns %in% colnames(rx) == TRUE) 
  {
    orxh[,k] <- rx[,cns]
  }
}
colnames(rxh) <- bin.name
orxh[is.na(rxh)] <- 0
nforms <- colSums(orxh)

###### calculations for SOM

lr.forms <- lm(SQSoll$Sobs ~ nforms)
summary(lr.forms)
lr.formm <- lm(N_obs ~ SQSoll$Sobs)
summary(lr.formm)
lr.forms <- lm(Dord$Sobs ~ SQSollr$S_SQS)
summary(lr.forms)

pdf(file="plots_1.pdf",width=11,height=5)
par(mfrow = c(1, 2)) 
plot(nforms, SQSoll$Sobs, xlab="Nfms", ylab="Sobs")
plot(N_obs, SQSoll$Sobs, xlab="Noccs", ylab="Sobs") # occurrences versus Sobs
dev.off()

pdf(file="plots_2.pdf",width=11,height=4)
par(mfrow = c(1, 3)) 
plot(SQSollr$Sobs, SQSollr$S_SQS, xlab="Sobs", ylab="D_SQS")
plot(Dordr$Sobs, Dordr$S_Hill, xlab="Sobs", ylab="D_Hill")
plot(Dordr$Sobs[1:54], mest, xlab="Sobs", ylab="D_CMR")
dev.off()

pdf(file="plots_3.pdf",width=11,height=4)
plot(1:54, smest, type="o", pch=21, bg= 136, ylab="p(sample")
arrows(1:54, smlci, 1:54, smuci, length=0.05, angle=90, code=3)
dev.off()

pdf(file="plots_4.pdf",width=11,height=4)
plot(2:54, mestr[1:53], type="o", pch=21, bg= 136, ylab="D_CMR") # dots:automated & squares:manually
lines(1:54,mest[1:54])
points(1:54, mest[1:54], type="p", pch=22, bg= 136)
arrows(2:54, mlci[2:54], 2:54, muci[2:54], length=0.05, angle=90, code=3)
arrows(2:54, mlcir[1:53],2:54, mucir[1:53], length=0.05, angle=90, code=3)
dev.off()

pdf(file="plots_5.pdf",width=11,height=4)
plot(1:54, mest, type="o", pch=21, bg= 136, ylab="n_gen") #dots:CMR; Squares: SQS; diamonds: Hill
lines(1:54,SQSoll$S_SQS[1:54]*2)
points(1:54, SQSoll$S_SQS[1:54]*2, type="p", pch=22, bg= 136)
lines(1:54,Dord$S_Hil[1:54]*2)
points(1:54, Dord$S_Hil[1:54]*2, type="p", pch=23, bg= 136)
arrows(1:54, mlci, 1:54, muci, length=0.05, angle=90, code=3)
arrows(1:54, (Dord$`+ci`[1:54])*2, 1:54, (Dord$`-ci`[1:54])*2, length=0.05, angle=90, code=3)
dev.off()

pdf(file="plots_finalCMR.pdf",width=11,height=4)
plot(1:54, mest[1:54], type="o", pch=21, bg= 136, ylab="D_CMR") # dots:automated & squares:manually
arrows(1:54, mlci[1:54], 1:54, muci[1:54], length=0.05, angle=90, code=3)
dev.off()



