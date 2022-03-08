require(googledrive)
require(dplyr)
drive_download(file = "WORKSPACE_masterdata_SNPinfo.RData")
load("WORKSPACE_masterdata_SNPinfo.RData")

#individuals in each set
ind.only.LD <- setdiff(masterdata.LD$treeID, masterdata.HD$treeID)
ind.name <- c(masterdata.HD$treeID, ind.only.LD)

#genotype changes
cur.geno <- apply(expand.grid(c("A", "T", "C", "G"), c("A", "T", "C", "G")), 1, paste, collapse = "") 
names(cur.geno) <- apply(expand.grid(c("T", "A", "G", "C"), c("T", "A", "G", "C")), 1, paste, collapse = "") 

# filetring low density data for individuasl corresponding to high density 
LD.HD.ind <- masterdata.LD %>%
  filter(treeID%in%masterdata.HD$treeID)

# filetring high density data for markers corresponding to low density 
HD.LD.snp <- masterdata.HD[,colnames(masterdata.LD)] %>% 
  filter(treeID%in%masterdata.LD$treeID)

# variables LD.HD.ind and HD.LD.snp contain same individials and snps
identical(names(LD.HD.ind), names(HD.LD.snp))

# SNP name
snp.names <- colnames(LD.HD.ind)[-c(1:4)]

# combined output
out <- matrix(NA, length(ind.name), length(snp.names))
for(i in 1:length(snp.names)){
  cat(".")
  if(i%%100 == 0) cat("\n")
  d1 <- HD.LD.snp %>% select(c("treeID", snp.names[i]))
  d2 <- LD.HD.ind %>% select(c("treeID", snp.names[i]))
  D <- dplyr::left_join(d1, d2, by = "treeID")
  id <- D[which(D[,2] != D[,3]),]
  if(nrow(id) <= 10){
    out[,i] <- c(masterdata.HD[,snp.names[i]], filter(masterdata.LD, treeID%in%ind.only.LD)[,snp.names[i]])
  } else{
    out[,i] <- c(masterdata.HD[,snp.names[i]], cur.geno[filter(masterdata.LD, treeID%in%ind.only.LD)[,snp.names[i]]])
  }
}
res<-cbind(rbind(masterdata.HD[,1:4],
                 filter(masterdata.LD, treeID%in%ind.only.LD)[,1:4]),
           out)
colnames(res) <- c(colnames(masterdata.HD)[1:4], snp.names)
res <- as_tibble(res)
masterdata.FULL <- bind_rows(res, masterdata.HD[,setdiff(colnames(masterdata.HD)[-c(1:4)], snp.names)])
save(masterdata.FULL, file = "mastedata_FULL.rda")






