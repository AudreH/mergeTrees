# Préparation des données

load("GSE87650.RData")

dat1.1 = exprs(gse_GSE87650[[1]])
dim(dat1.1)
dat2.1 = exprs(gse_GSE87650[[2]])
dim(dat2.1)

dat2.1 = apply(dat2.1, 2, as.numeric)

phen1.1 = pData(gse_GSE87650[[1]])
phen2.1 = pData(gse_GSE87650[[2]])

phen1 = pData(gse_GSE87650[[1]])
phen2 = pData(gse_GSE87650[[2]])

phen2 =  phen2[which(phen2$data_row_count!=0),]
dat2 = dat2.1[,which(colnames(dat2.1)%in%rownames(phen2))]

dat1 = dat1.1[,order(colnames(dat1.1))]
dat2 = dat2[,order(colnames(dat2))]
phen1 = phen1[order(rownames(phen1)),]
phen2 = phen2[order(rownames(phen2)),]

phen2$`cell type:ch1`[phen2$`cell type:ch1`=="wh blood"] = "WB"
phen2$`cell type:ch1`[phen2$`cell type:ch1`=="monocytes"] = "CD14"

phen1$ID = paste(phen1$`samplenumber:ch1`, phen1$`cell type:ch1`, sep = "_")
phen2$ID = paste(phen2$`patient_number:ch1`, phen2$`cell type:ch1`, sep = "_")

phen2 = phen2[-which(duplicated(phen2$ID)),]
dat2 = dat2[, which(colnames(dat2)%in%rownames(phen2))]

rownames(phen1) = colnames(dat1) = phen1$ID
rownames(phen2) = colnames(dat2) = phen2$ID

patient_id = intersect(phen1$ID , phen2$ID)
phen1 = phen1[which(phen1$ID%in%patient_id),]
phen2 = phen2[which(phen2$ID%in%patient_id),]

dat1 = dat1[, which(colnames(dat1)%in%rownames(phen1))]
dat2 = dat2[, which(colnames(dat2)%in%rownames(phen2))]
head(dat2)

dim(dat1)
dim(dat2)

dat1 = dat1[,order(colnames(dat1))]
dat2 = dat2[,order(colnames(dat2))]
phen1 = phen1[order(rownames(phen1)),]
phen2 = phen2[order(rownames(phen2)),]

save(dat1, dat2, phen1, phen2, file = "data_tablesformat.RData")