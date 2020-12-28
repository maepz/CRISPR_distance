#init

library(ggplot2)
library(reshape2)
library(gtools)
library(dplyr)


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.10")
# devtools::install_github("benjjneb/dada2")

library(dada2); packageVersion("dada2")
library("devtools")

setwd(utils::getSrcDirectory(function(x) {x})[1])



path<- "./Amplicons_fastq/V4_and_unmapped"
list.files(path)

#sort by names
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz"))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz"))
#replace names
sample.names <- sapply(strsplit(fnFs,"_"), `[`, 2)
sample.names 

#### IMPORT DATA ####
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

#sort by lib name

sample.names <- mixedsort(sample.names) 
fnFs <- mixedsort(fnFs)
fnRs <- mixedsort(fnRs)
sample.names
#### CHECK DATA QUAL ####

plotQualityProfile(fnFs[1]) # 1st Forward sample
plotQualityProfile(fnRs[1]) # 1st Reverse sample

#### FILTER AND TRIM ####

filt_path <- file.path(path, "/filtered_pairedend") 
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncQ=6,
                     truncLen = c(200,200),
                     trimLeft=c(19,20),
                     maxEE=c(2,2),multithread=T)


keep <- out[,"reads.out"] > 20 # Or other cutoff
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))[keep]
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))[keep]

#### DEREPLICATE ####

derepFs <- derepFastq(filtFs)
names(derepFs) <- sample.names[keep]

derepRs <- derepFastq(filtRs)
names(derepRs) <- sample.names[keep]

# apprentissage du taux d'erreur
errF <- learnErrors(filtFs)
errR <- learnErrors(filtRs)


dadaFs <- dada(derepFs, 
               err = errF, 
               multithread=TRUE,
               pool=TRUE)

dadaRs <- dada(derepRs, 
               err=errR,
               multithread=TRUE,
               pool=TRUE)

dadaFs[[1]]
dadaRs[[1]]

#save(dadaRs, file="data/dadaRs.rdata")
#save(dadaFs, file="data/dadaFs.rdata")

#### MERGE ####

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, 
                      minOverlap = 12, 
                      maxMismatch = 0)
#check results
head(mergers[[1]])
max(mergers[[1]]$nmatch) # Largest overlap 
min(mergers[[1]]$nmatch) # Smallest overlap    


#### TABLEAU DES ASV ###
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
seqtab[,1]
seqtab[,2] 
seqtab[,3]
#show on graph
hist(nchar(getSequences(seqtab)),xlab="Size", ylab="Frequency", main = "ASVs length", xlim=c(250,450), ylim=c(0,250)) 

#### CHIMERA REMOVAL ####
?removeBimeraDenovo

seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                    method = "pooled", 
                                    multithread = TRUE,
                                    verbose = TRUE) 
save(seqtab.nochim, file=paste0(path,"/seqtab.nochim.rdata"))
write.csv(seqtab.nochim, paste0(path,"/seqtab.nochim.csv"), col.names=NA)
length(seqtab.nochim)
round((sum(seqtab.nochim)/sum(seqtab)*100),2) # Percentage of the total sequence reads
hist(nchar(getSequences(seqtab.nochim)),xlab="Size", ylab="Frequency", main = "Non-chimeric ASVs length", xlim=c(250,450), ylim=c(0,250)) # Lenght of the non-chimeric sequences

seqtab.nochim.bin <- ifelse(seqtab.nochim>0,1,0) 

#### TABLEAU RESUME ####
out2 = out
rownames(out2) = sapply(strsplit(row.names(out2),"_"), `[`, 2)
out2 = subset(out2, rownames(out2) %in% rownames(seqtab))
nrow(out2)

out = out2
length(out[,1])
length(out[,2])

length(as.numeric(sapply(mergers, getN)))


getN <- function(x) sum(getUniques(x))
track <- data.frame(Input=as.numeric(out[,1]), # input
                    Filtered=as.numeric(out[,2]), # filtered
                    "Filt//In"=as.numeric(round(((out[,2]/out[,1])*100),2)),# % (Filtered / Input)
                    Merge = as.numeric(sapply(mergers, getN)), # Merged 
                    "Mer//In"=as.numeric(round(((sapply(mergers, getN)/out[,1])*100),2)),# % (Merged / Input)
                    Nonchim = as.numeric(rowSums(seqtab.nochim)),# Non-chimeric                       
                    "Nonchim//In"=as.numeric(round(((rowSums(seqtab.nochim)/out[,1])*100),2)),# % (Non-chimeric / Input)
                    ASV = as.numeric(rowSums(seqtab.nochim.bin))) # Number of ASVs per sample 
rownames(track) <- rownames(seqtab)#sample.names # Row names

head(track)

# graph that
gtrack<- track[,c(1,2,4,6)]
gtrack$ID <- rownames(gtrack)

lgtrack <- melt(gtrack, id.vars="ID")
bar_track <- ggplot(lgtrack ,aes(x=ID, y=as.numeric(value), fill=variable)) +
  geom_bar(stat="identity", position = "identity") + 
  theme_classic() + # Theme
  theme(axis.ticks.length=unit(0.3,"cm")) + # Ticks size
  theme(axis.text.x = element_text(angle=45) , legend.title = element_blank())+ # Changes the x labels orientation & delete legend title
  scale_x_discrete(name ="Sample ID", limits=rownames(track))+ # Changes x-axis title & sorts the x label names
  scale_y_continuous(name="Abundance", breaks=seq(from = 0, to = 1000, by = 100))+ #Changes y-axis title & sets the y breaks.
  ggtitle("Track")# Main title
bar_track 


#### FIGURE HAPLOTYPES ####
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz"))
meta=read.table('../Amplicons_metadata.txt',sep='\t',header=T)

str(as.data.frame(seqtab.nochim))
rownames(seqtab.nochim)
colnames(seqtab.nochim)=c(1:(ncol(seqtab.nochim)))
df = melt(seqtab.nochim, )

df %>% head(.)
colnames(df)[1]='Sample_id'
colnames(df)[2]='Haplo'

df$Sample_id = as.factor(substr(df$Sample_id,1,7))
df %>% head(.)
df = merge(df,meta,by='Sample_id',all.x=T)
df = transform(df, Genotype=match(Haplo, unique(Haplo)))
df$Genotype = as.character(df$Haplo)
# df$Genotype = as.factor(factor(df$Genotype, levels = c(24,11,15,32)))
# df$Genotype[is.na(df$Genotype)] <- 'Others'
library(plyr)
df$Genotype = mapvalues(df$Genotype, from = c("1", "2",'3','4'), to = c("V4_1", "V4_2",'V4_3','V4_4'))

df$Region = factor(df$Region, levels = c('MV','CB','MEF'))
df$Flow = factor(df$Flow, levels = c('H','L','B'))
df$Habitat = factor(paste(df$Region,df$Flow,sep='-'), levels=c('MV-H','MV-L','CB-H','CB-B','MEF-H','MEF-L','MEF-B'))
dcast(df,value~Sample)
df=merge(df,aggregate(value~Sample,data=df,FUN = 'sum'), by='Sample',all.x=T)

colnames(df)[4]='count'
colnames(df)[13]='total'

### Saved as svg: 540X286
## haplotypes with more than 10 reads in samples 
ggplot(subset(df,count >10))+
  theme_bw()+
  geom_bar(aes(x=Sample,y=count,fill=Genotype),stat='identity',position='stack')+
  scale_fill_brewer(palette = "Dark2")+
  xlab('Trophosome sample')+
  ylab('Reads count')+
  theme(axis.text.x = element_text(angle = 60, vjust=1, hjust = 1.1))+
  facet_grid(.~Habitat, scales='free',space = "free")

## High-confidence haplotypes only
ggplot(subset(df, !(is.na(Genotype))))+
  theme_bw()+
  geom_bar(aes(x=Sample,y=count,fill=Genotype),stat='identity',position='stack')+
  scale_fill_brewer(palette = "Dark2")+
  xlab('Trophosome sample')+
  ylab('Reads count')+
  theme(axis.text.x = element_text(angle = 60, vjust=1, hjust = 1.1))+
  facet_grid(.~Habitat, scales='free',space = "free")

ggplot(subset(df,!(is.na(Genotype))))+
  theme_bw()+
  geom_bar(aes(x=Sample,y=count,fill=Genotype),stat='identity',position='stack')+
  scale_fill_brewer(palette = "Dark2",na.value = "grey50")+
  xlab('Trophosome sample')+
  ylab('Reads count')+
  theme(axis.text.x = element_text(angle = 60, vjust=1, hjust = 1.1))+
  facet_grid(.~Habitat, scales='free',space = "free")


#### duplicates fisher tests ####


dupl=subset(df,  Duplicate=='y')
dupl=dcast(data=dupl,Sample~Genotype,value.var='value.x',sum)
rownames(dupl)=dupl$Sample
# same
dupl1=(dupl[c("R13L-5-T_bis","R13L-5-T"),c((2:(ncol(dupl)-4)))])
fisher.test(dupl1)
diversity(t(dupl2))  %>% mean(.)

# Same
dupl2=dupl[c("R15H-3-T_old","R15H-3-T"),c((ncol(dupl)-3):ncol(dupl))]
dupl2=(dupl[c("R15H-3-T_old","R15H-3-T"),c((2:(ncol(dupl)-4)))])

fisher.test(dupl2)
diversity(t(dupl2))  %>% mean(.)

# Same
dupl3=dupl[c("R16L-3-T_bis","R16L-3-T"),c((ncol(dupl)-3):ncol(dupl))]
dupl3=(dupl[c("R16L-3-T_bis","R16L-3-T"),c((2:(ncol(dupl)-4)))])

fisher.test(dupl3)
head(meta)
diversity(t(dupl3))  %>% mean(.)

