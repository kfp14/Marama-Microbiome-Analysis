#set working directory
setwd("~/Microbiome Analysis") #replace with whatever folder you're working in

#install updated R from inside RGui: if necessary, run once
install.packages("installr")
library(installr)
updateR()

#Help --> Check for updates from inside RStudio and download new version: if necessary, do once

#install BioConductor version 3.11 ONCE

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library("BiocManager") #load package to install others, use to install other packages

#install packages: run ONCE
BiocManager::install("knitr")
BiocManager::install("BiocStyle")
BiocManager::install("dada2")
BiocManager::install("phyloseq")
BiocManager::install("DECIPHER")
BiocManager::install("phangorn")
BiocManager::install("ggplot2")
BiocManager::install("gridExtra")

#load packages: run each session
library("knitr")
library("BiocStyle")
library("dada2")
library("ggplot2")
library("gridExtra")

set.seed(100)

#load the data
miseq_path <- file.path("c:/", "Users", "kfpar", "OneDrive", "Documents", "Microbiome Analysis", "Bacterial Sampled") #set the appropriate file path to get to folder containing fastq sequences
fns <- sort(list.files(miseq_path, full.names = TRUE)) #pulls out the fastq sequences from folder above
fnFs <- fns[grepl("R1", fns)] #separates forward reads
fnRs <- fns[grepl("R2", fns)] #separates reverse reads

#trim and filter the sequences
ii <- sample(length(fnFs), 3) #picks 3 seqs to sample for quality
for(i in ii) { print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd")) } #quality plot for forward reads; Phred>20 is good; mean is  green, median is solid orange, quartiles are dotted orange
for(i in ii) { print(plotQualityProfile(fnRs[i]) + ggtitle("Rev")) } #quality plot for reverse reads

filt_path <- file.path("c:/", "Users", "kfpar", "OneDrive", "Documents", "Microbiome Analysis", "filtered") #replace with your own file path
if(!file_test("-d", filt_path)) dir.create(filt_path) #creates a new folder called filtered for the filtered trimmed sequences
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))

for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[[i]], fnRs[[i]]),
                    c(filtFs[[i]], filtRs[[i]]),
                    trimLeft=10, truncLen=c(230, 240),
                    maxN=0, maxEE=2, truncQ=2,
                    compress=TRUE)
} #trims, filters, and populates the new folder filtered

for(i in ii) { print(plotQualityProfile(filtFs[i]) + ggtitle("Fwd")) } #check to see number of reads is reduced by filtering process; or check that file size went down

#NOTE ON MEMORY: Dada2 requires a lot of memory
#check how much is available to R (in MB) with command memory.limit()
#to allot more, close RStudio, right click on R, properties, shortcut, target, and add --max-mem-size=6GB
#also dedicate all available cores with --p-n-thread 0
#try not to be doing much else while this runs

#DADA2 to infer ribosomal sequence variants

derepFs <- derepFastq(filtFs) #demux and dereplicate forward reads
derepRs <- derepFastq(filtRs) #demux and dereplicate reverse reads
sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sam.names
names(derepRs) <- sam.names

#default max consist terminates after 10 loops without convergence; I cut it short by reducing to 5 and checking that the parameters still fit well
ddF <- dada(derepFs, err=NULL, selfConsist=TRUE, MAX_CONSIST=5) #will use parameter learning to estimate error rate and remove substitution and indel errors
ddR <- dada(derepRs, err=NULL, selfConsist=TRUE, MAX_CONSIST=5) #cutting it off at six rounds; will check error fit

plotErrors(ddF, nominalQ=TRUE) #check fit; fitted error rates black lines fit observations black points and decrease w increasing QS
plotErrors(ddR, nominalQ=TRUE) #even though I cut it off to save time the fit looks good

dadaFs <- dada(derepFs, err=ddF[[1]]$err_out, pool=FALSE) #pool=TRUE --> pooled interference, increases power to detect rare variants but requires more time and more processing power
dadaRs <- dada(derepRs, err=ddR[[1]]$err_out, pool=FALSE) #if you do run pooled interference, add multithread=TRUE to speed it up

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs) #merges forward and reverse and takes out nonoverlapping pairs

#make sequence table using dada package

seqtab.all <- makeSequenceTable(mergers[!grepl("mock",names(mergers))])
seqtab <- removeBimeraDenovo(seqtab.all,multithread=TRUE) #remove chimeras

ref_fasta <- "c://Users/kfpar/OneDrive/Documents/Microbiome Analysis/gg_13_8_train_set_97.fa.gz" #uses the greengenes training set; replace with your own file path
taxtab <- assignTaxonomy(seqtab, refFasta = ref_fasta)
colnames(taxtab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
View(taxtab) #allows you to see all assignments
write.csv(taxtab,"16Staxatable.csv") #makes an excel file with the taxonomic assignments
taxtab.final <- replace(taxtab, taxtab == "o__", NA) #remove cells that aren't completely empty; this will interfere with filtering specific taxa later
taxtab.final <- replace(taxtab.final, taxtab.final == "k__", NA)
taxtab.final <- replace(taxtab.final, taxtab.final == "p__", NA)
taxtab.final <- replace(taxtab.final, taxtab.final == "c__", NA)
taxtab.final <- replace(taxtab.final, taxtab.final == "f__", NA)
taxtab.final <- replace(taxtab.final, taxtab.final == "g__", NA)
taxtab.final <- replace(taxtab.final, taxtab.final == "s__", NA)

#multiple sequence alignment with DECIPHER package

library("DECIPHER")
seqs <- getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA) #multiple sequence alignment, MSA

#construct phylogenetic tree using phangorn
library("phangorn")
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "NNI", control = pml.control(trace = 1)) #if this doesn't work, here are other options:

#option 2 for phylo tree: modified fitGTR call - you could also try rearrangement = "NNI" or rearrangement = "non"
fitGTR <- optim.pml(fitGTR, rearrangement = "stochastic",ratchet.par=list(iter=5L, maxit=5L,prop=1/3), control = pml.control(trace = 1))

#option 3 for phylo tree: RAxML, better for modern large datasets. requires download of executable from Github
BiocManager::install("ips")
library("ips") #install and load ips package
out3 <- file.path("c://Users/Larry/Documents/Bacterial Analysis/seqs.fasta") #make a file path for a fasta file
writeXStringSet(alignment,out3,format="fasta") #write the alignment into a fasta file
align.raxml <- read.dna("seqs.fasta",format="fasta",as.matrix=TRUE) #read it back in to a DNA bin class object
exec <-"c://Users/Larry/Documents/Bacterial Analysis/raxmlHPC-PTHREADS-SSE3.exe" #file path to executable
alignment.rax.gtr <- raxml(align.raxml,
                           m="GTRGAMMAIX", # model
                           f="a", # best tree and bootstrap
                           p=1234, # random number seed
                           x=2345, # random seed for rapid bootstrapping
                           N=10, # number of bootstrap replicates
                           file="alignment", # name of output files
                           exec=exec,
                           threads=10
) #run this

#validate sample metadata and make sure it matches to sequence table

mimarks_path <- "c://Users/Larry/Documents/Bacterial Analysis/bactmetadata.csv" #replace with your file path
samdf <- read.csv(mimarks_path, header=TRUE)
samdf <- samdf[!duplicated(samdf$SampleID),] # Remove dupicate entries for reverse reads
all(rownames(seqtab) %in% samdf$SampleID) #should return TRUE
rownames(samdf) <- samdf$SampleID
keep.cols <- c("location","Description", "SampleID") #the metadata file can be as simple as this three columns
samdf <- samdf[rownames(seqtab), keep.cols]

#combine feature table, metadata, taxonomic assignments, and tree into a phyloseq object for downstream analysis

library("phyloseq")
bact <- phyloseq(tax_table(taxtab.final), sample_data(samdf),
                  otu_table(seqtab, taxa_are_rows = FALSE),phy_tree(fit$tree))

#normalize number of reads per sample using median sequencing depth
total = median(sample_sums(bact))
standf = function(x, t=total) round(t * (x / sum(x)))
bact.norm = transform_sample_counts(bact, standf)

#basic bar graph by taxa division

plot_bar(bact.norm, fill = "Kingdom", title="Kingdom Abundance with NAs") + 
  geom_bar(aes(color=Kingdom, fill=Kingdom), stat="identity", position="stack") #kingdom abundance with NA

#filter out just archaea

bact.arch <- subset_taxa(bact.norm, Kingdom %in% c("k__Archaea"))
total.arch = median(sample_sums(bact.arch)) #we need to renormalize
standf.arch = function(x, t=total.arch) round(t * (x / sum(x)))
bact.norm.arch = transform_sample_counts(bact.arch, standf.arch)

plot_bar(bact.norm.arch, fill = "Phylum", title="Archaeal Phylum Abundance") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") #phylum abundance for archaea

plot_bar(bact.norm.arch, fill = "Class", title="Archaeal Class Abundance") + 
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")

plot_bar(bact.norm.arch, fill = "Order", title="Archaeal Order Abundance") + 
  geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack")

#filter out just bacteria

justbact <- subset_taxa(bact.norm, Kingdom %in% c("k__Bacteria"))
total.bact = median(sample_sums(justbact)) #we need to renormalize
standf.bact = function(x, t=total.bact) round(t * (x / sum(x)))
justbact.norm = transform_sample_counts(justbact, standf.bact)

plot_bar(justbact.norm, fill = "Phylum", title="Bacterial Phylum Abundance") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

#biggest phylum was actinobacteria, let's subset to plot a descriptive breakdown of the phylum

bact.actin <- subset_taxa(justbact, Phylum %in% c("p__Actinobacteria")) #pull out the phylum
total.actin = median(sample_sums(bact.actin)) #we need to renormalize
standf.actin = function(x, t=total.actin) round(t * (x / sum(x)))
bact.norm.actin = transform_sample_counts(bact.actin, standf.actin)

plot_bar(bact.norm.actin, fill = "Order", title="Order Abundance in Phylum Actinobacteria") +
  geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack")

plot_bar(bact.norm.actin, fill = "Family", title="Family Abundance in Phylum Actinobacteria") +
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") #not super descriptive, pretty even split

#second biggest phylumwas proteobacteria, subset that

bact.proteo <- subset_taxa(justbact, Phylum %in% c("p__Proteobacteria")) #pull out the phylum
total.proteo = median(sample_sums(bact.proteo)) #we need to renormalize
standf.proteo = function(x, t=total.proteo) round(t * (x / sum(x)))
bact.norm.proteo = transform_sample_counts(bact.proteo, standf.proteo)

plot_bar(bact.norm.proteo, fill = "Order", title="Order Abundance in Phylum Proteobacteria") +
  geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack")

#alpha diversity metrics; dada2 removed singletons so warning message doesn't apply
plot_richness(bact.norm, measures=c("Chao1", "Shannon"), title="Alpha Diversity Metrics with NAs")

#phylogenetic trees

plot_tree(bact, method = "treeonly",
          ladderize = "left",
          title = "Bacterial and Archaeal Tree Before Agglomeration with NAs") #tree for original unfiltered data with NAs

#tree based on taxonomic agglomeration at genus rank; doesn't include NAs
length(get_taxa_unique(bact, taxonomic.rank = "Genus")) #tells you how many genera present after filtering
bact.glom = tax_glom(bact, "Genus", NArm = TRUE)
plot_tree(bact.glom, method = "treeonly",
          ladderize = "left", title = "Bacterial and Archaeal Tree Agglomerated By Genus") #has no NAs because they don't have taxonomic assignment