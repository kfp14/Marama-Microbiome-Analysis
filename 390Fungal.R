#save early and save often
#need to allocate more memory to R: right click on R, properties, target, and add --max-mem-size=6GB or however much you can give
#let R use all the cores while it works; also add --p-n-thread 0
#try to just have R running if possible
#make sure to save and then clear out things you don't need from the global environment when large processing steps are running

#set working directory
setwd("~/Fungal Analysis")

#install necessary packages first; see Bacterial workflow for installation
#load packages
library("BiocManager")

library("knitr")
library("BiocStyle")
library("dada2")
library("ggplot2")

set.seed(100)

#load the data
miseq_path <- file.path("c:/", "Users", "Larry", "Documents", "Fungal Analysis", "Fungal Samples") #replace with your file path to folder containing sequence reads
fns <- sort(list.files(miseq_path,full.names=TRUE)) #pull out sequences from folder
fnFs <- fns[grepl("R1",fns)] #separates forward reads
fnRs <- fns[grepl("R2",fns)] #separates reverse reads

#trim and filter sequences using dada package
ii <- sample(length(fnFs), 3) #picks 3 samples to test for quality

for(i in ii) { print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd")) } #quality plot for forward reads; keep bps where Phred>20; mean is green, median is solid orange, quartiles are dotted orange
for(i in ii) { print(plotQualityProfile(fnRs[i]) + ggtitle("Rev")) } #quality plot for reverse reads

filt_path <- file.path("c:/", "Users", "Larry", "Documents", "Fungal Analysis","filtered") #replace with your own file path to working directory
if(!file_test("-d", filt_path)) dir.create(filt_path) #created folder for filtered and trimmed sequences
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[[i]], fnRs[[i]]),
                    c(filtFs[[i]], filtRs[[i]]),
                    trimLeft=10, truncLen=c(245,245),
                    maxN=0, maxEE=2, truncQ=2,
                    compress=TRUE)
} #filters and trims the sequences, populating new folder filtered

#infer Ribosomal Sequence Variants with DADA2

derepFs <- derepFastq(filtFs) #demux and dereplicate forwards
derepRs <- derepFastq(filtRs) #demux and dereplicate reverse
sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sam.names
names(derepRs) <- sam.names

ddF <- dada(derepFs, err=NULL, selfConsist=TRUE, MAX_CONSIST=5) #if you have enough memory, can increase loops up to 20
ddR <- dada(derepRs, err=NULL, selfConsist=TRUE, MAX_CONSIST=5) #uses parameter learning to estimate sequencing error rate

plotErrors(ddF, nominalQ=TRUE) #plot error parameters to check fit to data
plotErrors(ddF, nominalQ=TRUE)

dadaFs <- dada(derepFs, err=ddF[[1]]$err_out,pool=TRUE,multithread=TRUE) #pooled interference, takes longer but increases power to detect rare variants
dadaRs <- dada(derepRs, err=ddF[[1]]$err_out,pool=TRUE,multithread=TRUE) #these steps bin the RSVs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs) #merge forward and reverse reads, pairs must overlap

#make the sequence table with dada2

seqtab.all <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
#save!

seqtab <- removeBimeraDenovo(seqtab.all, multithread=TRUE) #removes chimeric sequences from PCR erros
#save!

#assign taxonomy and make taxa table with dada2

ref_fasta <- "c://Users/Larry/Documents/Fungal Analysis/fungalits_warcup_trainingdata_V2/fungalits_warcup_trainingdata_V2/Warcup_v2.fasta" #replace with file directory to your training set; this is Warcup V2 ITS 
taxtab <- assignTaxonomy(seqtab, refFasta = ref_fasta,multithread=TRUE)
colnames(taxtab) <- c("Kingdom","Phlyum","Class","Order","Family","Genus")
#save!

#multiple sequence alignment with DECIPHER package

library(DECIPHER)
seqs <- getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
write.csv(taxtab,"ITS1taxatable.csv") #exports as excel file with taxonomic assignments

#construct phylogenetic tree with phangorn package

library(phangorn)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align) 
fitGTR <- update(fit, k=4, inv=0.2)
#save before this step!
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "NNI", control = pml.control(trace = 1)) #takes a while to run; try rearrangement = "stochastic" if you have lots of processing power, better results
detach("package:phangorn", unload=TRUE)
#save the tree!

#validate sample metadata and make sure it matches to sequence table

mimarks_path <- "c://Users/Larry/Documents/Fungal Analysis//fungalmetadata.csv" #replace with your file path
samdf <- read.csv(mimarks_path, header=TRUE)
samdf <- samdf[!duplicated(samdf$SampleID),] # Remove dupicate entries for reverse reads
all(rownames(seqtab) %in% samdf$SampleID) #should return TRUE
rownames(samdf) <- samdf$SampleID
keep.cols <- c("location","Description", "SampleID") #the metadata file can be as simple as this three columns
samdf <- samdf[rownames(seqtab), keep.cols]

#merge the sample-by-sequence feature table, metadata, taxa table, and tree into a phyloseq object to analyze data downsteam
#once you get the phyloseq object the data analysis is a breeze!

library(phyloseq)
fungi <- phyloseq(tax_table(taxtab), sample_data(samdf),
                  otu_table(seqtab, taxa_are_rows = FALSE),phy_tree(fitGTR$tree))

#normalize number of reads per sample using median sequencing depth
total = median(sample_sums(fungi))
standf = function(x, t=total) round(t * (x / sum(x)))
fungi.norm = transform_sample_counts(fungi, standf)

#remove unassigned taxa then normalize again
fungi.0 <- subset_taxa(fungi, !is.na(Phlyum) & !Phlyum %in% c("", "uncharacterized"))
total.0 = median(sample_sums(fungi.0))
standf.0 = function(x, t=total.0) round(t * (x / sum(x)))
fungi.norm.0 = transform_sample_counts(fungi.0, standf.0)

#note: because this is an understudied system, the NAs are likely not just artifacts
#I suggest repeating the graphic analysis twice, once with and once without NAs
#also note that here I had spelled "Phylum" wrong in assigning headers to taxtab

#basic bar graph by taxa division
plot_bar(fungi.norm, fill = "Phlyum", title="Phylum Abundance with NAs") + 
  geom_bar(aes(color=Phlyum, fill=Phlyum), stat="identity", position="stack") #phylum abundance with NA
plot_bar(fungi.norm.0, fill = "Phlyum", title="Phylum Abundance without NAs") + 
  geom_bar(aes(color=Phlyum, fill=Phlyum), stat="identity", position="stack") #phylum abundance without NA

plot_bar(fungi.norm, fill = "Class", title="Class Abundance with NAs") + 
  geom_bar(aes(color=Phlyum, fill=Phlyum), stat="identity", position="stack") #starts to get a little uninformative past class
plot_bar(fungi.norm.0, fill = "Class", title="Class Abundance without NAs") + 
  geom_bar(aes(color=Phlyum, fill=Phlyum), stat="identity", position="stack") #I'm not sure if you should make a new object by removing unassigned taxa at class level and then renormalizing to plot this, but it looks fine

#ascomycota was by far the most abundant; we'll keep just that and plot a bar graph according to genus

fungi.asco <- subset_taxa(fungi.norm, Phlyum %in% c("Ascomycota")) #pull out the ascomycota phylum
total.asco = median(sample_sums(fungi.asco)) #we need to renormalize
standf.asco = function(x, t=total.asco) round(t * (x / sum(x)))
fungi.norm.asco = transform_sample_counts(fungi.asco, standf.asco)
plot_bar(fungi.norm.asco, fill = "Order", title="Order Abundance in Phylum Ascomycota") +
  geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack") #plot by order; too many families to be informative
#obviously the NAs didn't pass the ascomycota filter so we won't have an NA plot

#alpha diversity metrics; dada2 removed singletons so warning message doesn't apply
plot_richness(fungi.norm, measures=c("Chao1", "Shannon"), title="Alpha Diversity Metrics with NAs")
plot_richness(fungi.norm.0, measures=c("Chao1", "Shannon"), title="Alpha Diversity Metrics without NAs")

#phylogenetic trees

plot_tree(fungi, method = "treeonly",
                   ladderize = "left",
                   title = "Fungal Tree Before Agglomeration with NAs") #tree for original unfiltered data with NAs

plot_tree(fungi.0, method = "treeonly",
          ladderize = "left",
          title = "Fungal Tree Before Agglomeration without NAs") #tree for original unfiltered data without NAs

#tree based on taxonomic agglomeration at genus rank; doesn't include NAs
length(get_taxa_unique(fungi.0, taxonomic.rank = "Genus")) #tells you how many genera present after filtering
fungi.0.glom = tax_glom(fungi.0, "Genus", NArm = TRUE)
plot_tree(fungi.0.glom, method = "treeonly",
                   ladderize = "left", title = "Fungal Tree Agglomerated By Genus") #has no NAs because they don't have taxonomic assignment

#tree based on height without taxonomy to include unassigned RSVs; kind of like OTU clustering methods
h1 = 2 #I had to set the tree distance pretty high to get anything meaningful out of this
fungi.height = tip_glom(fungi, h = h1)
plot_tree(fungi.height, method = "treeonly",
                   ladderize = "left", title = "Fungal Tree Agglomerated By Height")