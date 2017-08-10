## ----loadPS----------------------------------------------------------------
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")
ps <- readRDS(system.file("extdata", "MUClite.rds", package="decontam"))
ps

## ----see-meta--------------------------------------------------------------
sample_variables(ps)

## ----see-depths------------------------------------------------------------
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

## ----frequency-------------------------------------------------------------
contam.freq <- isContaminant(ps, method="frequency", conc="quant_reading", )
head(contam.freq)

## ----table-----------------------------------------------------------------
table(contam.freq)
head(which(contam.freq))

## ----plot-abundance, warning=FALSE-----------------------------------------
plot_frequency(ps, taxa_names(ps)[c(1,3)], conc="quant_reading")

## ----see-contams, warning=FALSE--------------------------------------------
set.seed(100)
plot_frequency(ps, taxa_names(ps)[sample(which(contam.freq),3)], conc="quant_reading")

## ----remove----------------------------------------------------------------
ps
ps.noncontam <- prune_taxa(!contam.freq, ps)
ps.noncontam

## ----prevalence------------------------------------------------------------
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
contam.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contam.prev)
head(which(contam.prev))

## ----prevalence-05---------------------------------------------------------
contam.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contam.prev05)

## ----see-prev-05-----------------------------------------------------------
# Make phyloseq object of presence-absence in negative controls
ps.neg <- prune_samples(sample_data(ps)$Sample_or_Control == "Control Sample", ps)
ps.neg.presence <- transform_sample_counts(ps.neg, function(abund) 1*(abund>0))
# Make phyloseq object of presence-absence in true positive samples
ps.pos <- prune_samples(sample_data(ps)$Sample_or_Control == "True Sample", ps)
ps.pos.presence <- transform_sample_counts(ps.pos, function(abund) 1*(abund>0))
# Make data.frame of prevalence in positive and negative samples
df.pres <- data.frame(prevalence.pos=taxa_sums(ps.pos.presence), prevalence.neg=taxa_sums(ps.neg.presence),
                      contam.prev=contam.prev)
ggplot(data=df.pres, aes(x=prevalence.neg, y=prevalence.pos, color=contam.prev)) + geom_point()

