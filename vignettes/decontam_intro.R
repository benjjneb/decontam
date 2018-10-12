## ----loadPS----------------------------------------------------------------
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")
ps <- readRDS(system.file("extdata", "MUClite.rds", package="decontam"))
ps

## ----see-meta-table--------------------------------------------------------
head(sample_data(ps))

## ----see-depths------------------------------------------------------------
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

## ----frequency-------------------------------------------------------------
contamdf.freq <- isContaminant(ps, method="frequency", conc="quant_reading")
head(contamdf.freq)

## ----table-----------------------------------------------------------------
table(contamdf.freq$contaminant)
head(which(contamdf.freq$contaminant))

## ----plot-abundance, warning=FALSE-----------------------------------------
plot_frequency(ps, taxa_names(ps)[c(1,3)], conc="quant_reading") + 
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

## ----see-contams, warning=FALSE--------------------------------------------
set.seed(100)
plot_frequency(ps, taxa_names(ps)[sample(which(contamdf.freq$contaminant),3)], conc="quant_reading") +
    xlab("DNA Concentration (PicoGreen fluorescent intensity)")

## ----remove----------------------------------------------------------------
ps
ps.noncontam <- prune_taxa(!contamdf.freq$contaminant, ps)
ps.noncontam

## ----prevalence------------------------------------------------------------
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

## ----prevalence-05---------------------------------------------------------
contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

## ----see-prev-05-----------------------------------------------------------
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                      contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

