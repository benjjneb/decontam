# plot_abundance function, which plots the abundance of selected taxa in each sample vs. each sample's DNA concentration after PCR amplification

plot_abundance <- function(ps, taxa_to_plot, conc, taxa_are_rows=TRUE, norm=TRUE, log=TRUE, returndf=FALSE){
  
  ot <- as(otu_table(ps), "matrix")
  if(taxa_are_rows){
    ot <- t(ot)
  }
  
  taxa_mismatch <- !(taxa_names(ps) %in% colnames(ot))
  if(sum(taxa_mismatch) > 0){
    stop("Error: 'taxa_are_rows' argument may not be correct")
  }
  
  if(norm){
    ot <- sweep(ot, 1, rowSums(ot), "/")
  }
  
  ot <- as(ot[,colnames(ot) %in% taxa_to_plot], "matrix")
  if(dim(ot)[2] == 1){
    colnames(ot) <- taxa_to_plot
  }
  
  st <- as(sample_data(ps), "data.frame")
  snames <- sample_names(ps)
  
  plot <- merge(st, ot, by.x = "row.names", by.y = "row.names", sort=FALSE)
  plot_melt <- melt(plot, id.vars = colnames(plot)[1:(dim(plot)[2]-length(taxa_to_plot))])
  
  colnames(plot_melt)[dim(plot_melt)[2]-1] <- "taxa_to_plot"
  colnames(plot_melt)[dim(plot_melt)[2]] <- "taxon_abundance"
  
  taxon_levels <- taxa_to_plot
  plot_melt$taxa_to_plot <- factor(plot_melt$taxa_to_plot, levels = taxon_levels)
  
  I <- which(colnames(plot_melt) == conc)
  colnames(plot_melt)[I] <- "DNA_conc"
  
  if(returndf == FALSE){
    if(log==TRUE){
      p1 <- ggplot(plot_melt, aes(log(DNA_conc),log(taxon_abundance)))
      return(p1 + geom_point())
    } else if(log==FALSE){
      p1 <- ggplot(plot_melt, aes(DNA_conc,taxon_abundance))
      return(p1 + geom_point())
    }
  } else if(returndf == TRUE){
    return(plot_melt)
  }
}

##Examples
#plot_abundance(fungal,"denovo7","qPCR_copies",taxa_are_rows=TRUE)

#plot_abundance(MUC,"Seq1","quant_reading", taxa_are_rows=FALSE)
#plot_abundance(MUC,c("Seq152","Seq1"),"quant_reading",taxa_are_rows=FALSE)
#plot_abundance(MUC,c("Seq152","Seq1"),"quant_reading",taxa_are_rows=FALSE) + facet_wrap(~taxa_to_plot)

#plot_abundance(cdiff, c("42372","9710"), conc="Fluorescence")
#p <- plot_abundance(cdiff, c("42372","9710"), conc="Fluorescence")
#p + facet_wrap(~OTU)
#p + geom_point(aes(color=Treatment)) + facet_wrap(~OTU)