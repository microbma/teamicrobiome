# load packages
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(colorspace)
library(DESeq2)
source("taxa_info.R")

# read files #######
tea.otutable <- read.table("otutab.txt", 
                           sep = "\t",
                           row.names = 1,
                           header = TRUE) # OTU-table
tea.taxa <- read.table("otu.sintax",sep = "\t",row.names = 1) 
tea.taxa.df <- taxa.df(tea.taxa$V2)
row.names(tea.taxa.df) <- row.names(tea.taxa)

tea.taxa.df <- tea.taxa.df[-grep("Chloroplast",tea.taxa.df[,1]), ]


# metadata #######

site.info <- readxl::read_xlsx("soilmeta.xlsx",sheet = 1)


## Latitude 
dms = site.info$lat
gsub(pattern = "''",replacement = "\"",dms)
deg <- unlist(strsplit(x = dms,"°"))[seq(1,length(dms)*2,2)]
left <- unlist(strsplit(x = dms,"°"))[seq(2,length(dms)*2,2)]
minute <- unlist(strsplit(x = left,"'"))#[seq(1,length(dms)*2,2)]
minute1 <- minute[minute != ""][seq(1,length(dms)*2,2)]
left1 <- minute[minute != ""][seq(2,length(dms)*2,2)]
second <- gsub("\"","", left1)
min.dec <- as.numeric(minute1)+as.numeric(second)/60
lat <- as.numeric(deg)+as.numeric(min.dec)/60

dms = site.info$long
gsub(pattern = "''",replacement = "\"",dms)
deg <- unlist(strsplit(x = dms,"°"))[seq(1,length(dms)*2,2)]
left <- unlist(strsplit(x = dms,"°"))[seq(2,length(dms)*2,2)]
minute <- unlist(strsplit(x = left,"'"))#[seq(1,length(dms)*2,2)]
minute1 <- minute[minute != ""][seq(1,length(dms)*2,2)]
left1 <- minute[minute != ""][seq(2,length(dms)*2,2)]
second <- gsub("\"","", left1)
min.dec <- as.numeric(minute1)+as.numeric(second)/60
long <- as.numeric(deg)+as.numeric(min.dec)/60


# meta data 

tea.meta <- data.frame(site = paste("site",rep(c(1,2,4:46),each=5),sep="-"),
                       sample = rep(LETTERS[1:5],45),
                       provence = rep(site.info$provence,each=5),
                       latitude = rep(lat,each=5),
                       longitude = rep(long,each=5),
                       altitude = rep(site.info$alt, each =5))

colnames(tea.otutable) <- colnames(tea.otutable)
rownames(tea.meta) <- paste("X",rep(c(1,2,4:46),each=5),LETTERS[1:5],sep = "")


# generate phyloseq object #######
tea.phylo <- phyloseq(otu_table(tea.otutable,taxa_are_rows = TRUE),
                      tax_table(tea.taxa.df),
                      sample_data(tea.meta))


# normolization #######
tea.deseq <- phyloseq_to_deseq2(tea.phylo,design = ~sample)
tea.dds <- DESeq(tea.deseq,
                 fitType = "parametric",
                 test = "Wald")
tea.phylo.norm <- phyloseq(otu_table(counts(tea.dds,normalized = TRUE),
                                     taxa_are_rows = TRUE),
                           tax_table(tea.taxa.df),
                           sample_data(tea.meta))+
        
# ordination #######

tea.pcoa <- ordinate(subset_samples(tea.phylo.norm),method = "PCoA")
pcoa.df <- plot_ordination(tea.phylo.norm,ordination = tea.pcoa,color="sample")

pcoa.df1 <- pcoa.df$data %>%
        arrange(sample)

ggplot(pcoa.df1,
       aes(x=Axis.1, y = Axis.2,color=sample))+
        geom_point(size = 3,alpha=.4)+
        stat_ellipse(level = .6,size=1)+
        theme_bw()

tea.nmds <- ordinate(subset_samples(tea.phylo.norm),method = "NMDS",n=3)
nmds.df <- plot_ordination(tea.phylo.norm,ordination = tea.nmds,color="sample")

# calculate PERMANOVA

tb <- data.frame((t(otu_table(tea.phylo.norm))))
tb_env <- data.frame(sample_data(tea.phylo.norm))

tb.ano <- anosim(tb,tb_env$sample)
plot(tb.ano)

tb.ado <- adonis(tb~sample, tb_env)
tb.ado1 <- adonis(tb[136:225,]~sample, tb_env[136:225,])

nmds.df1 <- nmds.df$data %>%
        arrange(sample)

ggplot(nmds.df1,
       aes(x=NMDS1, y = NMDS2,color=sample))+
        geom_point(size = 2,alpha=.5)+
        stat_ellipse(level = .6, 
                     alpha = 1, 
                     size=.5,
                     linetype = 1, 
                     #color = "grey",
                     aes(group = sample))+
        theme_bw()+xlim(-3,1.5)+
        geom_vline(xintercept = 0,color = "grey",linetype = 2) + 
        geom_hline(yintercept = 0,color = "grey",linetype = 2) + 
        geom_text(x = 0.5,y=-1,label = "Stress=0.09",color = "black")+
        scale_color_manual(values = pal_npg()(5)[c(1,5,2:4)],
                           name = "",
                           breaks = LETTERS[1:5],
                           label = c("New leaf","Old leaf","Root","Rhizosphere","Bulk soil"))+
        theme( panel.grid = element_blank(),
               panel.background = element_blank(),
               aspect.ratio = 1,
               legend.position = c(.2,.9),
               legend.background = element_blank(),
               legend.key = element_rect(fill = NA),
               legend.key.height =  unit(6,"pt")
        )
ggsave("nmds.pdf",height = 4,width = 6)

# alpha-diversity #######

shannon <- diversity(t(otu_table(tea.phylo)),"shannon")
shannon.df <- data.frame(sample_data(tea.phylo),shannon)

shannon.df$sample <- factor(shannon.df$sample, 
                            labels = c("New leaf","Old leaf","Root","Rhizosphere","Bulk soil"))

require(ggbeeswarm)
ggplot(shannon.df, aes(x = sample,y=shannon,color=sample,fill = sample))+
        geom_violin(draw_quantiles = .5)+
        geom_quasirandom(size = 2, alpha = .5,shape = 16)+
        scale_fill_aaas(alpha = .2)+
        scale_color_manual(values = pal_npg()(5)[c(1,5,2:4)],
                           name = "Sample type",
                           breaks = LETTERS[1:5],
                           label = c("New leaf","Old leaf","Root","Rhizosphere","Bulk soil"))+
        #facet_grid(~sample)+
        ylim(c(2,8))+xlab("")+ylab("Shannon–Wiener index")+
        theme_bw()+guides(color = FALSE, fill = FALSE)+
        theme( panel.grid = element_blank(),
               panel.background = element_blank(),
               aspect.ratio = 1,
               legend.key.height =  unit(6,"pt"),
               axis.text.x = element_text(angle = -30,hjust = 0)
        )

ggsave("shannon.pdf",width = 4,height = 4)

# composition ######
source("make_taxa_bar.R")
tea.bar <- make_bar(level =1, 
                    keep = 40, 
                    phyloseq = tea.phylo.norm,
                    sample.class = 1:45)


tea.bar.sample <- apply(tea.bar, 
                        MARGIN = 1, 
                        FUN = function(x){tapply(x, 
                                                 INDEX = rep(LETTERS[1:5],45),
                                                 sum)})

tea.bar.sample.per <- t(t(tea.bar.sample)/rowSums(t(tea.bar.sample)))

tea.bar.sample.per1 <- data.frame(tea.bar.sample.per) %>% 
        gather(key = "class",value = "per") %>%
        add_column(sample = rep(LETTERS[1:5],41))        

tea.bar.sample.per1$class <- factor(tea.bar.sample.per1$class, 
                                    levels = names(sort(tea.bar.sample.per[2,-61],
                                                        decreasing = TRUE))[-56])

ggplot(tea.bar.sample.per1[!is.na(tea.bar.sample.per1$class),],aes(x = class, y = per*100, fill = sample))+
        geom_bar(stat = "identity",width = 1)+
        scale_fill_manual(values = pal_npg()(5)[c(1,5,2:4)],
                          name = "Sample type",
                          breaks = LETTERS[1:5],
                          label = c("New leaf","Old leaf","Root","Rhizosphere","Bulk soil"))+
        theme_bw()+
        ylab("Relative abundance (%)")+xlab("Phylum")+
        theme(aspect.ratio = .5,panel.grid = element_blank(),
              axis.text.x = element_text(angle = -70,hjust = 0,vjust = .5,size = 9))

ggsave("composition.pdf",height = 4,width = 8)

# nestedness #######
tea.nest <- make_bar(level = 1, 
                     keep = length(unique(tea.taxa.df[,1])), 
                     phyloseq = tea.phylo,
                     sample.class = 1:225)
tea.nest[tea.nest != 0] <- 1

tea.nest.df <- gather(data.frame(tea.nest),
                      key = "sample",
                      value = "present") %>%
        add_column(taxa = rep(row.names(tea.nest),225)) %>%
        add_column(pos = rep(rep(LETTERS[1:5],each = 60),45))

tea.nest.df$sample <- factor(tea.nest.df$sample, 
                             levels = colnames(tea.nest)[order(colSums(tea.nest),
                                                               decreasing = FALSE)][c(1:209,215)])
tea.nest.df$taxa <- factor(tea.nest.df$taxa, 
                           levels = rownames(tea.nest)[order(rowSums(tea.nest))])

ggplot(tea.nest.df[!is.na(tea.nest.df$sample),],#[grep(pattern = "C",tea.nest.df$sample),], 
       aes(x = sample,y=taxa,color = pos,size = factor(present)))+
        geom_point(shape=15)+
        xlab("225 samples\n(sorted by richness)")+
        ylab("60 phyla\n(sorted by prevalence)")+
        scale_size_discrete(range = c(0,.3),
                            limit = c(.5,1))+
        scale_color_manual(values = pal_npg()(5)[c(1,5,2:4)],
                           name = "Sample type",
                           breaks = LETTERS[1:5],
                           label = c("New leaf","Old leaf","Root","Rhizosphere","Bulk soil"))+
        #scale_color_manual(values = qualitative_hcl(225,alpha = .6)[sample(225)])+
        theme_minimal()+
        #theme_bw()+
        guides(size = FALSE)+
        theme(aspect.ratio = 60/225,
              panel.background = element_rect(fill = "white"),
              legend.position = "bottom",
              legend.key.width = unit(.6,units = "pt"),
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 8),
              axis.text = element_blank(),
              axis.title = element_text(size=8))
ggsave("nestedness.pdf",width = 5,height = 2)

sample.dis.df <- data.frame(ord = 1:225,
                            sample = rep(LETTERS[1:5],45)[order(colSums(tea.nest))])
ggplot(sample.dis.df[1:209,], aes(x = ord,
                                  fill = sample,
                                  color = sample))+
        geom_density(bw = "ucv")+
        scale_fill_manual(values = pal_npg(alpha = .1)(5)[c(1,5,2:4)],
                          breaks = LETTERS[1:5])+
        scale_color_manual(values = pal_npg(alpha = 1)(5)[c(1,5,2:4)],
                           breaks = LETTERS[1:5])+
        theme_void()+guides(fill=FALSE,color=FALSE)+
        theme(aspect.ratio = 1.5/10)

ggsave("nest_density.pdf",width = 5,height = 2)

## nestedness of random matrix with FALCON
nestedness.value.g <- function(tl = 1){
        tea.nest <- make_bar(level = tl, 
                             keep = length(unique(tea.taxa.df[,tl])), 
                             phyloseq = tea.phylo,
                             sample.class = 1:225)
        tea.nest[tea.nest != 0] <- 1
        
        tea.null <- permatfull(m = tea.nest,
                               fixedmar = "none",
                               shuffle = "both",
                               mtype = "prab",
                               times = 1)
        c(mean(mapply(tea.null$perm,
                      FUN = function(x=tea.null$perm[[1]]){
                              a <- nestednodf(x)
                              a$statistic[3]})),
          nestednodf(tea.null$orig)$statistic[3])
}

nodf.taxa <- NULL
for(i in 1:5){
        nodf.taxa <- c(nodf.taxa,nestedness.value.g(tl = i))
        print(i)
}

nestedness.value <- function(tl = 1, sl = 1){
        tea.nest <- make_bar(level = tl, 
                             keep = length(unique(tea.taxa.df[,tl])), 
                             phyloseq = tea.phylo,
                             sample.class = 1:225)
        tea.nest[tea.nest != 0] <- 1
        
        tea.null <- permatfull(m = tea.nest[,seq(sl,225,5)],
                               fixedmar = "none",
                               shuffle = "both",
                               mtype = "prab",
                               times = 1)
        c(mean(mapply(tea.null$perm,
                      FUN = function(x=tea.null$perm[[1]]){
                              a <- nestednodf(x)
                              a$statistic[3]})),
          nestednodf(tea.null$orig)$statistic[3])
}



nodf.Newleaf <- NULL
for(i in 1:5){
        nodf.Newleaf <- c(nodf.Newleaf,nestedness.value(tl = i,sl = 1))
        print(i)
}

nodf.Oldleaf <- NULL
for(i in 1:5){
        nodf.Oldleaf <- c(nodf.Oldleaf,nestedness.value(tl = i,sl = 2))
        print(i)
}

nodf.root <- NULL
for(i in 1:5){
        nodf.root <- c(nodf.root,nestedness.value(tl = i,sl = 3))
        print(i)
}

nodf.rhizo <- NULL
for(i in 1:5){
        nodf.rhizo <- c(nodf.rhizo,nestedness.value(tl = i,sl = 4))
        print(i)
}

nodf.bulk <- NULL
for(i in 1:5){
        nodf.bulk <- c(nodf.bulk,nestedness.value(tl = i,sl = 5))
        print(i)
}

nodf.whole <- data.frame(nodf = c(nodf.taxa,nodf.Newleaf,nodf.Oldleaf,nodf.root,nodf.rhizo,nodf.bulk),
                         mat = rep(c("null","exp"),5),
                         taxa = rep(c("Phylum","Class","Order","Family","Genus"),each=2),
                         group = rep(c("Whole dataset","New leaf","Old leaf","Root","Rhizosphere","Bulk soil")[c(1,6:2)],each = 10))
nodf.whole$taxa <- factor(nodf.whole$taxa, 
                          levels = c("Phylum","Class","Order","Family","Genus"))                     

nodf.whole <- nodf.whole %>%
        mutate(gl = paste(mat,group))


nodf.whole$group <- factor(nodf.whole$group,
                           c("Whole dataset","New leaf","Old leaf","Root","Rhizosphere","Bulk soil")[c(1,6:2)])

ggplot(nodf.whole,#[nodf.whole$group == "whole",], 
       aes(x = taxa, 
           y = nodf, 
           #alpha = mat, 
           color = mat,
           group = mat))+
        #geom_bar(stat = "identity",position = "dodge")+
        geom_line()+
        geom_point()+
        scale_color_aaas(alpha = .5,name = "Matrix",labels = c("Experiment","Null model"))+
        facet_wrap(~group)+
        #scale_color_aaas(alpha = .5)+
        theme_bw()+xlab("")+ylab("NODF")+
        theme(axis.text.x = element_text(angle = -60,hjust = 0),
              aspect.ratio = .5,
              panel.grid = element_blank(),
              legend.position = "right")

ggsave("nest_taxa_group.pdf",width = 4,height = 4)

# enriched plots #######
tea.deseq.de <- phyloseq_to_deseq2(subset_samples(phylo, sample %in% c("D", "E")),
                                   design = ~sample)
tea.dds.de <- DESeq(tea.deseq.de,
                    fitType = "parametric",
                    test = "Wald")
pv <- rep(0,25641)
pv[results(tea.dds.de)$pvalue < 0.05] <- 1
dds.de.df <- data.frame(results(tea.dds.de)[,],
                        p = pv)
ggplot(dds.de.df, aes(y = padj, x = -log2FoldChange, color = factor(p)))+
        #geom_bin2d(binwidth = c(0.5, 0.5))+
        geom_point(alpha=.1,size = 1)+
        xlab("Fold changed (log2)")+ylab("Adjusted P (-log10)")+
        scale_y_log10(breaks = c(1e-26,1e-18,1e-10,1e-2),
                      labels = expression(26,18,10,2))+
        scale_x_continuous(breaks = seq(-6,6,2))+
        geom_vline(xintercept = 0,linetype = 2)+
        scale_color_manual(values = c("grey","darkred"))+
        guides(color=FALSE)+
        theme_bw()+
        theme(panel.grid = element_blank(),
              axis.title = element_text(size=8),
              aspect.ratio = 1)

ggsave("dds_de.pdf",width = 1.5,height = 1.5)

tea.deseq.cd <- phyloseq_to_deseq2(subset_samples(phylo, sample %in% c("C", "D")),
                                   design = ~sample)
tea.dds.cd <- DESeq(tea.deseq.cd,
                    fitType = "parametric",
                    test = "Wald")
pv <- rep(0,25641)
pv[results(tea.dds.cd)$padj < 0.05] <- 1
dds.cd.df <- data.frame(results(tea.dds.cd)[,],
                        p = pv)

ggplot(dds.cd.df, aes(y = padj, x = -log2FoldChange, color = factor(p)))+
        #geom_bin2d(binwidth = c(0.5, 0.5))+
        geom_point(alpha=.1,size = 1)+
        xlab("Fold changed (log2)")+ylab("Adjusted P (-log10)")+
        scale_y_log10(breaks = c(1e-120,1e-100,1e-80,1e-60,1e-40,1e-20),
                      labels = seq(120,20,-20))+
        scale_x_continuous(breaks = seq(-10,10,2))+
        geom_vline(xintercept = 0,linetype = 2)+
        scale_color_manual(values = c("grey","darkred"))+
        guides(color=FALSE)+
        theme_bw()+
        theme(panel.grid = element_blank(),
              axis.title = element_text(size=8),
              aspect.ratio = 1)

ggsave("dds_cd.pdf",width = 1.5, height = 1.5)

tea.deseq.bc <- phyloseq_to_deseq2(subset_samples(phylo, sample %in% c("B", "C")),
                                   design = ~sample)
tea.dds.bc <- DESeq(tea.deseq.bc,
                    fitType = "parametric",
                    test = "Wald")
pv <- rep(0,25641)
pv[results(tea.dds.bc)$padj < 0.05] <- 1
dds.bc.df <- data.frame(results(tea.dds.bc)[,],
                        p = pv)
ggplot(dds.bc.df, aes(y = padj, x = -log2FoldChange, color = factor(p)))+
        #geom_bin2d(binwidth = c(0.5, 0.5))+
        geom_point(alpha=.1,size = 1)+
        xlab("Fold changed (log2)")+ylab("Adjusted P (-log10)")+
        scale_y_log10(breaks = c(1e-180,1e-150,1e-120,1e-90,1e-60,1e-30),
                      labels = seq(180,30,-30))+
        scale_x_continuous(breaks = seq(-10,10,2),limits=c(-10,10))+
        geom_vline(xintercept = 0,linetype = 2)+
        scale_color_manual(values = c("grey","darkred"))+
        guides(color=FALSE)+
        theme_bw()+
        theme(panel.grid = element_blank(),
              axis.title = element_text(size=8),
              aspect.ratio = 1)
ggsave("dds_bc.pdf",width = 1.5, height = 1.5)

tea.deseq.ab <- phyloseq_to_deseq2(subset_samples(phylo, sample %in% c("A", "B")),
                                   design = ~sample)
tea.dds.ab <- DESeq(tea.deseq.ab,
                    fitType = "parametric",
                    test = "Wald")
pv <- rep(0,25641)
pv[results(tea.dds.ab)$padj < 0.05] <- 1
dds.ab.df <- data.frame(results(tea.dds.ab)[,],
                        p = pv)
ggplot(dds.ab.df, aes(y = padj, x = -log2FoldChange, color = factor(p)))+
        #geom_bin2d(binwidth = c(0.5, 0.5))+
        geom_point(alpha=.1,size = 1)+
        xlab("Fold changed (log2)")+ylab("Adjusted P (-log10)")+
        scale_y_log10(breaks = c(1e-30,1e-25,1e-20,1e-15,1e-10,1e-5),
                      labels = seq(30,5,-5))+
        scale_x_continuous(breaks = seq(-10,10,2))+
        geom_vline(xintercept = 0,linetype = 2)+
        scale_color_manual(values = c("grey","darkred"))+
        guides(color=FALSE)+
        theme_bw()+
        theme(panel.grid = element_blank(),
              axis.title = element_text(size=8),
              aspect.ratio = 1)
ggsave("dds_ab.pdf",width = 1.5, height = 1.5)

tea.deseq.ac <- phyloseq_to_deseq2(subset_samples(phylo, sample %in% c("A", "C")),
                                   design = ~sample)
tea.dds.ac <- DESeq(tea.deseq.ac,
                    fitType = "parametric",
                    test = "Wald")
pv <- rep(0,25641)
pv[results(tea.dds.ac)$padj < 0.05] <- 1
dds.ac.df <- data.frame(results(tea.dds.ac)[,],
                        p = pv)
ggplot(dds.ac.df, aes(y = padj, x = -log2FoldChange, color = factor(p)))+
        #geom_bin2d(binwidth = c(0.5, 0.5))+
        geom_point(alpha=.1,size = 1)+
        xlab("Fold changed (log2)")+ylab("Adjusted P (-log10)")+
        scale_y_log10(breaks = c(1e-180,1e-150,1e-120,1e-90,1e-60,1e-30),
                      labels = seq(180,30,-30))+
        scale_x_continuous(breaks = seq(-10,10,2),limits=c(-10,10))+  
        geom_vline(xintercept = 0,linetype = 2)+
        scale_color_manual(values = c("grey","darkred"))+
        guides(color=FALSE)+
        theme_bw()+
        theme(panel.grid = element_blank(),
              axis.title = element_text(size=8),
              aspect.ratio = 1)
ggsave("dds_ac.pdf",width = 1.5, height = 1.5)

## manhatton graph #####

### rhizo V.S. soil
sign.otuname.de <- row.names(dds.de.df)[dds.de.df$pvalue < 0.05]
sign.otuname.de.na <- sign.otuname.de[!is.na(sign.otuname.de)]

dds.de <- data.frame(phylum = as.character(tea.taxa.df[sign.otuname.de.na,2]),
                     genus = as.character(tea.taxa.df[sign.otuname.de.na,5]),
                     fc = dds.de.df[sign.otuname.de.na,2],stringsAsFactors = FALSE)

dom.phyla <- names(sort(table(dds.de$phylum),decreasing = TRUE))[1:11]
dds.de[!as.character(dds.de$phylum) %in% dom.phyla,1] <- "Others"

dom.genus <- names(sort(table(dds.de$genus),decreasing = TRUE))[1:30]
dds.de[!as.character(dds.de$genus) %in% dom.genus,2] <- "Others"
dds.de$genus <- factor(dds.de$genus, levels = unique(dds.de$genus[order(dds.de$phylum)]))
dds.de$genus <- gsub("_genera_incertae_sedi", "", dds.de$genus)

set.seed(1234)
ggplot(dds.de[dds.de$genus != "Others",], 
       aes(x = genus, y = -fc, color = genus))+
        geom_quasirandom(alpha=.5,size = .6)+
        #geom_violin(fill=NA,width =.8)+
        geom_hline(yintercept = 0,linetype = 2)+
        xlab("")+ylab("Fold changed (log2)")+
        coord_flip()+
        scale_y_continuous(breaks = seq(-5,5,2),labels = seq(-5,5,2))+
        scale_color_manual(values = qualitative_hcl(31)[sample(31)])+
        theme_bw()+guides(color=FALSE)+
        theme(axis.line.x = element_blank(),
              aspect.ratio = 3,
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.text.x = element_text(size = 7),
              axis.text.y = element_text(size = 6,
                                         face = "italic"),
              axis.title = element_text(size = 7))
ggsave("enrich_de.pdf",height = 3,width = 2.5)

### root V.S. rhizosphere
sign.otuname.cd <- row.names(dds.cd.df)[dds.cd.df$p == 1]
sign.otuname.cd.na <- sign.otuname.cd[!is.na(sign.otuname.cd)]

dds.cd <- data.frame(phylum = as.character(tea.taxa.df[sign.otuname.cd.na,2]),
                     genus = as.character(tea.taxa.df[sign.otuname.cd.na,5]),
                     fc = dds.cd.df[sign.otuname.cd.na,2],
                     stringsAsFactors = FALSE)

dom.phyla.cd <- names(sort(table(dds.cd$phylum),decreasing = TRUE))[1:11]
dds.cd[!as.character(dds.cd$phylum) %in% dom.phyla.cd,1] <- "Others"

dom.genus.cd <- names(sort(table(dds.cd$genus),decreasing = TRUE))[1:30]
dds.cd[!as.character(dds.cd$genus) %in% dom.genus.cd,2] <- "Others"
dds.cd$genus <- factor(dds.cd$genus, 
                       levels = unique(dds.cd$genus[order(dds.cd$phylum)]))
dds.cd$genus <- gsub("_genera_incertae_sedi", "", dds.cd$genus)

set.seed(1234)
ggplot(dds.cd[dds.cd$genus != "Others",], 
       aes(x = genus, y = -fc, color = genus))+
        geom_quasirandom(alpha=.5,size = .6)+
        #geom_violin(fill=NA,width =.8)+
        geom_hline(yintercept = 0,linetype = 2)+
        xlab("")+ylab("Fold changed (log2)")+
        coord_flip()+
        scale_y_continuous(breaks = seq(-10,10,4),labels = seq(-10,10,4))+
        scale_color_manual(values = qualitative_hcl(31)[sample(31)])+
        theme_bw()+guides(color=FALSE)+
        theme(axis.line.x = element_blank(),
              aspect.ratio = 3,
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.text.x = element_text(size = 7),
              axis.text.y = element_text(size = 6,
                                         face = "italic"),
              axis.title = element_text(size = 7))
ggsave("enrich_cd.pdf",height = 3,width = 2.5)

### root V.S. old leaf
sign.otuname.bc <- row.names(dds.bc.df)[dds.bc.df$p == 1]
sign.otuname.bc.na <- sign.otuname.bc[!is.na(sign.otuname.bc)]

dds.bc <- data.frame(phylum = as.character(tea.taxa.df[sign.otuname.bc.na,2]),
                     genus = as.character(tea.taxa.df[sign.otuname.bc.na,5]),
                     fc = dds.bc.df[sign.otuname.bc.na,2],stringsAsFactors = FALSE)

dom.phyla.bc <- names(sort(table(dds.bc$phylum),decreasing = TRUE))[1:11]
dds.bc[!as.character(dds.bc$phylum) %in% dom.phyla.bc,1] <- "Others"

dom.genus.bc <- names(sort(table(dds.bc$genus),decreasing = TRUE))[1:30]
dds.bc[!as.character(dds.bc$genus) %in% dom.genus.bc, 2] <- "Others"
dds.bc$genus <- factor(dds.bc$genus, 
                       levels = unique(dds.bc$genus[order(dds.bc$phylum)]))
dds.bc$genus <- gsub("_genera_incertae_sedi", "", dds.bc$genus)

set.seed(1234)
ggplot(dds.bc[dds.bc$genus != "Others",], 
       aes(x = genus, y = -fc, color = genus))+
        geom_quasirandom(alpha=.5,size = .6)+
        #geom_violin(fill=NA,width =.8)+
        geom_hline(yintercept = 0,linetype = 2)+
        xlab("")+ylab("Fold changed (log2)")+
        coord_flip()+
        scale_y_continuous(breaks = seq(-10,10,4),labels = seq(-10,10,4))+
        scale_color_manual(values = qualitative_hcl(31)[sample(31)])+
        theme_bw()+guides(color=FALSE)+
        theme(axis.line.x = element_blank(),
              aspect.ratio = 3,
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.text.x = element_text(size = 7),
              axis.text.y = element_text(size = 6,
                                         face = "italic"),
              axis.title = element_text(size = 7))
ggsave("enrich_bc.pdf",height = 3,width = 2.5)

### old leaf V.S. new leaf
sign.otuname.ab <- row.names(dds.ab.df)[dds.ab.df$p == 1]
sign.otuname.ab.na <- sign.otuname.ab[!is.na(sign.otuname.ab)]

dds.ab <- data.frame(phylum = as.character(tea.taxa.df[sign.otuname.ab.na,2]),
                     genus = as.character(tea.taxa.df[sign.otuname.ab.na,5]),
                     fc = dds.ab.df[sign.otuname.ab.na,2],stringsAsFactors = FALSE)

dom.phyla.ab <- names(sort(table(dds.ab$phylum),decreasing = TRUE))[1:11]
dds.ab[!as.character(dds.ab$phylum) %in% dom.phyla.ab,1] <- "Others"

dom.genus.ab <- names(sort(table(dds.ab$genus),decreasing = TRUE))[1:30]
dds.ab[!as.character(dds.ab$genus) %in% dom.genus.ab, 2] <- "Others"
dds.ab$genus <- factor(dds.ab$genus, 
                       levels = unique(dds.ab$genus[order(dds.ab$phylum)]))
dds.ab$genus <- gsub("_genera_incertae_sedi", "", dds.ab$genus)

set.seed(1234)
ggplot(dds.ab[dds.ab$genus != "Others",], 
       aes(x = genus, y = -fc, color = genus))+
        geom_quasirandom(alpha=.5,size = .6)+
        #geom_violin(fill=NA,width =.8)+
        geom_hline(yintercept = 0,linetype = 2)+
        xlab("")+ylab("Fold changed (log2)")+
        coord_flip()+
        scale_y_continuous(breaks = seq(-5,5,2),labels = seq(-5,5,2))+
        scale_color_manual(values = qualitative_hcl(31)[sample(31)])+
        theme_bw()+guides(color=FALSE)+
        theme(axis.line.x = element_blank(),
              aspect.ratio = 3,
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.text.x = element_text(size = 7),
              axis.text.y = element_text(size = 6,
                                         face = "italic"),
              axis.title = element_text(size = 7))
ggsave("enrich_ab.pdf",height = 3,width = 2.5)

### root V.S. new leaf
sign.otuname.ac <- row.names(dds.ac.df)[dds.ac.df$p == 1]
sign.otuname.ac.na <- sign.otuname.ac[!is.na(sign.otuname.ac)]

dds.ac <- data.frame(phylum = as.character(tea.taxa.df[sign.otuname.ac.na,2]),
                     genus = as.character(tea.taxa.df[sign.otuname.ac.na,5]),
                     fc = dds.ac.df[sign.otuname.ac.na,2],stringsAsFactors = FALSE)

dom.phyla.ac <- names(sort(table(dds.ac$phylum),decreasing = TRUE))[1:11]
dds.ac[!as.character(dds.ac$phylum) %in% dom.phyla.ac,1] <- "Others"

dom.genus.ac <- names(sort(table(dds.ac$genus),decreasing = TRUE))[1:30]
dds.ac[!as.character(dds.ac$genus) %in% dom.genus.ac, 2] <- "Others"
dds.ac$genus <- factor(dds.ac$genus, 
                       levels = unique(dds.ac$genus[order(dds.ac$phylum)]))
dds.ac$genus <- gsub("_genera_incertae_sedi", "", dds.ac$genus)

set.seed(1234)
ggplot(dds.ac[dds.ac$genus != "Others",], 
       aes(x = genus, y = -fc, color = genus))+
        geom_quasirandom(alpha=.5,size = .6)+
        #geom_violin(fill=NA,width =.8)+
        geom_hline(yintercept = 0,linetype = 2)+
        xlab("")+ylab("Fold changed (log2)")+
        coord_flip()+
        scale_y_continuous(limits = c(-10,10),
                           breaks = seq(-10,10,4),
                           labels = seq(-10,10,4))+
        scale_color_manual(values = qualitative_hcl(31)[sample(31)])+
        theme_bw()+guides(color=FALSE)+
        theme(axis.line.x = element_blank(),
              aspect.ratio = 3,
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.text.x = element_text(size = 7),
              axis.text.y = element_text(size = 6,
                                         face = "italic"),
              axis.title = element_text(size = 7))
ggsave("enrich_ac.pdf",height = 3,width = 2.5)


# core-OTUs #######

comm.01 <- comm
comm.01[comm!=0] <- 1
comm.a <- comm.01[seq(1,225,5),]
comm.b <- comm.01[seq(2,225,5),]
comm.c <- comm.01[seq(3,225,5),]
comm.d <- comm.01[seq(4,225,5),]
comm.e <- comm.01[seq(5,225,5),]

otusum.a <- colSums(comm.a)
otusum.b <- colSums(comm.b)
otusum.c <- colSums(comm.c)
otusum.d <- colSums(comm.d)
otusum.e <- colSums(comm.e)

core.a <- names(otusum.a[otusum.a == 45])
core.b <- names(otusum.a[otusum.b == 45])
core.c <- names(otusum.a[otusum.c == 45])
core.d <- names(otusum.a[otusum.d == 45])
core.e <- names(otusum.a[otusum.e == 45])

core.otu <- c(core.a,core.b,core.c,core.d,core.e)
core.tb <- table(core.otu)
g.core.class <- tea.taxa.df[names(core.tb[core.tb == 5]),2]
a.core.class <- tea.taxa.df[core.a[core.a %in% names(core.tb[core.tb == 1])],2]
b.core.class <- tea.taxa.df[core.b[core.b %in% names(core.tb[core.tb == 1])],2]
c.core.class <- tea.taxa.df[core.c[core.c %in% names(core.tb[core.tb == 1])],2]
d.core.class <- tea.taxa.df[core.d[core.d %in% names(core.tb[core.tb == 1])],2]
e.core.class <- tea.taxa.df[core.e[core.e %in% names(core.tb[core.tb == 1])],2]

length(g.core.class)
length(table(g.core.class))
sort(table(g.core.class))

length(e.core.class)
length(table(e.core.class))
sort(table(e.core.class))
length(d.core.class)
length(table(d.core.class))
sort(table(d.core.class))
length(c.core.class)
length(table(c.core.class))
sort(table(c.core.class))
length(b.core.class)
length(table(b.core.class))
sort(table(b.core.class))
length(a.core.class)
length(table(a.core.class))
sort(table(a.core.class))

par(mfrow=c(2,3),mar = rep(0,4))
pie(table(g.core.class),
    lty = 0,
    labels = NA,
    col = rainbow(20,s=.9, v = .9, alpha = .7),
    radius = 1)
draw.circle(x = 0,y=0,
            radius = .5,
            col = "white",
            border = FALSE)

pie(table(a.core.class),
    lty = 0,
    labels = NA,
    col = rainbow(20,s=.9, v = .9, alpha = .7),
    radius = 1)
draw.circle(x = 0,y=0,
            radius = .5,
            col = "white",
            border = FALSE)

pie(table(b.core.class),lty = 0,
    labels = NA,
    col = rainbow(20,s=.9, v = .9, alpha = .7),
    radius = 1)
draw.circle(x = 0,y=0,
            radius = .5,
            col = "white",
            border = FALSE)

pie(table(c.core.class),
    lty = 0,
    labels = NA,
    col = rainbow(20,s=.9, v = .9, alpha = .7),
    radius = 1)
draw.circle(x = 0,y=0,
            radius = .5,
            col = "white",
            border = FALSE)

pie(table(d.core.class),
    lty = 0,
    labels = NA,
    col = rainbow(20,s=.9, v = .9, alpha = .7),
    radius = 1)
draw.circle(x = 0,y=0,
            radius = .5,
            col = "white",
            border = FALSE)

pie(table(e.core.class),
    lty = 0,
    labels = NA,
    col = rainbow(20,s=.9, v = .9, alpha = .7),
    radius = 1)
draw.circle(x = 0,y=0,
            radius = .5,
            col = "white",
            border = FALSE)

# distance-decay ####### 
## calculate geographic distance
require(geosphere)
geodist <- distm(x = sample_data(tea.phylo)[seq(1,225,5),5:4],
                 fun = distVincentyEllipsoid)
dist.bc.1 <- vegdist(t(otu_table(tea.phylo.norm))[seq(1,225,5),])
dist.bc.2 <- vegdist(t(otu_table(tea.phylo.norm))[seq(2,225,5),])
dist.bc.3 <- vegdist(t(otu_table(tea.phylo.norm))[seq(3,225,5),])
dist.bc.4 <- vegdist(t(otu_table(tea.phylo.norm))[seq(4,225,5),])
dist.bc.5 <- vegdist(t(otu_table(tea.phylo.norm))[seq(5,225,5),])

geo.decay.df <- data.frame(geo = c(geodist[lower.tri(geodist)]),
                           bc = c(dist.bc.1,dist.bc.2,dist.bc.3,dist.bc.4,dist.bc.5),
                           sample = rep(c("New leaf","Old leaf","Root","Rhizosphere","Bulk soil"),each= 990))

geo.decay.df$sample <- factor(geo.decay.df$sample,
                              levels = c("New leaf","Old leaf","Root","Rhizosphere","Bulk soil"))

ggplot(geo.decay.df, aes(x = geo/1000000, y = 1 - bc, color = 1-bc))+
        geom_point(size = .5,alpha = .2)+
        theme_bw() +
        scale_color_continuous_sequential(palette = "a",limits=c(0,.8),
                                          breaks = seq(0.1,.8,.2),
                                          name = "Bray-Curtis\nsimilarity")+
        geom_smooth(method = "lm", 
                    se = TRUE,
                    color="darkred",
                    size =.5, 
                    linetype =2)+
        #guides(color = FALSE)+
        xlab("Geographic distance (1000 km)")+
        ylab("Bray-Curtis similarity")+
        facet_grid(~sample)+
        theme(aspect.ratio = 1,
              axis.text = element_text(size=6),
              axis.title = element_text(size =7),
              strip.text = element_text(size =7),
              legend.title = element_text(size =7),
              legend.key.width = unit(3,"mm"),
              legend.text = element_text(size = 6),
              legend.key.height = unit(4,"mm"))

ggsave("dis_delay.pdf",width = 5,height = 2)

require(vegan)
corlog.1 <- mantel.correlog(D.eco = dist.bc.1,
                            D.geo = geodist)
corlog.2 <- mantel.correlog(D.eco = dist.bc.2,
                            D.geo = geodist)
corlog.3 <- mantel.correlog(D.eco = dist.bc.3,
                            D.geo = geodist)
corlog.4 <- mantel.correlog(D.eco = dist.bc.4,
                            D.geo = geodist)
corlog.5 <- mantel.correlog(D.eco = dist.bc.5,
                            D.geo = geodist)

corlog.df <- data.frame(x = c(corlog.1$mantel.res[1:5,1],
                              corlog.2$mantel.res[1:5,1],
                              corlog.3$mantel.res[1:5,1],
                              corlog.4$mantel.res[1:5,1],
                              corlog.5$mantel.res[1:5,1]),
                        y = c(corlog.1$mantel.res[1:5,3],
                              corlog.2$mantel.res[1:5,3],
                              corlog.3$mantel.res[1:5,3],
                              corlog.4$mantel.res[1:5,3],
                              corlog.5$mantel.res[1:5,3]),
                        p = c(corlog.1$mantel.res[1:5,5],
                              corlog.2$mantel.res[1:5,5],
                              corlog.3$mantel.res[1:5,5],
                              corlog.4$mantel.res[1:5,5],
                              corlog.5$mantel.res[1:5,5]),
                        sample = rep(c("Young leaf",
                                       "Old leaf",
                                       "Root",
                                       "Rhizosphere",
                                       "Bulk soil"),
                                     each = 5)
                        )
                
corlog.df$p1 <- 0
corlog.df$p1[corlog.df$p < 0.05] <- 1

corlog.df$sample <- factor(corlog.df$sample,
       levels = c("Young leaf","Old leaf","Root","Rhizosphere","Bulk soil"))


ggplot(corlog.df, aes(x = x, y = y))+
    geom_line(color = "black",size = 1)+
    geom_point(aes(color = factor(p1)))+
        facet_grid(~sample)+
    scale_x_continuous(name = "Distance (km)",
                       breaks = seq(10^5,10^6,2*10^5),
                       labels = seq(200,1000,200))+
    scale_color_manual(values = c("grey","#EE0000FF"))+
    ylab("Mantel r")+
    geom_hline(yintercept = 0,linetype = 2)+
    theme_bw()+guides(color = FALSE)+
    theme(panel.grid = element_blank(),
          aspect.ratio = 1,
          axis.text.x = element_text(angle = -90,
                                     hjust = 0,
                                     vjust = 0.5))

ggsave("corrlog.pdf",height = 3,width = 6)

# stioch #######
require(NST)
bulk.phylo <- subset_samples(tea.phylo.norm,sample == "E")
comm <- t(otu_table(tea.phylo.norm))
sample = as.matrix(sample_data(tea.phylo.norm))
nst <- tNST(comm,sample[,-c(1)],dist.method = "bray",rand = 1000)

ggplot(data = nst$index.grp,aes(x = group, y = NST.i.bray,fill = group))+
        geom_bar(stat = "identity")+
        geom_text(size = 3,aes(angle = seq(-72,-360,-72), 
                               label = paste(c("New leaf","Old leaf","Root","Rhizosphere","Bulk soil"),
                                             "\n(",round(NST.i.bray*100,1),"%)",sep="")))+
        geom_text(label = "NST",x=0,y=-.3,size = 3)+
        coord_polar(start = pi/6)+ylim(-.3,1)+
        scale_fill_manual(values = pal_npg(alpha = .7)(5)[c(1,5,2:4)])+
        theme_void()+guides(fill=FALSE)

ggsave("pie_nst.pdf",height = 4,width = 4)

# environmental driving #######

# source tracker #######

st.rhizo <- c(.042,.539,.419)
st.root <- c(.007,.017,0.016,.953)
st.old <- c(.140,.006,.849)
st.new <- c(.056,.004,.936)

par(mfcol=c(1,4),mar=rep(1,4),bg = "white")
pie(st.rhizo, 
    border = FALSE, cex=.6,
    label = c("Root","Bulk soil","Unknown"),
    col= c(pal_aaas(alpha = .5)(6)[1:2],grey(.7)))
draw.circle(x=0,y=0,radius = .4,col = "white",border = FALSE)
text(x=0,y=0,label = "Rhizosphere",cex=.7,font = 2)

pie(st.root,
    border = FALSE, cex=.6,
    label = c("New leaf","Old leaf","Rhizosphere","Unknown"),
    col= c(pal_aaas(alpha = .5)(6)[3:5],grey(.7)))
draw.circle(x=0,y=0,radius = .4,col = "white",border = FALSE)
text(x=0,y=0,label = "Root",cex=.7,font = 2)

pie(st.old,
    border = FALSE, cex=.6,
    label = c("New leaf","Root","Unknown"),
    col= c(pal_aaas(alpha = .5)(6)[c(3,1)],grey(.7)))
draw.circle(x=0,y=0,radius = .4,col = "white",border = FALSE)
text(x=0,y=0,label = "Old leaf",cex=.7,font = 2)

pie(st.new,
    border = FALSE, cex=.6,
    label = c("Old leaf","Root","Unknown"),
    col= c(pal_aaas(alpha = .5)(6)[c(4,1)],grey(.7)))
draw.circle(x=0,y=0,radius = .4,col = "white",border = FALSE)
text(x=0,y=0,label = "New leaf",cex=.7,font = 2)


# potential interactions #######
## prepare abundance table

### bulk soil
require(WGCNA)
require(RMThreshold)
require(igraph)
require(ggraph)

otus <- t(as.matrix(otu_table(tea.phylo.norm)))
make.net <- function(i){
        otus.e <- otus[seq(i,225,5),]
        otus.e <- otus.e[,order(colSums(otus.e),decreasing = TRUE)[1:1000]]
        cor.mat <- cor(otus.e,method = "spearman")
        sp.mat <- neten::Network_Enhancement(abs(cor.mat))
        net.mat <- sp.mat
        net.mat[sp.mat < 7] <- 0
        g.e <- graph_from_adjacency_matrix(net.mat,weighted = TRUE,mode = "undirected")
        edge.value <- cor.mat[lower.tri(cor.mat,diag = FALSE) & net.mat != 0]
        edge.dir <- rep("positive", length(E(g.e)))
        edge.dir[edge.value < 0] <- "negative"
        V(g.e)$phylum <- tea.taxa.df[colnames(otus.e),1]
        V(g.e)$class <- tea.taxa.df[colnames(otus.e),2]
        V(g.e)$order <- tea.taxa.df[colnames(otus.e),3]
        V(g.e)$family <- tea.taxa.df[colnames(otus.e),4]
        V(g.e)$genus <- tea.taxa.df[colnames(otus.e),5]
        E(g.e)$direction <- edge.dir
        return(g.e)
}

g.a <- make.net(i =1)
g.b <- make.net(i =2)
g.c <- make.net(i =3)
g.d <- make.net(i =4)
g.e <- make.net(i =5)

g.taxa.name <- c(unique(c(names(sort(table(V(g.a)$class),decreasing = TRUE)[1:6]),
                          names(sort(table(V(g.b)$class),decreasing = TRUE)[1:6]),
                          names(sort(table(V(g.c)$class),decreasing = TRUE)[1:6]),
                          names(sort(table(V(g.d)$class),decreasing = TRUE)[1:6]),
                          names(sort(table(V(g.e)$class),decreasing = TRUE)[1:6]))),
                 "Others")

g.a.t <- factor(V(g.a)$class,levels = g.taxa.name)
g.a.t[is.na(g.a.t)] <- "Others"
V(g.a)$taxa <- g.a.t

g.b.t <- factor(V(g.b)$class,levels = g.taxa.name)
g.b.t[is.na(g.b.t)] <- "Others"
V(g.b)$taxa <- g.b.t

g.c.t <- factor(V(g.c)$class,levels = g.taxa.name)
g.c.t[is.na(g.c.t)] <- "Others"
V(g.c)$taxa <- g.c.t

g.d.t <- factor(V(g.d)$class,levels = g.taxa.name)
g.d.t[is.na(g.d.t)] <- "Others"
V(g.d)$taxa <- g.d.t

g.e.t <- factor(V(g.e)$class,levels = g.taxa.name)
g.e.t[is.na(g.e.t)] <- "Others"
V(g.e)$taxa <- g.e.t

g.a.p <- subgraph.edges(g.a, E(g.a)[E(g.a)$direction != "negative"])
g.b.p <- subgraph.edges(g.b, E(g.b)[E(g.b)$direction != "negative"])
g.c.p <- subgraph.edges(g.c, E(g.c)[E(g.c)$direction != "negative"])
g.d.p <- subgraph.edges(g.d, E(g.d)[E(g.d)$direction != "negative"])
g.e.p <- subgraph.edges(g.e, E(g.e)[E(g.e)$direction != "negative"])

g.a.n <- subgraph.edges(g.a, E(g.a)[E(g.a)$direction == "negative"])
g.b.n <- subgraph.edges(g.b, E(g.b)[E(g.b)$direction == "negative"])
g.c.n <- subgraph.edges(g.c, E(g.c)[E(g.c)$direction == "negative"])
g.d.n <- subgraph.edges(g.d, E(g.d)[E(g.d)$direction == "negative"])
g.e.n <- subgraph.edges(g.e, E(g.e)[E(g.e)$direction == "negative"])

write_graph(g.a.p, "g_a_p.graphml","graphml")
write_graph(g.b.p, "g_b_p.graphml","graphml")
write_graph(g.c.p, "g_c_p.graphml","graphml")
write_graph(g.d.p, "g_d_p.graphml","graphml")
write_graph(g.e.p, "g_e_p.graphml","graphml")

write_graph(g.a.n, "g_a_n.graphml","graphml")
write_graph(g.b.n, "g_b_n.graphml","graphml")
write_graph(g.c.n, "g_c_n.graphml","graphml")
write_graph(g.d.n, "g_d_n.graphml","graphml")
write_graph(g.e.n, "g_e_n.graphml","graphml")

g.a.p.exp <- read_graph("g_a_p_export.graphml", format = "graphml")
g.b.p.exp <- read_graph("g_b_p_export.graphml", format = "graphml")
g.c.p.exp <- read_graph("g_c_p_export.graphml", format = "graphml")
g.d.p.exp <- read_graph("g_d_p_export.graphml", format = "graphml")
g.e.p.exp <- read_graph("g_e_p_export.graphml", format = "graphml")

### make graphics ####
ga <- ggraph(g.a.p.exp) +
        geom_edge_link(width = .05, 
                       color = "#00640066") +
        geom_node_point(aes(color = taxa), alpha = .9,size = .5) +
        scale_color_manual(values = c(qualitative_hcl(11,palette = "Dark 3"),grey(.3))[c(1,10,12,2:9)])  +
        scale_edge_width_continuous(range = c(.3,.8))+
        theme_void()+
        theme(legend.position = "none")
ggsave("gg_ap.pdf",ga,height = 4,width = 4)

gb <- ggraph(g.b.p.exp) +
        geom_edge_link(width = .05, 
                       color = "#00640066") +
        geom_node_point(aes(color = taxa), alpha = .9,size = .5) +
        scale_color_manual(values = c(qualitative_hcl(11,palette = "Dark 3"),grey(.3))[c(1,10,11,12,2:9)])  +
        scale_edge_width_continuous(range = c(.3,.8))+
        theme_void()+
        theme(legend.position = "none")
ggsave("gg_bp.pdf",gb,height = 4,width = 4)

gc <- ggraph(g.c.p.exp)+ 
        geom_edge_link(width = .05, 
                       color = "#00640066") +
        geom_node_point(aes(color = taxa), alpha = .9,size = .5) +
        scale_color_manual(values = c(qualitative_hcl(11,palette = "Dark 3"),grey(.3))[c(1,10,11,12,2:9)])  +
        scale_edge_width_continuous(range = c(.3,.8))+
        theme_void()+
        theme(legend.position = "none")
ggsave("gg_cp.pdf",gc,height = 4,width = 4)

gd <- ggraph(g.d.p.exp)+ 
        geom_edge_link(width = .05, 
                       color = "#00640066") +
        geom_node_point(aes(color = taxa), alpha = .9,size = .5) +
        scale_color_manual(values = c(qualitative_hcl(11,palette = "Dark 3"),grey(.3))[c(1,10,11,12,2:6,7:9)])  +
        scale_edge_width_continuous(range = c(.3,.8))+
        theme_void()+
        theme(legend.position = "none")
ggsave("gg_dp.pdf",gd,height = 4,width = 4)

ge <- ggraph(g.e.p.exp) +
        geom_edge_link(width = .05, 
                       color = "#00640066") +
        geom_node_point(aes(color = taxa), alpha = .9,size = .5) +
        scale_color_manual(values = c(qualitative_hcl(11,palette = "Dark 3"),grey(.3))[c(1,10,11,12,2:6,7:9)])  +
        scale_edge_width_continuous(range = c(.3,.8))+
        theme_void()+
        theme(legend.position = "none")
ggsave("gg_ep.pdf",ge,height = 4,width = 4)

ggraph(g.a.n) +
        geom_edge_link(width = 1, 
                       color = "darkred") +
        geom_node_point(aes(color = taxa), alpha = .9,size = 2) +
        scale_color_manual(values = c(qualitative_hcl(11,palette = "Dark 3"),grey(.3))[c(1,10,12,2:9)])  +
        scale_edge_width_continuous(range = c(.3,.8))+
        theme_void()+
        theme(legend.position = "none")
ggsave("gg_an.pdf",height = 2,width = 2)

ggraph(g.b.n) +
        geom_edge_link(width = 1, 
                       color = "darkred") +
        geom_node_point(aes(color = taxa), alpha = .9,size = 2) +
        scale_color_manual(values = c(qualitative_hcl(11,palette = "Dark 3"),grey(.3))[c(1,12,2:9)])  +
        scale_edge_width_continuous(range = c(.3,.8))+
        theme_void()+
        theme(legend.position = "none")
ggsave("gg_bn.pdf",height = 2,width = 2)

ggraph(g.c.n) +
        geom_edge_link(width = 1, 
                       color = "darkred") +
        geom_node_point(aes(color = taxa), alpha = .9,size = 2) +
        scale_color_manual(values =c(qualitative_hcl(11,palette = "Dark 3"),grey(.3))[c(1,11,12,2,4:8)])  +
        scale_edge_width_continuous(range = c(.3,.8))+
        theme_void()+
        theme(legend.position = "none")
ggsave("gg_cn.pdf",height = 2,width = 2)


ggraph(g.d.n) +
        geom_edge_link(width = 1, 
                       color = "darkred") +
        geom_node_point(aes(color = taxa), alpha = .9,size = 2) +
        scale_color_manual(values = c(qualitative_hcl(11,palette = "Dark 3"),grey(.3))[c(1,10,12,2,4,5,8,9)])  +
        scale_edge_width_continuous(range = c(.3,.8))+
        theme_void()+
        theme(legend.position = "none")
ggsave("gg_dn.pdf",height = 2,width = 2)


ggraph(g.e.n) +
        geom_edge_link(width = 1, 
                       color = "darkred") +
        geom_node_point(aes(color = taxa), alpha = .9,size = 2) +
        scale_color_manual(values = c(qualitative_hcl(11,palette = "Dark 3"),grey(.3))[c(1,10,11,12,2:5,8,9)])  +
        scale_edge_width_continuous(range = c(.3,.8))+
        theme_void()+
        theme(legend.position = "none")
ggsave("gg_en.pdf",height = 2,width = 2)




# calculate topological properties  ####
source("/Users/bma/Documents/Github/phyloseq_ext/network/netpro.R")
netpro.df <- rbind(netpro(g.a),netpro(g.b),netpro(g.c),netpro(g.d),netpro(g.e))

### circos ####
# library
library(circlize)
source("circlize_plot.R")
color <- c(qualitative_hcl(11,palette = "Dark 3",alpha = .3),grey(.3,alpha = .3))
graph = g.e
df1 <- maker.df1(graph)
df2 <- maker.df2(graph)
colors = color[as.numeric(as.character(df1$region))]
#circos.par(mar = c(0.5, 0.5, 0.5, 0.5), bg=rgb(1,1,1,0.1) )
factors = c(1:12)
circos.par(cell.padding = c(0, 0, 0, 0)) 
circos.initialize(factors=factor(df1$region,levels = df1$region), 
                  xlim = cbind(0, df1$xmax)) 
circos.trackPlotRegion(ylim = c(0, 1), 
                       track.height = 0.05, 
                       bg.col = colors, 
                       bg.border = NA ) 
sum1=rep(1,12)
sum2=rep(1,12)
linewidth = 20
for(k in 1:dim(df2)[1]){
        
        i<-match(df2$orig[k],df1$region)
        j<-match(df2$dest[k],df1$region)
        
        circos.link(sector.index1 = df1$region[i],
                    point1 = c(sum1[i],sum1[i]+df2$Freq[k]/linewidth),
                    sector.index2 = df1$region[j],
                    point2 = c(sum2[j],sum2[j]+df2$Freq[k]/linewidth),
                    col = colors[i],
                    h = 30)
        sum1[i]=sum1[i]+df2$Freq[k]/linewidth
        sum1[j]=sum1[j]+df2$Freq[k]/linewidth
        sum2[j]=sum2[j]+df2$Freq[k]/linewidth
        sum2[i]=sum2[i]+df2$Freq[k]/linewidth
}






