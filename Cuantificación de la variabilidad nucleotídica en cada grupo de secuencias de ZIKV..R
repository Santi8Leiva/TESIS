library(VariantAnnotation)
library(tidyverse)
vcf1 <- readVcf("/home/santi/Escritorio/ZIKA/Experimentos/4-Bioinformatica/15.norm/a.csq.n.cons.asian.fasta.african7.fastq.sorted.bam.vcf.gz")
vcf2 <- readVcf("/home/santi/Escritorio/ZIKA/Experimentos/4-Bioinformatica/15.norm/a.csq.n.cons.asian.fasta.asian.pre.fastq.sorted.bam.vcf.gz")
vcf3 <- readVcf("/home/santi/Escritorio/ZIKA/Experimentos/4-Bioinformatica/15.norm/a.csq.n.cons.asian.fasta.brain.7seq.fastq.sorted.bam.vcf.gz")
vcf4 <- readVcf("/home/santi/Escritorio/ZIKA/Experimentos/4-Bioinformatica/13.split.asian/concat.asian.epi.vcf")

epi <- rowRanges(vcf4) 
mcols(epi) <- cbind(mcols(epi), VariantAnnotation::info(vcf4))
epi$DV <- as.numeric(geno(vcf4)$DV)
epi.miss <- subset(epi, grepl("missense", epi$BCSQ))
count.epi.miss <- sum(as.numeric(epi.miss$DV))/558
epi.syno <- subset(epi, grepl("synonymous", epi$BCSQ))
count.epi.syno <- sum(as.numeric(epi.syno$DV))/558
count.epi.all <- sum(as.numeric(epi$DV))/558
r4 <- cbind("Asian Epidemic", count.epi.all, count.epi.miss, (count.epi.all-count.epi.miss))

br <- rowRanges(vcf3) 
mcols(br) <- cbind(mcols(br), VariantAnnotation::info(vcf3))
br$DV <- as.numeric(geno(vcf3)$DV)
br.miss <- subset(br, grepl("missense", br$BCSQ))
count.br.miss <- sum(as.numeric(br.miss$DV))/7
br.syno <- subset(br, grepl("synonymous", br$BCSQ))
count.br.syno <- sum(as.numeric(br.syno$DV))/7
count.br.all <- sum(as.numeric(br$DV))/7
r3 <- cbind("Brain", count.br.all, count.br.miss, (count.br.all-count.br.miss))

pre <- rowRanges(vcf2) 
mcols(pre) <- cbind(mcols(pre), VariantAnnotation::info(vcf2))
pre$DV <- as.numeric(geno(vcf2)$DV)
pre.miss <- subset(pre, grepl("missense", pre$BCSQ))
count.pre.miss <- sum(as.numeric(pre.miss$DV))/3
pre.syno <- subset(pre, grepl("synonymous", pre$BCSQ))
count.pre.syno <- sum(as.numeric(pre.syno$DV))/3
count.pre.all <- sum(as.numeric(pre$DV))/3
r2 <- cbind("Asian Pre Epidemic", count.pre.all, count.pre.miss, (count.pre.all-count.pre.miss))

af <- rowRanges(vcf1) 
mcols(af) <- cbind(mcols(af), VariantAnnotation::info(vcf1))
af$DV <- as.numeric(geno(vcf1)$DV)
af.miss <- subset(af, grepl("missense", af$BCSQ))
count.af.miss <- sum(as.numeric(af.miss$DV))/7
af.syno <- subset(af, grepl("synonymous", af$BCSQ))
count.af.syno <- sum(as.numeric(af.syno$DV))/7
count.af.all <- sum(as.numeric(af$DV))/7
r1 <- cbind("African", count.af.all, count.af.miss, (count.af.all-count.af.miss))

sum.data <- as_tibble(rbind(r1, r2, r3, r4))
sum.data <- dplyr::rename(sum.data, lineage = V1, total = count.af.all, missense = count.af.miss, syno =V4)
sum <- dplyr::select(sum.data, -total)
su <- pivot_longer(sum, cols = c(missense, syno), names_to = "snp.csq", values_to = "n")
su$n <- as.numeric(su$n)

pB <-  ggplot(su, aes(x= lineage, y=n, fill = snp.csq, colour= snp.csq)) + 
  geom_bar(stat = "identity", position = "stack", alpha = 0.8, width = 0.6) +
  scale_fill_manual("Tipo de Variación", labels= c("missense"="No Sinónima", "syno"= "Sinónima"), values = c("red2", "cyan2")) +
  scale_color_manual("Tipo de Variación", labels= c("missense"="No Sinónima", "syno"= "Sinónima"), values = c("red4", "cyan4")) +
  scale_y_continuous("Variaciones por secuencia") + 
  scale_x_discrete(limits =c("African", "Asian Pre Epidemic", "Brain", "Asian Epidemic"),
                   labels = c("Asian Pre Epidemic"= "Asian\nPre-Epidemic", "Asian Epidemic"= "Asian\nEpidemic")) +
  geom_text(aes(label=round(digit=1, n), y = c(1280, 600, 280, 112, 60, 15, 50, 15)),
            fontface= "bold", size=3, colour= "gray10")
      pB + theme_bw()  + theme (axis.text.x = element_text(face = "bold", size = 10),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(face = "bold.italic", size = 12),
                          axis.text.y = element_text(face = "bold"),
                          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                          legend.position = "bottom",
                          legend.text = element_text(size = 11, face = "bold"),
                          legend.title = element_text(face = "bold", size = 12))

      
      ggsave(filename = "bar.mix.600.png", plot=last_plot(), path = "/home/santi/Imágenes/", 
       dpi = 600, type = "cairo")
