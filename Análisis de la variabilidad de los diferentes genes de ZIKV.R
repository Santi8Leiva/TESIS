library(tidyverse)
#estas primeras lineas son para importar y acamodar los datos.
# OJO! ajusten la direccion donde descargaron los 2 archivos de datos y ponganlo entre ""
  snp.af<-read.table("/home/santi/Escritorio/ZIKA/Experimentos/4-Bioinformatica/10.short/a.csq.cons.asian.fasta.african7.fastq.sorted.bam.vcf.gz",
                     sep="\t", header=F, blank.lines.skip=TRUE, comment.char = "#")
  snp.as<-read.table("/home/santi/Escritorio/ZIKA/Experimentos/4-Bioinformatica/10.short/a.csq.cons.asian.fasta.asian.pre.fastq.sorted.bam.vcf.gz",
                     sep="\t", header=F, blank.lines.skip=TRUE, comment.char = "#")
  snp.br<-read.table("/home/santi/Escritorio/ZIKA/Experimentos/4-Bioinformatica/10.short/a.csq.cons.asian.fasta.brain.7seq.fastq.sorted.bam.vcf.gz",
                     sep="\t", header=F, blank.lines.skip=TRUE, comment.char = "#")
  snpsm <- rbind(snp.af, snp.as, snp.br)
colnames(snpsm)<-c("chr","start","lineage","ref","alt","qual", "filter", "info", "format", "geno")

#crea columna para indicar el tipo de muatación.
snpsm$SNP.csq <- ifelse(grepl("missense", snpsm$info), "missense", "synonymous")

# creo una nueva columna con los nombres de los genes en los que ocurren las mutaciones.
snpsm <- snpsm %>% mutate(gen = case_when(between(start, 0, 107) ~ "5' UTR",
                                          between(start, 108, 473) ~ "C",
                                          between(start, 474, 977) ~ "prM",
                                          between(start, 978, 2489) ~ "E",
                                          between(start, 2490, 3545) ~ "NS1",
                                          between(start, 3546, 4223) ~ "NS2A",
                                          between(start, 4224, 4613) ~ "NS2B",
                                          between(start, 4614, 6464) ~ "NS3",
                                          between(start, 6465, 6845) ~ "NS4A",
                                          between(start, 6846, 6914) ~ "p2K",
                                          between(start, 6915, 7667) ~ "NS4B",
                                          between(start, 7668, 10376) ~ "NS5",
                                          between(start, 10376, 10808) ~ "3' UTR"))

# acá hace un resumen de las mutacione sy crea una nueva tabla donde cuenta cuantas mutaciones hay en cada gen.
c <- snpsm %>% group_by(gen, lineage, SNP.csq) %>% tally()

# acá crea la columna long que indica la longitud del gen 
# y despues calcula la cantidad de mutaciones cada 100 nucleotidos por gen.
c <- c %>% mutate(long = case_when(gen=="5' UTR"~ 107,
                                   gen=="C"~ 366,
                                   gen=="prM"~ 504,
                                   gen=="E"~ 1512,
                                   gen=="NS1"~ 1056,
                                   gen=="NS2A"~ 678,
                                   gen=="NS2B"~ 390,
                                   gen=="NS3"~ 1851,
                                   gen=="NS4A"~ 381,
                                   gen=="p2K"~ 69,
                                   gen=="NS4B"~ 753,
                                   gen=="NS5"~ 2708,
                                   gen=="3' UTR"~ 432)) %>% mutate(n.f.100 = n/(long/100))

#elimina los extremos no codificantes
c.c <- c [c(-1:-5),]
data.segm<-data.frame(x = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf ), 
                      y = c(4.84, 0.52, 1.84, 0.46, 17.52, 1.63), 
                      xend =c(Inf, Inf, Inf, Inf, Inf, Inf), 
                      yend = c(4.84, 0.52, 1.84, 0.46, 17.52, 1.63),
                      lineage=c("Asian Pre-Epidemic", "Asian Pre-Epidemic", "Brain", "Brain", "African", "African"),
                      colo = c("cyan4", "red4", "cyan4", "red4", "cyan4", "red4"))
# hace el grafico.
p0 <- ggplot(c.c, aes(x=gen, y=n.f.100, size=n, colour=SNP.csq)) +
  geom_point(alpha=0.7) +
  facet_wrap(vars(lineage), scales="free_y", ncol = 1, nrow = 3, strip.position = "right") +
  scale_x_discrete(name = "Genes de ZIKV", limits=c("C", "prM", "E", "NS1", "NS2A", "NS2B", "NS3", "NS4A", "p2K", "NS4B", "NS5")) +
  scale_y_continuous(name = "Tasa de Variaciones por cada 100 nucleótidos") +
  scale_color_manual(name = "Tipo de\nVariación", values = c("red2","cyan2"), labels= c("missense"="No Sinónima", "synonymous"= "Sinónima")) +
  scale_size_continuous(name = "Cantidad Total de\nVariaciones", breaks = c(5, 25, 50, 100, 200, 400))
p0 + geom_segment(data=data.segm, aes(x=x, y=y, xend=xend, yend=yend), colour=c("cyan4", "red4", "cyan4", "red4", "cyan4", "red4"),  linetype =3, inherit.aes = FALSE) +
  theme_bw() + 
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold.italic", size = 14),
        axis.title.y = element_text(face = "bold.italic", size = 14),
        axis.text.y = element_text(face = "bold", size= 12),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 12, face = "plain"),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_rect(fill = "gray80"))

ggsave(filename = "point.600.png", plot=last_plot(), path = "/home/santi/Escritorio/ZIKA/Experimentos/4-Bioinformatica/10.short/graph", 
       dpi = 600, type = "cairo") 