#Mapping read of IBD patients to vOTU. Input: csv.file, gff.file with gene annoatation and their position in each vOTU, depth.files with mapping reads of IBD patients to each vOTU. Output: plot with gene annotation, abundance and presence of selected vOTU. 
library(GenomicRanges)
library(rtracklayer)
library(gggenes)
annotation <- left_join(annotation,kegg, by="ko_id")
gff_file <- "genes.gff"
gff <- read.table(gff_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
fasta_file <- "reference_configs.fasta"
fasta_sequences <- readDNAStringSet(fasta_file)
contig_name <- "vOTU_242"
selected_contig <- fasta_sequences[names(fasta_sequences) == contig_name]
contig_granges <- GRanges(seqnames = names(selected_contig), ranges = IRanges(start = 1, end = width(selected_contig)))
gr <- makeGRangesFromDataFrame(gff, keep.extra.columns = TRUE)
selected_gr <- subsetByOverlaps(gr, contig_granges)
selected_gr <- as.data.frame(selected_gr)
selected_gr_gene <- left_join(selected_gr,annotation,by= "gene")
#gene_annotation
gene_annotation <- ggplot(selected_gr_gene,
                          aes(xmin = start, xmax = end, y = 1, fill = annotation, forward = orientation)) +
  geom_segment(aes(x = 0, xend = max(end), y = 1, yend = 1), color = "black") +
  geom_gene_arrow() +
  theme_void()
UC <- read.csv("UC_aligned.bam.depth", sep = "\t") 
UC$category <- "UC"
CD <- read.csv("CD_aligned.bam.depth", sep = "\t") 
CD$category <- "CD"
HY <- read.csv("HY_aligned.bam.depth", sep = "\t") 
HY$category <- "UC"
all <- bind_rows(CD,UC,HY)
all_1 <- all %>%
  mutate(read = case_when(
    category == "UC" ~ read / 65,
    category == "CD" ~ read / 37,
    category == "HY" ~ read / 53))
#select contig
all_sel <- all_1 %>% as_tibble() %>%
  filter(vOTU == contig_name)

# abundance plotting abundance and presence of UC, same to CD and HY.
uc_abundance <- ggplot(all_sel[all_sel$category == "UC", ], aes(x = position, y = read)) +
  geom_point(size = 1, color = "white") +
  geom_area(aes(y = read), fill = "#858585", alpha = 1) + 
  theme_bw()

uc_presence <- ggplot(all_sel[all_sel$category == "UC", ], aes(x = position, y = presence)) +
  geom_col(size = 1, fill = "#99CC00") +
  theme_bw()

