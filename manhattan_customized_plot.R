# installing the packages required
library(igraph)
library(dplyr)
library(martini)
library(snpStats)

# Install the package if not already installed
if (!"IRanges" %in% installed.packages()) {
  install.packages("IRanges")
}

# Load the package
library(IRanges)

# load the files which required or created 
genotypes <- read.csv('/Users/Dell/Downloads/autism_vcf/genotypes.csv',row.names = 1)
genotypes <- as.matrix(genotypes)
genotypes <- new("SnpMatrix", genotypes)

fam <- read.csv('/Users/Dell/Downloads/autism_vcf/fam.csv')
fam <- fam[, -1]

map <- read.csv('/Users/Dell/Downloads/autism_vcf/map.csv')
map <- map[,-1]

# aggregating three files into one list
gwas <- list(map=map, genotypes=genotypes, fam=fam)
snp_map <- read.csv('/Users/Dell/Downloads/autism_vcf/snp_gene.csv')

ppi <- read.csv('/Users/Dell/Downloads/autism_vcf/gene_gene.csv') #this interaction was downloaded from string database.
gi_net <- get_GI_network(gwas,snpMapping = snp_map,ppi = ppi)
cones_gi <- scones.cv(gwas,gi_net)

plot_ideogram(gwas,cones_gi) # this was default function in Martini where you can still plot but cannot specify the varaints required in your plot





gwas2bed <- function(gwas) {
  
  map <- sanitize_map(gwas)
  bed <- subset(map, select = c("chr", "pos"))
  colnames(bed) <- c("chr", "start")
  bed$chr <- paste0("chr", bed$chr)
  bed$chr <- ifelse(bed$chr == "chr23", "chrX", bed$chr)
  bed$end <- bed$start
  
  return(bed)
  
}
snp_test <- function(gwas, covars, score, family, link) {
  
  genotypes <- gwas[['genotypes']]
  phenotypes <- gwas[['fam']][['affected']]
  
  if (score == 'chi2') {
    
    tests <- single.snp.tests(phenotypes, snp.data = genotypes)
    c <- chi.squared(tests, df=1)
    
  } else if (score == 'glm') {
    
    if (ncol(covars) && nrow(covars) == nrow(genotypes)) {
      covars <- as.matrix(covars)
      tests <- snp.rhs.tests(phenotypes ~ covars, snp.data = genotypes,
                             family = family, link = link)
    } else {
      tests <- snp.rhs.tests(phenotypes ~ 1, snp.data = genotypes,
                             family = family, link = link)
    }
    
    c <- chi.squared(tests)
    names(c) <- names(tests)
  } else if (score == 'r2') {
    
    c <- apply(genotypes, 2, cor, phenotypes, use="pairwise.complete.obs")
    c[is.na(c)] <- 0.01
    c <- c^2
    
  }
  
  c[is.na(c)] <- 0
  
  return(c)
}

sanitize_map <- function(gwas) {
  
  map <- gwas[["map"]]
  colnames(map) <- c("chr","snp","cm","pos","allele.1", "allele.2")
  
  return(map)
  
}

group_snps <- function(bed, chr_col, pos_col, threshold) {
  
  check_installed("IRanges", "group_snps")
  
  # make row numbers to reorder at the end
  bed[['id']]  <- 1:nrow(bed)
  
  ir_bed <- by(bed, bed[,chr_col], function(chr) {
    ir_chr <- IRanges::IRanges(start = chr[,pos_col], width = threshold)
    reduced_chr <- IRanges::reduce(ir_chr)
    ol <- as.matrix(IRanges::findOverlaps(ir_chr, reduced_chr))
    ir_chr[ol[,'queryHits'],] <- reduced_chr[ol[,'subjectHits'],]
    ir_chr <- as.data.frame(ir_chr)[,c('start','end')]
    ir_chr <- data.frame(chr_range   = chr[,chr_col],
                         start_range = ir_chr[['start']],
                         end_range   = ir_chr[['end']])
    unique(cbind(ir_chr, chr[,c(chr_col, pos_col)]))
    
  }) %>% do.call(rbind, .)
  
  bed <- merge(bed, ir_bed, by = c(chr_col, pos_col))
  bed[order(bed[['id']]), c('chr_range','start_range','end_range')]
  
}
check_installed <- function(pkgs, fn = "This function") {
  installed <- unlist(lapply(pkgs, requireNamespace, quietly = TRUE))
  if (!all(installed)) {
    stop(paste0(fn, " requires the following packages to be installed:\n", 
                paste(pkgs[!installed], collapse = '\n')), call. = FALSE)
}}




check_installed(c("circlize", "IRanges"))

# Initialize the circos plot with the ideogram of the 'hg19' species
circlize::circos.initializeWithIdeogram(species = 'hg19')

bed <- gwas2bed(gwas)
bed[['snp']] <- sanitize_map(gwas)[['snp']]
bed[['c']] <- snp_test(gwas, covars, 'chi2')
bed[['selected']] <- bed[['snp']] %in% names(V(cones_gi))
# bring the selected snps to the front
bed <- bed[with(bed, order(selected)),]

circlize::circos.genomicTrackPlotRegion(
  bed,
  ylim = c(0, 1.1 * max(bed$c, na.rm = TRUE)),
  panel.fun = function(region, value, ...) {
    # color according to selection/non-selection
    col = ifelse(value[[3]], "red3", "#bad7df")
    col[value$snp == 'rs35649919'] <- "red1"  # replace '#your_color' with the color you want for 'rsid123456'
    
    # assign different sizes and shapes based on color
    size = ifelse(col == "red3", 1.3, 0.9)  # larger size for '#e84545'
    size[value$snp == 'rs35649919'] <- 2  # replace '#your_size' with the size you want for 'rsid123456'
    shape = ifelse(col == "red3", 17, 16)  # different shape for '#e84545'
    shape[value$snp == 'rs35649919'] <- 15  # replace '#your_shape' with the shape you want for 'rsid123456'
      
      circlize::circos.genomicPoints(region, value, col = col,
                                     cex = size, pch = shape)
  }, 
  track.height = 0.5)


# Create edges
edges <- as_ids(E(cones_gi))
edges <- do.call(rbind, strsplit(edges, '|', fixed = TRUE))
edges <- as.data.frame(edges)

# Get positional information and remove self interactions
edges <- merge(edges, bed, by.x = 'V1', by.y = 'snp')
edges <- merge(edges, bed, by.x = 'V2', by.y = 'snp')
edges <- edges[edges[['chr.x']] != edges[['chr.y']],]

# Group into regions
edges <- cbind(group_snps(edges, "chr.x", "start.x", 50000),
               group_snps(edges, "chr.y", "start.y", 50000))
edges <- unique(edges)

x1 <- edges[,c(1,2,3)]
x2 <- edges[,c(4,5,6)]

# Add links between regions
circlize::circos.genomicLink(x1, x2, col=sample(1:5, nrow(x1), replace=TRUE))

# Clear the circos plot
circlize::circos.clear()


    
  



