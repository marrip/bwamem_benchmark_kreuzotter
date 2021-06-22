# Load/install libraries

packages <- c('ggplot2', 'tidyr')
for(package in packages){
  if(!require(package, character.only = T)){
    install.packages(package)
  }
}

# Functions

get_data <- function(input){
  bm_data <- data.frame()
  for(bm in input){
    data <- read.delim(file = bm, sep = "\t")
    data$file <- bm
    bm_data <- rbind(bm_data, data)
  }
  bm_data <- separate(bm_data, "file", c("folder", "unit", "tool", "sample"), sep="/")
  bm_data$sample <- gsub(".tsv", "", bm_data$sample)
  bm_data <- separate(bm_data, "sample", c("sample", "threads"), sep="_")
  return(bm_data)
}

plot_data <- function(data, output){
  ggplot(data, aes(x=threads, y=s, color=tool)) +
    geom_point() +
    facet_grid(~sample) +
    labs(x="thread count", y="time [s]") +
    expand_limits(x=0, y=0)
  
  ggsave(snakemake@output[[1]])
}

# Run plotting workflow

sink(file = snakemake@log[[1]])
data <- get_data(snakemake@input)
plot_data(data, snakemake@output[[1]])
