library(optparse)
library(dplyr)

option_list <- list(
  make_option(c("-d", "--datapath"), type='character', action='store', default='./',
              help="Path to where the CSV output files are stored"),
  make_option(c("-f", "--filename"), type='character', action='store', default='cell_filter_info.csv',
              help="Name of the files to search for"),
  make_option(c("-o", "--output"), type='character', action='store', default='Project_Cell_Filters.csv',
	      help="Name of the output file")
)
opt <- parse_args(OptionParser(option_list=option_list))

filenames <- Sys.glob(file.path(opt$datapath, '*', opt$filename))
filters <- read.csv(filenames[1]) 
filters$Sample <- basename(dirname(filenames[1]))

if(length(filenames) > 1) {
  for (filename in filenames[2:length(filenames)]) {
    data <- read.csv(filename)
    data$Sample <- basename(dirname(filename))
    filters <- full_join(filters, data)
  }
}

write.csv(filters[,c("Sample", setdiff(colnames(filters), "Sample"))], opt$output, row.names=FALSE)
