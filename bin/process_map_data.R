#!/usr/bin/env Rscript

# this function takes reads_loc.txt as input and outputs tsv file
# with desired info

# import library
library(dplyr)
library(tidyr)

# setup variable to control command line args
args <- commandArgs(TRUE)

# read and process read_locs.txt
maps <- read.table(args[1])
#maps <- read.table("../read_locs.txt")

# remove columns that are not needed
maps <- maps[ , c(1:4)]

# cleanup first column
maps$V1 <- substring(maps$V1, 3)

# sgRNA reads are 20 bases adding this to starting position gives end position
maps <- mutate(maps, End = maps$V4 + 20)

# split V1 into 2 columns to compare
maps <- separate(maps, V1, c("a", "b"), sep = "\\|")

# check if sgRNA matches the target gene
maps$Match <- mapply(grepl, maps$b, maps$a)

# change colnames
colnames(maps) <- c("sgRNA", "Target_gene",
                    "FLAG(strand)",
                    "Chr", "Start",
                    "End", "Match")

#export table
write.table(
            maps,
            "./reads_mapped_table.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE
)
