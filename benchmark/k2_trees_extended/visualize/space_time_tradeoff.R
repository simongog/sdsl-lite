require(tikzDevice)
library(gdata)
library(xtable)
library(ggplot2)
library(reshape2)

source("../../basic_functions.R")
scale_data <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarking_results/scale2.log")