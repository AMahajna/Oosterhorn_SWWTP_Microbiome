###############################################################################
# install packages

if( !require("corrplot") ){
  install.packages("corrplot")
  library("corrplot")
}

if( !require("plotly") ){
  install.packages("plotly", type = "source")
  library("plotly")
}

if( !require("ggcorrplot") ){
  install.packages("ggcorrplot")
  library("ggcorrplot")
}

if( !require("BiocManager") ){
  install.packages("BiocManager")
  library("BiocManager")
}

if( !require("Biostrings") ){
  BiocManager::install("Biostrings")
  library("Biostrings")
}

if( !require("GGally") ){
  install.packages("GGally")
  library("GGally")
}

if( !require("miaViz") ) {
  BiocManager::install("miaViz")
  library("miaViz")
}

if( !require("mia") ) {
  BiocManager::install("mia")
  library("mia")
}

if( !require("scater") ) {
  BiocManager::install("scater")
  library("scater")
}

if( !require("devtools") ) {
  install.packages("devtools")
  library("devtools")
}

if( !require("tidyverse") ) {
  install.packages("tidyverse")
  library("tidyverse")
}

if( !require("readxl") ) {
  install.packages("readxl")
  library("readxl")
}

if( !require("readr") ) {
  install.packages("readr")
  library("readr")
}

if( !require("ggplot2") ) {
  install.packages("ggplot2")
  library("ggplot2")
}

if( !require("forecast") ) {
  install.packages("forecast")
  library("forecast")
}

if( !require("FeatureTerminatoR") ){
  install.packages("FeatureTerminatoR")
  library("FeatureTerminatoR")
}

if( !require("DataExplorer") ){
  install.packages("DataExplorer")
  library("DataExplorer")
}

if( !require("caret") ){
  install.packages("caret")
  library("caret")
}

if( !require("mlbench") ){
  install.packages("mlbench")
  library("mlbench")
}

if( !require("randomForest") ){
  install.packages("randomForest")
  library("randomForest")
}

if( !require("psych") ){
  install.packages("psych")
  library("psych")
}

if( !require("heatmaply") ){
  install.packages("heatmaply")
  library("heatmaply")
}

if( !require("rio") ){
  install.packages("rio")
  library("rio")
}

if( !require("car") ){
  install.packages("car")
  library("car")
}

if( !require("FSelector") ){
  install.packages("FSelector")
  library("FSelector")
}

if( !require("data.table") ){
  install.packages("data.table")
  library("data.table")
}

if( !require("dplyr") ){
  install.packages("dplyr")
  library("dplyr")
}

if( !require("lubridate") ){
  install.packages("lubridate")
  library("lubridate")
}

if( !require("reshape2") ){
  install.packages("reshape2")
  library("reshape2")
}






