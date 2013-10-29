## root directory
rootPath <- "."

dataPath <- file.path(rootPath, "data")       # path to raw data
outputPath <- file.path(rootPath, "output")   # path to output direcotry
plotPath <- file.path(outputPath, "plots")    # path for plot files
csvPath <- file.path(outputPath, "csv")       # path for csv files

## create non-existing directories
dir.create(dataPath, showWarnings=FALSE, recursive=TRUE)
dir.create(plotPath, showWarnings=FALSE, recursive=TRUE)
dir.create(csvPath, showWarnings=FALSE, recursive=TRUE)
