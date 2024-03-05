suppressMessages(library(bc3net))

# read the input from the Matlab file
args = commandArgs(trailingOnly=TRUE)

# read in the generated tmp file
data = as.matrix(read.csv(args[1], header = FALSE, sep = ","))

# run bnet 
bnet <- bc3net(data, igraph = FALSE)
# write the output to a file that is then feed to matlab 
write.table(bnet, paste(args[2],'bnet.csv', sep="/"), col.names = FALSE, row.names = FALSE, sep = "\t")
