# An example of a driver program for compute_statistics.R

source("compute_statistics.R")

# Place filenames (with full paths) into the file FILENAMES
d = scan("FILENAMES",what="character",quiet=TRUE)

## An example of a function which extracts galaxy ID number from its filename.
##
## If you use such a function, state id.extract=TRUE and
## id.extract.function=d2id below.
##
#d2id = function(name) {
#  splt = strsplit(name,"\\.")
#  splt = strsplit(splt[[1]],"_")
#  return(as.integer(splt[[1]][3]))
#}

# See the README.md file for descriptions of arguments
result = compute_statistics(filename.input=d,
         filename.increm="output_increm.Rdata",
         filename.output="output.Rdata",
         verbose=TRUE,delta.file=200,load=FALSE)

print(result$finish)
print(result$file.number)

q()
