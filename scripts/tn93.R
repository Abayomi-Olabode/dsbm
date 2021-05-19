library(seqinr)
library(ape)
#library(strataG)
library(parallel)

setwd('/Users/abayomi/Desktop/new_protocol/')
mylist = c("200", "100", "50")

for (i in mylist) {
  
  for (j in seq(1, 10, 1)) {
    
    files <- list.files(paste0("./ref_rec_10/", i, "/", j, "/mafft/"))
    
    mclapply(files, function(f) {
      my_align <- read.alignment(file = paste0("./ref_rec_10/", i, "/", j, "/mafft/", f), format='fasta')
      y <- as.DNAbin(my_align)
      tn93 <- dist.dna(y, model = "TN93", variance = FALSE,
                       gamma = FALSE, pairwise.deletion = FALSE,
                       base.freq = NULL, as.matrix = FALSE)
      
      write.table(data.frame(as.matrix(tn93)), 
                  file = paste0("./ref_rec_10/", i, "/", j, "/tn93/", f, ".tn93.csv"),
                  append = F, 
                  quote = FALSE, 
                  sep = ",", 
                  col.names = NA)
    },
    mc.cores = 8
    )
  }
}
