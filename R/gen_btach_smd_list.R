library(readr)
library(ggplot2)
library(dplyr)
library(stringr)
library(plyr)

test <- read_csv("C:/Users/srwisots/Downloads/figuring out alphas - adding_3.csv")

#LF.path <- "/home3/sadie/BUSTED-SRV/Min_Sim_Codon.LF.bf"
LF.path <- "/home3/sadie/BUSTED-SRV/median_Codon.LF.bf"
#LF.path <- "/home3/sadie/BUSTED-SRV/16_seq.LF.bf"
#LF.path <- "/home3/sadie/BUSTED-SRV/lots_seq_Codon.LF.bf"
reps <- 100


#sim.path <- "/home3/sadie/BUSTED-SRV/o1_cv0/"
#sim.path <- "/home3/sadie/BUSTED-SRV/6_seq_tree/"
sim.path <- "/home3/sadie/BUSTED-SRV/32_seq_tree/"
#sim.path <- "/home3/sadie/BUSTED-SRV/16_seq_tree/"
#sim.path <- "/home3/sadie/BUSTED-SRV/63_seq_tree/"
  
lengths <- c("_100","_300", "_500", "", "_5000")

input.100 <- paste0( "sbatch --wrap=\"(" ,
       "echo ", LF.path, "; ",
       "echo ", test$`alpha 1`, "; ", 
       "echo ", test$`alpha 1 prop`, "; ", 
       "echo ", test$`alpha 2`, "; ",
       "echo ", test$`alpha 2 prop`, "; ",
       "echo ", test$`alpha 3`, "; ",
       "echo ", test$`omega 1`, "; ", 
       "echo ", test$`omega 1 prop`, "; ", 
       "echo ", test$`omega 2`, "; ",
       "echo ", test$`omega 2 prop`, "; ",
       "echo ", test$`omega 3`, "; ",
       "echo ", reps, "; ",
       "echo ", sim.path, test$sim_names, lengths[1], ")|~/bin/hyphy/HYPHYMP ~/BUSTED-SRV/HBL/BUSTED-SRV-sim",
       lengths[1], ".bf\"\n"
       )

input.300 <- paste0( "sbatch --wrap=\"(" ,
                     "echo ", LF.path, "; ",
                     "echo ", test$`alpha 1`, "; ", 
                     "echo ", test$`alpha 1 prop`, "; ", 
                     "echo ", test$`alpha 2`, "; ",
                     "echo ", test$`alpha 2 prop`, "; ",
                     "echo ", test$`alpha 3`, "; ",
                     "echo ", test$`omega 1`, "; ", 
                     "echo ", test$`omega 1 prop`, "; ", 
                     "echo ", test$`omega 2`, "; ",
                     "echo ", test$`omega 2 prop`, "; ",
                     "echo ", test$`omega 3`, "; ",
                     "echo ", reps, "; ",
                     "echo ", sim.path, test$sim_names, lengths[2], ")|~/bin/hyphy/HYPHYMP ~/BUSTED-SRV/HBL/BUSTED-SRV-sim",
                     lengths[2], ".bf\"\n"
)

input.500 <- paste0( "sbatch --wrap=\"(" ,
                     "echo ", LF.path, "; ",
                     "echo ", test$`alpha 1`, "; ", 
                     "echo ", test$`alpha 1 prop`, "; ", 
                     "echo ", test$`alpha 2`, "; ",
                     "echo ", test$`alpha 2 prop`, "; ",
                     "echo ", test$`alpha 3`, "; ",
                     "echo ", test$`omega 1`, "; ", 
                     "echo ", test$`omega 1 prop`, "; ", 
                     "echo ", test$`omega 2`, "; ",
                     "echo ", test$`omega 2 prop`, "; ",
                     "echo ", test$`omega 3`, "; ",
                     "echo ", reps, "; ",
                     "echo ", sim.path, test$sim_names, lengths[3], ")|~/bin/hyphy/HYPHYMP ~/BUSTED-SRV/HBL/BUSTED-SRV-sim",
                     lengths[3], ".bf\"\n"
)

input.1000 <- paste0( "sbatch --wrap=\"(" ,
                     "echo ", LF.path, "; ",
                     "echo ", test$`alpha 1`, "; ", 
                     "echo ", test$`alpha 1 prop`, "; ", 
                     "echo ", test$`alpha 2`, "; ",
                     "echo ", test$`alpha 2 prop`, "; ",
                     "echo ", test$`alpha 3`, "; ",
                     "echo ", test$`omega 1`, "; ", 
                     "echo ", test$`omega 1 prop`, "; ", 
                     "echo ", test$`omega 2`, "; ",
                     "echo ", test$`omega 2 prop`, "; ",
                     "echo ", test$`omega 3`, "; ",
                     "echo ", reps, "; ",
                     "echo ", sim.path, test$sim_names, lengths[4], ")|~/bin/hyphy/HYPHYMP ~/BUSTED-SRV/HBL/BUSTED-SRV-sim",
                     lengths[4], ".bf\"\n"
)

input.5000 <- paste0( "sbatch --wrap=\"(" ,
                      "echo ", LF.path, "; ",
                      "echo ", test$`alpha 1`, "; ", 
                      "echo ", test$`alpha 1 prop`, "; ", 
                      "echo ", test$`alpha 2`, "; ",
                      "echo ", test$`alpha 2 prop`, "; ",
                      "echo ", test$`alpha 3`, "; ",
                      "echo ", test$`omega 1`, "; ", 
                      "echo ", test$`omega 1 prop`, "; ", 
                      "echo ", test$`omega 2`, "; ",
                      "echo ", test$`omega 2 prop`, "; ",
                      "echo ", test$`omega 3`, "; ",
                      "echo ", reps, "; ",
                      "echo ", sim.path, test$sim_names, lengths[5], ")|~/bin/hyphy/HYPHYMP ~/BUSTED-SRV/HBL/BUSTED-SRV-sim",
                      lengths[5], ".bf\"\n"
)



#cat(input.5000, file = "32_seq_5000_o1_cmds.txt")
cat(input.100, input.300, input.500,input.1000, input.5000, file = "32_seq_add_3.txt")

