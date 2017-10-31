############Generate Codon Sequences
#make sure to exclude stop codons
#be able to specify length
# number of sequences
install.packages("phangorn")


Codons <- c("TTT","TTC","TTA","TTG",
            "TCT","TCC","TCA","TCG",
            "TAT","TAC","TGT","TGC",
            "TGG","CTT","CTC","CTA",
            "CTG","CCT","CCC","CCA",
            "CCG","CAT","CAC","CAA",
            "CAG","CGT","CGC","CGA",
            "CGG","ATT","ATC","ATA",
            "ATG","ACT","ACC","ACA",
            "ACG","AAT","AAC","AAA",
            "AAG","AGT","AGC","AGA",
            "AGG","GTT","GTC","GTA",
            "GTG","GCT","GCC","GCA",
            "GCG","GAT","GAC","GAA",
            "GAG","GGT","GGC","GGA","GGG")

codon_freq <-c(
  0.07036383359263253,
  0.01731390948893221,
  0.02752186314732952,
  0.02683887396014298,
  0.03277754029233277,
  0.008065327554164121,
  0.01282049217861427,
  0.01250233575564504,
  0.02568910605290558,
  0.006321128829145523,
  0.01004794686512546,
  0.009798594594693457,
  0.0480003900261009,
  0.01181110189585198,
  0.01877470424600504,
  0.01830878666171929,
  0.04037071961144739,
  0.009933725177657846,
  0.01579046171273139,
  0.01539860181019696,
  0.01880586689400103,
  0.004627420943956824,
  0.007355660835935235,
  0.007173120984303035,
  0.01473893113235911,
  0.003626700061094885,
  0.005764933842450205,
  0.005621869855169922,
  0.0275398622849727,
  0.006776530762951509,
  0.01077184516823543,
  0.010504528463092,
  0.06367704098677834,
  0.01566853975795236,
  0.02490641453407825,
  0.0242883309498638,
  0.02966263590361272,
  0.00729886600851827,
  0.01160213939812641,
  0.01131421791131674,
  0.02324782740683193,
  0.005720421400950371,
  0.00909307369561169,
  0.008867417787828581,
  0.04343883280676464,
  0.01068867315950123,
  0.01699051274990056,
  0.01656887209167851,
  0.008547076168746835,
  0.01324911050085923,
  0.01618075534481809,
  0.003981479109348442,
  0.006328883909291891,
  0.006171824801270541,
  0.003120448885891702,
  0.004960206546687033,
  0.00483711286587539,
  0.02369557204540649,
  0.005830594621348161,
  0.009268202963569734,
  0.009038201005676032
)

num_files = 30

for( n in seq(from = 1, to = num_files, by =1)){
  num_seqs =5
  list_of_seqs <- list()
  for(i in seq(from = 1, to = num_seqs, by = 1)){
   seq <- sample(Codons,200,replace = TRUE, prob = codon_freq)
   seq_collapse <- paste(seq, collapse = "")
   list_of_seqs[[i]]<- seq_collapse
  }
  
  seq_names <- paste("seq", 1:num_seqs, sep = "_")
  
  file_name <- paste("R_codon_sim/Codon_sim_test_",n,".fasta", sep = "")
  write.fasta(list_of_seqs,seq_names,file_name )
  codon_dat <- read.dna(file_name, format = "fasta") 
  codon_PhyDat <- phyDat(codon_dat, type = "CODON", levels = NULL)
  dna_dist <- dist.ml(codon_PhyDat, model = "JC69")
  codon_NJ <- NJ(dna_dist)
  tree_name <- paste("R_codon_sim/Codon_sim_test_",n,".tree", sep = "")
  write.tree(codon_NJ, tree_name)
  }