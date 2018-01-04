############Generate Codon Sequences
#make sure to exclude stop codons
#be able to specify length
# number of sequences
install.packages("phangorn")
library("phangron")

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
  0.007018776800925262,
  0.02691762554590063,
  0.01618762252653746,
  0.009919614450785725,
  0.006259573964688113,
  0.02400601598217675,
  0.01443664948911825,
  0.008846635548760697,
  0.004845431924341126,
  0.01858265704223181,
  0.01117516986135365,
  0.006848033197273988,
  0.01267654409474975,
  0.04861566008390369,
  0.02923630663390879,
  0.01791571858670948,
  0.006389892528032194,
  0.02450579911950536,
  0.01473720723179166,
  0.009030814350967822,
  0.005698714468360446,
  0.02185507054906077,
  0.01314312184541552,
  0.008053974644669984,
  0.00441127993190415,
  0.01691764601625296,
  0.01017387166195978,
  0.006234447596796208,
  0.01154072236370785,
  0.04425968397720488,
  0.02661672576376049,
  0.01631046542418203,
  0.007014597672624323,
  0.02690159822178083,
  0.01617798407337897,
  0.009913708102334778,
  0.006255846881828861,
  0.02399172229203252,
  0.01442805360237588,
  0.008841368074663921,
  0.004842546851591323,
  0.0185715925348177,
  0.01116851592863004,
  0.006843955733329454,
  0.01266899621202216,
  0.04858671329065352,
  0.02921889871798068,
  0.01790505118857146,
  0.02667354891493221,
  0.009829667955634254,
  0.006202815030955515,
  0.02378834048572881,
  0.01430574460065521,
  0.008766418332109661,
  0.01841415806680862,
  0.01107383856801175,
  0.006785938386247646,
  0.01256159903134907,
  0.04817483566923342,
  0.02897120527075193,
  0.01775326710202835
)

num_files = 1

for( n in seq(from = 1, to = num_files, by =1)){
  num_seqs =31
  list_of_seqs <- list()
  for(i in seq(from = 1, to = num_seqs, by = 1)){
   seq <- sample(Codons,321,replace = TRUE, prob = codon_freq)
   seq_collapse <- paste(seq, collapse = "")
   list_of_seqs[[i]]<- seq_collapse
  }
  
  seq_names <- paste("seq", 1:num_seqs, sep = "_")
  
  file_name <- paste("simulations/Codon_sim_test_",n,".fasta", sep = "")
  write.fasta(list_of_seqs,seq_names,file_name )
  codon_dat <- read.dna(file_name, format = "fasta") 
  codon_PhyDat <- phyDat(codon_dat, type = "CODON", levels = NULL)
  dna_dist <- dist.ml(codon_PhyDat, model = "JC69")
  codon_NJ <- NJ(dna_dist)
  tree_name <- paste("simulations/Codon_sim_test_",n,".tree", sep = "")
  write.tree(codon_NJ, tree_name)
  }