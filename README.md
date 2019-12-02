# BUSTED[S]
A repo with data, scripts and information associated with the forthcoming BUSTED[S] publication. 

## HBL
includes **BUSTED-SRV-sim.bf** the script used to generate simulation replicates.

**batch_script.sh** an example of the slurm batch script used to run simulation replicates through BUSTED an BUSTED[S]

## Data

**to find the simulation data from the paper please visit: http://data.hyphy.org/web/busteds/**

includes original examples of data and results jsons from BUSTED[S]

**paper_LFs** includes likelihood function batch files (.LF.bf) used to generate the 16 sequence and 31 sequence simulation replicates in the paper as well as newick trees for both when omega = 1.

**BUSTED-SRV sim parameters - Inputs for sims.csv** shows the naming convention and parameter values used for each simulation of 100 replicants. Due to size limitations on github actual sequences will be hosted on silverback (TBD). 

# How to use BUSTED-SRV-sim.bf and analyze the data

Note: These small tutorials assume a wokring installation of [HyPhy](https://github.com/veg/hyphy) and a basic working knowledge of the HyPhy command line interface, R, and bash.

## Create your own simulated data set

1. Generate or use pre-existing Likelihood Function file (.LF.bf)
      * pre-existing can be found in the **paper_LFs** folder
2. Run `HYPHYMP BUSTED-SRV-sim.bf` You will be prompted for
      * .LF.bf
      * Synonymous rate categories and weights
      * Non-synonymous rate categories and weights
      * number of replicants 
      * the path and name you wish to save them to (eg: `/home/sadie/simulations/Sim_A` will give you files Sim_A_replicant_1 etc.)

## Use the pre-exisiting parameters

if you would like to use the paramters from the paper to regenerate the simulated data sets do the following:

1. Download the **BUSTED-SRV sim parameters - Inputs for sims.csv** file, the **paper_LFs** and **gen_batch_cmd_list.R**
2. using **gen_batch_cmd_list.R** change the relevant file paths 
3. Run `R gen_batch_cmd_list.R` to get a .txt file of bash commands
4. Run the bash commands to generate all of the replicates (may not want to run all at once)

## Analyzing the data

* There are several possible ways to run the analysis on a large number of files but the basic command used for this paper is:
> `HYPHYMP busted --alignment {input} --srv Yes --grid-size 2000 --starting-points 10 --output {output}`

* you can use bash, snakemake, or your favorite tool to loop over all the replicates in a folder and analyze them

* We've provided a snakemake example in **Snakefile** and an older bash example in **batch_script.sh**


