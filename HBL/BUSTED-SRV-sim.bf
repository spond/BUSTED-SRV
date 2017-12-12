RequireVersion ("2.3.3");
LoadFunctionLibrary ("GrabBag");

SetDialogPrompt ("Load a BUSTED likelihood function fit (BUSTED.lf.bf) extension:");
ExecuteAFile (PROMPT_FOR_FILE);


fit_info = report_loaded_fit  ();


Export (lfstring, ^(lfname));
fprintf("/home3/sadie/BUSTED-SRV/data/LF1.txt", lfstring);


fprintf (stdout, "\n\nLoaded the following information\n\n", fit_info);

prompt_for_alphas (fit_info["alpha rate count"]);
prompt_for_omegas (fit_info["omega rate count"]);

sim_info = report_loaded_fit();

Export (lfstring,^(lfname));
fprintf("/home3/sadie/BUSTED-SRV/data/LF2.txt", lfstring);

fprintf (stdout, "\n\nUsing this information for the simulations\n\n", sim_info, "\n\n");

replicates = prompt_for_a_value ("How many replicates", 100, 1, 10000, 1);

SetDialogPrompt ("Set simulation information and repliates (with .xx extensions) to");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, sim_info);
sim_path = LAST_FILE_PATH;

DATA_FILE_PRINT_FORMAT = 9; // FASTA

GetString (tree, Tree, 0);
//fprintf(stdout, tree, "tree");
//GetInformation ( eqf, "busted.test_pi");
//fprintf(stdout, eqf);

charactersUniversalCode = {{"A","C","G","T"} {"3","TAA,TAG,TGA","",""}};


for (repl = 1; repl < replicates; repl += 1) {
    fprintf (stdout, "Simulating replicate ", repl, " / ", replicates, "\n");
   // DataSet simfile = SimulateDataSet (^(lfname));
    DataSet simfile = Simulate (^(tree), ^("busted.test_pi"), charactersUniversalCode, 1000);
    DataSetFilter simnucs = CreateFilter (simfile, 1);
    fprintf (sim_path + "_replicate." + repl, CLEAR_FILE, simnucs);
}

function report_loaded_fit () {
    GetInformation (omega.weights, "^busted.test.bsrel_mixture_aux_[0-9]+$");
    GetInformation (srv.weights,   "^busted.test.rv_gdd_weights_[0-9]+$");
    GetString (lfname, LikelihoodFunction,0);
    //fprintf (stdout, lfname, " lf ");    
    
    omega.rates = Columns (omega.weights) + 1;
    srv.rates   = Columns (srv.weights) + 1;
    
    omega.distribution = {omega.rates, 2};
    srv.distribution          = {srv.rates, 2};
    
    extract_distribution (omega.distribution, omega.rates, "busted.test.omega", "busted.test.bsrel_mixture_aux_", "1", );
    extract_distribution (srv.distribution, srv.rates, "busted.test.rv_gdd_rates_", "busted.test.rv_gdd_weights_", "busted.test.rv_gdd_norm");
   
    return { "sequences" : busted_srv.codon_data.species,
             "sites" : busted_srv.codon_data.sites, 
             "alpha rate count" : srv.rates,
             "alpha distribution" :  srv.distribution,
             "omega rate count" : omega.rates,
             "omega distribution" :  omega.distribution
           }
    
}

function extract_distribution (matrix, rates, rate_prefix, weight_prefix, normalizer) {
    left_over_weight = 1;
    for (k = 1; k <= rates; k += 1) {
        matrix [k-1][0] = Eval (rate_prefix + k + "/" + normalizer);
        if (k < omega.rates) {
            matrix [k-1][1] = left_over_weight * Eval (weight_prefix + k);
            left_over_weight = left_over_weight * (1-Eval (weight_prefix + k));
        } else {
            matrix [k-1][1] = left_over_weight;    
        }
    }

}

function prompt_for_alphas (rate_count) {
    leftover_weight = 1;
    current_rate     = 1; //alphas start at 1 now (this changes occasionally)
    
    fprintf (stdout, "\nSpecify the alpha distribution to be used for simulation (remember that it will be normalized to have mean 1)\n");
    fprintf (stdout, rate_count);
    while (current_rate <= rate_count && leftover_weight > 0) {
        Eval ("busted.test.rv_gdd_rates_" + current_rate + " = " + prompt_for_a_value ("alpha rate [" + current_rate + "]", 2*(current_rate) / rate_count, 0, 1000,  0));
        if (current_rate < rate_count) {
            new_weight = prompt_for_a_value ("alpha weight [" + current_rate + "]", leftover_weight / (rate_count - current_rate + 1), 0, leftover_weight,  0);
            Eval ("busted.test.rv_gdd_weights_" + current_rate + " = " + new_weight / leftover_weight);
            leftover_weight = leftover_weight - new_weight;
        }
        //prompt_for_a_value ("
        current_rate += 1;
    }       
}

function prompt_for_omegas (rate_count) {
    leftover_weight = 1;
    current_rate     = 1; //omegas start at 1
    
    fprintf (stdout, "\nSpecify the omega distribution to be used for simulation (remember that only the last omega value is permitted to be in the [1, infty) range)\n");
    
    while (current_rate <= rate_count && leftover_weight > 0) {
        if (current_rate < rate_count) {
            Eval ("busted.test.omega" + current_rate + " = " + prompt_for_a_value ("omega [" + current_rate + "]", (current_rate) / rate_count, 0, 1,  0));
            new_weight = prompt_for_a_value ("omega weight [" + current_rate + "]", leftover_weight / (rate_count - current_rate + 1 ), 0, leftover_weight,  0);
            Eval ("busted.test.bsrel_mixture_aux_" + (current_rate) + " = " + new_weight / leftover_weight);
            leftover_weight = leftover_weight - new_weight;
        } else {
			fprintf(stdout,current_rate);
            Eval ("busted.test.omega" + current_rate + " = " + prompt_for_a_value ("omega [" + current_rate + "]", 3, 1, 1000,  0));
        
        }
        //prompt_for_a_value ("
        current_rate += 1;
    }      
}