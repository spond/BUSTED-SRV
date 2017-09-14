####Newly formatted .jsons

library(jsonlite)

library(stringr)

library(dplyr)

test <-fromJSON("../data/sim_replicate.6.BUSTED.json")

filepath <- "../data/sim_replicate.6.BUSTED.json"

library(jsonlite)

library(stringr)

library(dplyr)



compile <- function(cur.dir,csv){
  jsons <- list.files(path = cur.dir,
                      pattern = '*.json', recursive = TRUE) # list all of the json files in a given directory and the subdirectories in it
  
  
  df<- NULL
  #create a table with 78 variables to fill 
  #step thru list of json files and read info from them
  #increments of two because want one line for each rep that includes BUSTED and BUSTED-SRV info
  for (i in  seq(from=1, to=length(jsons), by=2)){
    filepath = paste(cur.dir,jsons[i], sep="") #file path of the current json
    
    test = filepath %>% readLines() %>% gsub(x=.,pattern="nan",replacement ='"NA"') %>% fromJSON() #read the JSON in
    #have to account for weird behavior caused by nan vs NA 
    
    FILE = test$input$`file name` #get name of file (useful for matching later)
    Sites = test$input$`number of sites` #get number of nucleotide sites
    
    tree_string = test$input$trees$`0` # get tree string
    #from that can get branch names and lengths:
    #x= tree_string %>% str_replace_all("\\(",":") %>% str_replace_all("\\)",":") %>%     str_replace_all(",",":") %>% str_split(":")
    #x= unlist(x)
    #x =x[x !=""]
    #br_len = matrix(x,ncol = 2,  byrow = T)
    #colnames(br_len) = c("Branch", "Length")
    
    
    Sequences = test$input$`number of sequences` #number of non Node named branch is the numb of seqs started with
    
    #now we get to the trickier part
    #reading in BUSTED and BUSTED-SRV results
    #check if jsons[i] and jsons[i+1] are matched
    # j.i = str_extract(jsons[i],".*?(?=\\.BUSTED)" )
    # j.ii = str_extract(jsons[i+1],".*?(?=\\.BUSTED)" )
    # if(j.i == j.ii){
    #   js =jsons[c(i,i+1)]
    #   j.srv = js[which(str_detect(js,"SRV"))]
    #   j.busted <- js[which(!str_detect(js,"SRV"))]
    # }    else {
    #   print(c("Issue with:", FILE))
    # }
    # 
    #getting BUSTED-SRV results for a file
    filepath = paste(cur.dir,j.srv, sep="")
    
    test = filepath %>% readLines() %>% gsub(x=.,pattern="nan",replacement ='"NA"') %>% fromJSON() #read in json     
    
    
    BUSTED.SRV.P = test$`test results`$p
    BUSTED.SRV.LR =test$`test results`$LR
    BUSTED.SRV.AICc = test$fits$`Unconstrained model`$`AIC-c`
    BUSTED.SRV.treelength = test$fits$`Unconstrained model`$`tree length`
    num_omega_rate = length(test$fits$`Unconstrained model`$`rate distributions`$FG[,1])
    num_alpha_rate = length(test$fits$`Unconstrained model`$`rate distributions`$SRV[,1])
    #OMEGA values for BUSTED.SRV
    srv.omega.rates = as.numeric(test$fits$`Unconstrained model`$`rate distributions`$FG[,1])
    names(srv.omega.rates) = paste0("BUSTED.SRV.omega",seq(from = 1, to =num_omega_rate,by=1),".MLE")
    srv.omega.props = as.numeric(test$fits$`Unconstrained model`$`rate distributions`$FG[,2])
    names(srv.omega.props) = paste0("BUSTED.SRV.omega",seq(from = 1, to =num_omega_rate,by=1),".prop")
    #ALPHA values for BUSTED.SRV
    srv.alpha.rates = as.numeric(test$fits$`Unconstrained model`$`rate distributions`$SRV[,1])
    names(srv.alpha.rates) = paste0("BUSTED.SRV.alpha",seq(from = 1, to =num_alpha_rate,by=1),".MLE")
    srv.alpha.props = as.numeric(test$fits$`Unconstrained model`$`rate distributions`$SRV[,2])
    names(srv.alpha.props) = paste0("BUSTED.SRV.alpha",seq(from = 1, to =num_alpha_rate,by=1),".prop")
    
    mom2 = sum(srv.alpha.rates^2*srv.alpha.props)
    mean= sum(srv.alpha.rates*srv.alpha.props)
    CV.SRV = sqrt(mom2-mean^2)/mean
    
    
    mom2 = sum(srv.alpha.rates^2*srv.alpha.props)
    mean= sum(srv.alpha.rates*srv.alpha.props)
    CV.SRV = sqrt(mom2-mean^2)/mean
    
    
    
    filepath = paste(cur.dir,j.busted, sep="")
    
    test = filepath %>% readLines() %>% gsub(x=.,pattern="nan",replacement ='"NA"') %>% fromJSON()
    #print(filepath)
    BUSTED.P = test$`test results`$p
    BUSTED.LR = test$`test results`$LR
    BUSTED.AICc = test$fits$`Unconstrained model`$`AIC-c`
    BUSTED.treelength = test$fits$`Unconstrained model`$`tree length`
    num_omega_rate = length(test$fits$`Unconstrained model`$`rate distributions`$FG[,1])
    
    #OMEGA values for BUSTED
    busted.omega.rates = as.numeric(test$fits$`Unconstrained model`$`rate distributions`$FG[,1])
    names(busted.omega.rates) = paste0("BUSTED.omega",seq(from = 1, to =num_omega_rate,by=1),".MLE")
    busted.omega.props =as.numeric( test$fits$`Unconstrained model`$`rate distributions`$FG[,2])
    names(busted.omega.props) = paste0("BUSTED.omega",seq(from = 1, to =num_omega_rate,by=1),".prop")
    
    
    
    
    #print(FILE)
    x<- c(FILE, BUSTED.LR, BUSTED.SRV.LR, CV.SRV, BUSTED.P, BUSTED.SRV.P,BUSTED.AICc,BUSTED.SRV.AICc,
          BUSTED.treelength ,BUSTED.SRV.treelength, Sites, Sequences)
    x[2:length(x)] <- as.numeric(x[2:length(x)])
    names(x) <- c("FILE", "BUSTED.LR","BUSTED.SRV.LR","CV.SRV", "BUSTED.P","BUSTED.SRV.P","BUSTED.AICc","BUSTED.SRV.AICc",
                  "BUSTED.treelength", "BUSTED.SRV.treelength","Sites","Sequences")
    df <-rbind(df, c(x, busted.omega.rates , busted.omega.props,
                     srv.omega.rates, srv.omega.props,srv.alpha.rates,srv.alpha.props))
    
    
  }
  
  write.csv(file = csv, x = df, row.names= F)
  
  # return(as.data.frame(df,stringAsFactors = FALSE))
}

#can't mix and match rate categories yet
simulation_inputs <- function(dir,csv){
  require("stringr")
  require("jsonlite")
  require("dplyr")
  list = list.files(path = dir, recursive = T, pattern ="^([^.]+)$")
  #set up the empty data frame
  
  setup.tab <- NULL
  #loop thru each file to get info in correct format
  for(i in seq(from = 1, to= length(list))){
    x=readLines(paste(dir,list[i], sep = "/"))
    #making this a readable json
    
    x1 = x[2:(length(x)-1)]  %>% gsub(x=.,pattern="\\{",replacement ='\\[') %>% gsub(x=.,pattern ="\\}", replacement ="\\]")
    x1 = c(x[1],x1,x[length(x)])
    num_rates=as.numeric(str_extract(x1[4], "[0-9]+"))
    
    end_1 = 6+(num_rates-2)
    start_2 = end_1+5
    end_2 = start_2+num_rates-2
    x1[c(6:end_1,start_2:end_2)]  = x1[c(6:end_1,start_2:end_2)] %>%  gsub(x=.,pattern ="\\]", replacement ="\\],")
    
    r= fromJSON(x1)
    omega_rates = r$`omega distribution`[,1]
    names(omega_rates)= paste0("True.omega",seq(from=1, to = r$`omega rate count`, by =1), ".value" )
    omega_weights = r$`omega distribution`[,2]
    names(omega_weights) = paste0("Ture.omega",seq(from=1, to = r$`omega rate count`, by =1), ".prop" )
    
    Alpha_rates = r$`alpha distribution`[,1]
    names(Alpha_rates)= paste0("True.alpha",seq(from=1, to = r$`alpha rate count`, by =1), ".value" )
    Alpha_weights = r$`alpha distribution`[,2]
    names(Alpha_weights) = paste0("True.alpha",seq(from=1, to = r$`alpha rate count`, by =1), ".prop" )
    
    mom2 = sum(Alpha_rates^2*Alpha_weights)
    
    mean = sum(Alpha_rates*Alpha_weights)
    
    CV.SRV = sqrt(mom2-mean^2)/mean
    x<- c(r$sites,list[i],CV.SRV)
    names(x)<-c("Sites","FILE","True.CV")
    
    setup.tab = rbind(setup.tab,c(x, omega_rates,omega_weights,Alpha_rates,Alpha_weights))
  }
  
  
  #return(as.data.frame(setup.tab,stringsAsFactors = FALSE))
  write.csv(file = csv, x = setup.tab, row.names= F)
}

add_truth <- function(results, truth){
  
  n <- truth %>% select(-FILE, -Sites) %>% colnames()
  results[,n]<-NA
  for( i in seq(from = 1, to = nrow(truth), by = 1)){
    temp = which(str_detect(results$FILE,truth$FILE[i]))
    results[temp,n]<-truth[i,3:length(truth)]
  }
  return(results)
}




process_dat <- function(dir, basename){
  require("dplyr")
  temp = paste(dir,basename,"_Truth.csv", sep = "")
  simulation_inputs(dir,temp)
  truth = read.csv(temp, as.is = T)
  temp = paste(dir,basename,"_results.csv", sep = "")
  compile(dir,temp)
  dat = read.csv(temp, as.is = T)
  dat = add_truth(dat, truth)
  #dat = mutate(dat, Cat = str_extract(dat$FILE, "YesYes|YesNo|NoNo|NoYes"))
  dat$True.CV = round(dat$True.CV, 3)
  write.csv(file = paste(dir,basename,"_processed.csv", sep = ""),x = dat, row.names = F)
}



dir= "final_results/1000_Codons/"
basename = "All_5_8_17"
all.dat = process_dat(dir, basename)