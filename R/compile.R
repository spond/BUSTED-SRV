####Newly formatted .jsons


library(jsonlite)

library(stringr)

library(dplyr)





compile <- function(dir,csv){
  #function reads BUSTED.JSONs and BUSTED_SRV.JSONS into a *.csv format
  #the csv can then be more easily used in R 
  srv.jsons <- list.files(path = dir,
                      pattern = '*BUSTED_SRV.json', recursive = TRUE) # list all of the json files in a given directory and the subdirectories in it
  busted.jsons <- list.files(path = dir,
                          pattern = '*BUSTED.json', recursive = TRUE)
  
  df.SRV<- NULL
  df.BUSTED <- NULL
  #create a table with 78 variables to fill 
  #step thru list of json files and read info from them
  #increments of two because want one line for each rep that includes BUSTED and BUSTED-SRV info
  for (i in  seq(from=1, to=length(srv.jsons), by=1)){
    filepath = paste(dir,srv.jsons[i], sep="") #file path of the current json
    
   test = filepath %>% readLines() %>% gsub(x=.,pattern="nan",replacement ='"NA"') %>% fromJSON() #read the JSON in
      #have to account for weird behavior caused by nan vs NA 
      
      FILE = test$input$`file name` #get name of file (useful for matching later)
      Sites = test$input$`number of sites` #get number of nucleotide sites
      
      tree_string = test$input$trees$`0` # get tree string
      
      
      Sequences = test$input$`number of sequences` #number of non Node named branch is the numb of seqs started with
      
      
      
      BUSTED.SRV.P = test$`test results`$`p-value`
      BUSTED.SRV.LR =test$`test results`$LRT
      BUSTED.SRV.AICc = test$fits$`Unconstrained model`$`AIC-c`
      
      #TO-DO: GET TREE LENGTH
      temp <- test$`branch attributes`$`0`
      temp <- temp %>% unlist() %>% t() %>% as.data.frame()
      
      unconstrained.bls <- temp %>% select(contains("uncon"))
      BUSTED.SRV.treelength <- unconstrained.bls %>% t() %>%   as.numeric() %>% sum()
      
      
      #get rates and weights
      temp <- test$fits$`Unconstrained model`$`Rate Distributions`$Test
      temp <- temp %>% unlist() %>% t() %>% as.data.frame() #turn into data.frame for easier manip
      
      
      srv.omega.rates <- temp %>% select(contains("omega"))
      srv.omega.props <- temp %>% select(contains("prop"))
      srv.alpha.rates <- temp %>% select(contains("SRV_rate"))
      srv.alpha.props <- temp %>% select(contains("weight"))
      names(srv.omega.rates) <- paste("srv.omega.rate",seq(1,length(srv.omega.rates)), sep = ".")
      names(srv.omega.props) <- paste("srv.omega.prop",seq(1,length(srv.omega.props)), sep = ".")
      names(srv.alpha.rates) <- paste("srv.alpha.rate",seq(1,length(srv.alpha.rates)), sep = ".")
      names(srv.alpha.props) <- paste("srv.alpha.prop",seq(1,length(srv.alpha.props)), sep = ".")
      
      
      mom2 = sum(srv.alpha.rates^2*srv.alpha.props)
      mean= sum(srv.alpha.rates*srv.alpha.props)
      CV.SRV = sqrt(mom2-mean^2)/mean
      
      
   
    x<- c(FILE,Sites, Sequences,  BUSTED.SRV.LR, CV.SRV,  BUSTED.SRV.P, BUSTED.SRV.AICc,
            BUSTED.SRV.treelength)
    x[2:length(x)] <- as.numeric(x[2:length(x)])
    names.SRV <- c("FILE", "Sites","Sequences","BUSTED.SRV.LR","CV.SRV", "BUSTED.SRV.P", "BUSTED.SRV.AICc",
                   "BUSTED.SRV.treelength")
    names(x)<- names.SRV
    df.SRV <-rbind.data.frame(df.SRV, c(x, srv.omega.rates, srv.omega.props,srv.alpha.rates,srv.alpha.props), 
                              stringsAsFactors = F)
  }   
  for (i in  seq(from=1, to=length(busted.jsons), by=1)){ 
   
      filepath = paste(dir,busted.jsons[i], sep="")
      test = filepath %>% readLines() %>% gsub(x=.,pattern="nan",replacement ='"NA"') %>% fromJSON() #read the JSON in
      #have to account for weird behavior caused by nan vs NA 
      
      FILE = test$input$`file name` #get name of file (useful for matching later)
      Sites = test$input$`number of sites` #get number of nucleotide sites
      
      tree_string = test$input$trees$`0` # get tree string
      
      
      Sequences = test$input$`number of sequences` #number of non Node named branch is the numb of seqs started with
      
      
      
      BUSTED.P = test$`test results`$`p-value`
      BUSTED.LR =test$`test results`$LRT
      BUSTED.AICc = test$fits$`Unconstrained model`$`AIC-c`
      
      #TO-DO: GET TREE LENGTH
      temp <- test$`branch attributes`$`0`
      temp <- temp %>% unlist() %>% t() %>% as.data.frame()
      
      unconstrained.bls <- temp %>% select(contains("uncon"))
      BUSTED.treelength <- unconstrained.bls %>% t() %>%   as.numeric() %>% sum()
      
      
      #BUSTED.SRV.treelength = test$fits$`Unconstrained model`$`tree length`
      
      
      #get rates and weights
      temp <- test$fits$`Unconstrained model`$`Rate Distributions`$Test
      temp <- temp %>% unlist() %>% t() %>% as.data.frame() #turn into data.frame for easier manip
      
      
      busted.omega.rates <- temp %>% select(contains("omega"))
      busted.omega.props <- temp %>% select(contains("prop"))
      names(busted.omega.rates) <- paste("busted.omega.rate",seq(1,length(busted.omega.rates)), sep = ".")
      names(busted.omega.props) <- paste("busted.omega.prop",seq(1,length(busted.omega.props)), sep = ".")
      
      
      
    
    
    #print(FILE)
    x<- c(FILE, Sites, Sequences, BUSTED.LR, BUSTED.P, BUSTED.AICc,BUSTED.treelength)
    x[2:length(x)] <- as.numeric(x[2:length(x)])
    names.BUSTED <- c("FILE","Sites","Sequences","BUSTED.LR","BUSTED.P","BUSTED.AICc", "BUSTED.treelength")
    names(x) <- names.BUSTED
    df.BUSTED <-rbind.data.frame(df.BUSTED, c(x, busted.omega.rates,busted.omega.props), stringsAsFactors = F)
    
    
  }

  temp <- full_join(df.BUSTED,df.SRV,by=c("FILE","Sites","Sequences"))
  write.csv(file = csv, x = temp, row.names= F)
  
  #return(as.data.frame(df,stringAsFactors = FALSE))
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
    num_rates <- sapply(x1, function (x) str_detect(x,"rate count")) %>% which(.==TRUE) %>% names(.) %>%
      str_extract(.,"[0-9]+") %>% as.numeric()
    
    
    
    
    end_1 = 6+(num_rates[1]-2)
    start_2 = end_1+5
    end_2 = start_2+num_rates[2]-2
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



