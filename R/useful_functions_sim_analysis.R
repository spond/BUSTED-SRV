#potentially useful function from old analaysis
library(jsonlite)
library(ggplot2)
library(knitr)
library(stringr)
library(reshape2)
library(dplyr)
library(xtable)
library(tidyr)
library(readr)


pwr_tab <- function(dat) {
  
  
  A1.dat = dat
  A1.basic = A1.dat %>% group_by(True.omega.3.value, True.CV, Sites) %>% summarise(num_reps = n())
  
  A1.BUSTED = A1.dat %>% group_by(True.omega.3.value,True.CV, Sites) %>% filter(BUSTED.P <
                                                                               0.05) %>% tally()
  A1.BUSTED = rename(A1.BUSTED, "BUSTED_PWR" = n)
  
  A1.SRV = A1.dat %>%  group_by(True.omega.3.value, True.CV, Sites) %>% filter(BUSTED.SRV.P <
                                                                              0.05) %>%  tally()
  A1.SRV = rename(A1.SRV, "SRV_PWR" = n)
  
  A1.pwr.tab = full_join(A1.BUSTED, A1.SRV, by = c("True.omega.3.value", "True.CV", "Sites")) %>% 
    full_join(., A1.basic, by = c("True.omega.3.value", "True.CV", "Sites"))
  
  
  
  A1.means = A1.dat %>% group_by(True.omega.3.value, True.CV, Sites) %>%   summarise(
    "$BUSTED \\omega_3$ MLE" = mean(busted.omega.3.rate, na.rm = T),
    "SRV $\\omega_3$ MLE" = mean(srv.omega.3.rate, na.rm = T),
    "Mean CV" = mean(CV.SRV, na.rm = T)
  )
  A1.pwr.tab = full_join(A1.pwr.tab, A1.means, by = c("True.omega.3.value", "True.CV", "Sites"))
  A1.pwr.tab = replace(A1.pwr.tab, is.na(A1.pwr.tab), 0)
  A1.pwr.tab$BUSTED_PWR = A1.pwr.tab$BUSTED_PWR / A1.pwr.tab$num_reps
  A1.pwr.tab$SRV_PWR = A1.pwr.tab$SRV_PWR / A1.pwr.tab$num_reps
  return(A1.pwr.tab)
}


plot_pwr_lines <- function(pwr.tab){
  cl = colors()
  pwr.tab %>% ggplot(aes(color = factor(True.CV))) + geom_point(aes(x = True.omega.3.value, y = BUSTED_PWR)) + 
    geom_smooth(aes(x = True.omega.3.value, y = BUSTED_PWR))+ geom_point(aes(x = True.omega.3.value, y = SRV_PWR)) + 
    geom_smooth(aes(x = True.omega.3.value, y = SRV_PWR), linetype = 'dashed') + 
    labs(x = expression("True "*omega[3]), y = "Power", color = "CV of SRV")
}

plot_pwr_tile <- function(pwr.tab){
  temp=   pwr.tab %>% select(one_of(
    c("True.omega.3.value", "True.CV", "BUSTED_PWR", "SRV_PWR", "Sites")
  )) %>% filter(True.omega.3.value >= 1.1  && True.CV != 1.031) %>% melt(id.vars = c("True.omega.3.value", "True.CV", "Sites"))
  temp %>% ggplot(aes(y=True.omega.3.value, x = True.CV))+ geom_tile(aes(fill = value))+ facet_grid(Sites~variable)+
    labs(y =expression("True "*omega[3]), x = "True CV of SRV", color = "Power")
}

var = "omega.1"

mean_med_tab <- function(dat,var){
  
 
  omega.1.bias = dat %>% group_by(True.omega.3.value, True.CV) %>% summarise_at(vars(matches(var)), 
                                                                                funs(mean(.,na.rm=T),median(.,na.rm=T)))
 vars <- c(paste0("True.", var,".value_mean"), paste0("True.", var,".prop_mean"), paste0("True.", var, ".value_median"),
           paste0("True.", var, ".prop_median"))
   
  
  t <- omega.1.bias %>% melt(id.vars =c(vars, "True.omega.3.value","True.CV"))
  t  <- ungroup(t) %>% mutate( analysis = str_extract(t$variable,"srv\\.|busted\\.|True"))
  t = t %>% arrange(True.omega.3.value)
  t = mutate(t, 
             stat = str_extract(t$variable,"mean|median"), idk = str_extract(t$variable, "rate|prop"),
             thing = interaction(stat,idk))
  
  # glimpse(b)
  c = spread(select(t, -variable, - stat, -idk), key = thing, value = value)
 #  
 #  #rename some shit
 #  
 #  c =c %>% rename(True.value.Mean =True.omega.1.value_mean, True.Prop.Mean =True.omega.1.prop_mean, 
 #                  True.value.Median =True.omega.1.value_median,True.Prop.Median =True.omega.1.prop_median)
 #  c = c %>% select(True.omega.3.value, True.CV,analysis, mean.rate, True.value.Mean, median.rate,mean.prop,True.Prop.Mean,
 #                   median.prop)
 #  # cnames = c(expression("True "*omega[3]), "True CV of SRV", "Method", expression("Mean MLE "*omega[1]),expression("True "*omega[1]), expression("Median MLE "*omega[1]), expression("Mean Prop "*omega[1]), expression("True prop "*omega[1]), expression("Median prop "*omega[1])
 #  #            )
 #  
 #  cnames = c("True $\\omega_3$", "True CV \n of SRV", "Method", "Mean MLE $\\omega_1$", "True $\\omega_1$", "Median MLE $\\omega_1$", "Mean prop $\\omega_1$", "True prop $\\omega_1$", "Median prop $\\omega_1$"           )
 #  
 #  e= c %>% arrange(True.omega.3.value) %>% filter(True.omega.3.value %in% c(1.1,1.5,2,2.5,2.9))
 #  colnames(e) <- cnames
 #  add.to.row <- list(pos = list(0), command = NULL)
 #  command <- paste0("\\hline\n\\endhead\n",
 #                    "\\hline\n",
 #                    "\\multicolumn{", dim(e)[2] , "}{l}",
 #                    "{\\footnotesize Continued on next page}\n",
 #                    "\\endfoot\n",
 #                    "\\endlastfoot\n")
 #  add.to.row$command <- command
 #  
 #  e.xtab = e%>% xtable(caption = "Table for $\\omega_1$ values.") 
 #  align(e.xtab) = "r|lp{1cm}l|p{1.5cm}p{1.5cm}p{1.5cm}|p{1.5cm}p{1.5cm}p{1.5cm}|"
 # a.thing<-  e.xtab%>% print( include.rownames = F, sanitize.text.function = function(x){x},  add.to.row = add.to.row,tabular.environment = 'longtable', floating = FALSE, type = 'latex')
return(c)
}

mean_med_plot <- function(mean.med.tab, var){
  
  d = mean.med.tab %>% select(contains("value"), contains( "rate"), True.CV, analysis) %>% 
    melt(id.vars = c("True.omega.3.value","True.CV", "analysis"))
  
  d  %>%  ggplot(aes(x = True.omega.3.value, y = value,color = variable )) +
    geom_point() + geom_smooth()+ facet_grid(True.CV~analysis) +
    labs(x = expression("True "*omega[3]))+
    scale_color_discrete(name = "Measure", breaks = c("mean.rate", paste0("True.",var,".value_mean"), "median.rate"), 
                         labels = c("Mean MLE", "True value", "Median MLE"))
}
