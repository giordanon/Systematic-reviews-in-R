library("tidyverse")
library("revtools")
library("readxl")
library("janitor")
library("bayestestR")
library("mi")
library("metafor")
library("multidplyr")
knitr::opts_chunk$set(echo = TRUE,  message = FALSE,  warning = FALSE,  tidy = TRUE)
# Function for cleaning doi's

doi_clean = function(data){
  data_w = data %>% 
    mutate(doi = str_replace(doi,pattern ="https://doi.org/", replacement = ""))
  return(data_w)
}


# Load bibliography file

load_bibliography <- function(path, # path to .ris files 
                              files #list of .ris files to read
                              ) {
  
  article_search = NULL
  for(i in 1:length(files)){
    # Read .ris file
    bibliography = read_bibliography(paste0(path, files[i])) %>% 
      dplyr::select(label, title, journal, year, doi, author, abstract, n1 )
    # Clean DOI URL's
    data = doi_clean(data = bibliography)
    # Row binding
    article_search = rbind(article_search, data)
  }
  return(article_search)
}


# Function for passing column names as argument of a another function

rename_fun = function(input_df,  new_cols) { 
  cols = "VAL"
  rename_(input_df, .dots = setNames(cols, new_cols)) 
}



# Multiple imputation function
multiple.imputation = function(n.imp, # Number of imputed datasets you want to generate
                               df.variables, # data.frame with the variables to consider the distribution from
                               impute.var, # name of the variable you want to imputate
                               var.name, # name of the output variable
                               seed = 1 # set seed for replicability of results
                               ){
  set.seed(seed)
  # TRANSFORMATIONS	
  data  =  df.variables %>% 
    rename(VAR_SD = impute.var) %>% 
    mutate(log_VAR =  log(VAR_SD)) %>% 
    unnest() %>% 
    #rename(log_VAR =  names(data[5])) %>% 
    select(-VAR_SD) %>% as.data.frame()
  # IMPUTATION
  mdf <- missing_data.frame(data) # warnings about missingness patterns
  mdf <- change(mdf, y = paste0(names(data[4])), what = "transformation", to = "identity")
  
  IMP<-	mi(mdf,max.minutes=60,n.iter= 20, n.chains = n.imp)
  round(mipply(IMP, mean, to.matrix = TRUE), 3)
  
  IMP.cmp <- mi::complete(y=IMP ,m=n.imp)
  lapply(IMP.cmp , summary)
  
  # Back Transformation
  data_VAR <- data.frame(VAR = exp(data$log_VAR))
  for (k in 1:n.imp) {# Second loop for imputed dataset 
    
    #Back-transformation
    data_VAR <- cbind(data_VAR, as.data.frame(exp(IMP.cmp[[k]]$log_VAR)))
    
  }# End of second loop
  
  n.imp = n.imp + 1
  #Calculate mean across all imputations
  df.iter = data_VAR[,2:n.imp] %>% 
    rowid_to_column("ID") %>% 
    pivot_longer(names_to = "VAR", 
                 values_to = "VAL", cols = -ID) %>% 
    
    #Remove extreme imputed values
    mutate(VAL = case_when(VAL > quantile(VAL, probs= 0.95) + 1.5 * IQR(VAL) ~ mean(VAL),
                           VAL < quantile(VAL, probs= 0.05) - 1.5 * IQR(VAL) ~ mean(VAL),
                           TRUE~VAL)) %>% 
    
    group_by(ID) %>% 
    summarise(VAL = mean(VAL))
  
  df.iter = rename_fun(df.iter, var.name)
  
  return(df.iter[,2])
}

# Function for calculating the pooled sample variance (treated + control)

pooled.var = function(sd.treated, # column that contains SD values of treated treatment
                      sd.control, # column that contains SD values of control treatment
                      n.control, # column that contains reps of control treatment
                      n.treated, # column that contains reps of treated treatment
                      m.treated, # column that contains mean values of treated treatment
                      m.control){# column that contains Smean values of control treatment
  var = ( (sd.treated^2)/(n.treated*(m.treated^2)) ) + ( (sd.control^2)/(n.control*(m.control^2)) ) 
  return(var)
}

# Function for back transforming effect sizes

trans = function(x){
  out = (exp(x)-1)*100
  return(out)
}

# Bootstrap meta-analytic model

bootstrap_rma = function(data, #input data
                         response_variable,  # name of the response variable in between quotations
                         moderator, # name of the moderator variable if any just type in NA
                         boot_num, # number of  bootstrap samples
                         cores = 16, # number of cores available in your computer
                         seed = 1) # seed for replicability
{
  set.seed(seed)
  cluster <- new_cluster(cores)
  cluster_library(cluster, c("metafor","tidyverse"))
  
  
  
  if (is.na(moderator)) {
    
    system.time(
      data %>%
        rename(#MOD = moderator, 
          RV = response_variable) %>%  
        #Run bootstrapped models
        modelr::bootstrap(n =  boot_num, id = 'boot_num') %>%
        group_by(boot_num) %>%
        partition(cluster) %>% 
        mutate(fit = strap %>% map(~rma(yi = RV,
                                        vi = VAR,
                                        weights = W,
                                        data = .x)
        )
        ) %>% 
        mutate(mod_val = fit %>% map(~data.frame(ESTIM = coef(.x),
                                                 PVAL = summary(.x)$pval)
        )
        )  %>%
        dplyr::select(-fit, -strap) %>%
        collect() %>% 
        unnest(mod_val) %>% 
        saveRDS(paste0("output/",response_variable,"_mod.RData"))
    )
    
  }
  
  if (!is.na(moderator)) {
    
    system.time(
      data %>%
        rename(MOD = moderator, 
               RV = response_variable) %>%  
        #Run bootstrapped models
        modelr::bootstrap(n =  boot_num, id = 'boot_num') %>%
        group_by(boot_num) %>%
        partition(cluster) %>% 
        mutate(fit = strap %>% map(~rma(yi = RV,
                                        vi = VAR,
                                        mods = ~ 0 + MOD, 
                                        weights = W,
                                        data = .x)
        )
        ) %>% 
        mutate(mod_val = fit %>% map(~data.frame(MOD = levels(as.factor(pluck(.x, "data")$MOD)),
                                                 ESTIM = coef(.x),
                                                 PVAL = summary(.x)$pval)
        )
        )  %>%
        dplyr::select(-fit, -strap) %>%
        collect() %>% 
        unnest(mod_val) %>% 
        saveRDS(paste0("output/",response_variable,"_",moderator,"_mod.RData"))
    )
  }
}

# Summarise bootstraps

summarise_bootstraps = function(data){  # data output as obteined from bootstrap_rma() function
  
  if ("MOD" %in% colnames(data)) {
    
    boot = data %>% 
      group_by(MOD) %>% 
      summarise_at(vars(ESTIM, PVAL), list(q975 = ~quantile(., 0.975, na.rm=T),
                                           q025 = ~quantile(., 0.025, na.rm=T),
                                           #hdi = ~hdi(.),
                                           q500 = ~quantile(.,0.500, na.rm=T))) %>% 
      unnest() %>% 
      #dplyr::select(MOD,ESTIM_mean,CILB_LB, CIUB_UB, CI_low3, CI_high3) %>% 
      ungroup()
    
  }
  
  if ("MOD" %nin% colnames(data)) {
    
    boot = data %>% 
      ungroup() %>%  
      summarise_at(vars(ESTIM, PVAL), list(q975 = ~quantile(., 0.975, na.rm=T),
                                           q025 = ~quantile(., 0.025, na.rm=T),
                                           #hdi = ~hdi(.),
                                           q500 = ~quantile(.,0.500, na.rm=T))) %>% 
      unnest() %>% 
      #dplyr::select(ESTIM_mean,CILB_LB, CIUB_UB, CI_low3, CI_high3) %>% 
      ungroup()
    
  }
  
  
  return(boot)
  
}


