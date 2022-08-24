#=======================#
# 
#  Genomics Exercise
#  
#  Part 1
#=======================#


## Part 1 - Calculate gene-gene correlations for all genes in the file  & suppress any warnings

#libraries
library(tidyverse) # data manipulation and plottin
library(data.table) # to save data frame as data.table
options(scipen=50)



# 1a. load data 

nrv<-read.csv(file="GTEX.nerve_samples.PMC5359387.txt",sep="")


# 1b. transpose data to make genes column name

t_nrv<-nrv %>% 
  select(-Description) %>% 
  t(.) %>% 
  data.frame()

col_nm<-t_nrv[1,]

# 1c. Fix variable types and set sample IDs as a column

nrv_cln<-t_nrv %>% 
  set_names(., nm=col_nm) %>% 
  .[-1,] %>% 
  mutate_if(is.character,as.numeric) %>% 
  mutate(samp=rownames(.))

setDT(nrv_cln)

# 1d. Create a matrix of all possible combination of genes

Mx=combn(names(nrv_cln)[1:618],2)
dim(Mx)

# 1e. Apply function to calculate Pearson Correlation and pvalue

corr_prilm<-nrv_cln[ ,map(1:ncol(Mx),function(x) c(pearson=ifelse(sum(eval(as.name(Mx[1,x])))>0 & sum(eval(as.name(Mx[2,x]))>0),
                                                                  cor.test(eval(as.name(Mx[1,x])),eval(as.name(Mx[2,x])),method="pearson")$estimate,"Col Sum 0!"),
                                                   pvalue=ifelse(sum(eval(as.name(Mx[1,x])))>0 & sum(eval(as.name(Mx[2,x]))>0),
                                                                 cor.test(eval(as.name(Mx[1,x])),eval(as.name(Mx[2,x])),method="pearson")$p.value,"Col Sum 0!")))]


#1f. Transpose and cbind matrix and create numeric columns for correlation coefficient and p-value

corr_intrm<-as.data.frame(cbind(t(Mx),t(corr_prilm))) %>% 
  set_names(c("g_1","g_2","crr_est","crr_p")) %>% 
  mutate(pearson_corr=as.numeric(crr_est),
         pvalue=as.numeric(crr_p))


# 1g. Finalize correlation df 
corr_fin<-corr_intrm %>% select(-crr_est,-crr_p)


#1h Examine correlations with error messages - ensure estimates and plaues are NA.
na_corr=filter(corr_intrm, is.na(pearson_corr))
head(na_corr)

