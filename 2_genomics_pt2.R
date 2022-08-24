#=======================#
# 
#  Genomics Exercise
#  Part 2
#=======================#

##  Part 2. Find Top 10 most correlated genes &
##          Create a function that takes the input data and can return either correlated or anti-correlated gene pairs, 
##          for an arbitrary number of pairs


#libraries


library(tidyverse) # data manipulation and plottin
library(data.table) # to save data frame as data.table
options(scipen=50)


#2a. Most Anti correlated 10 gene pair
(acorr10<-corr_fin %>%
    slice_min(pearson_corr,n=10))


# 2b. Function

GetTopCorr<-function(df,gene1,gene2,corr,npair){
  
  CorrPair<-  df %>%
    dplyr:: slice_max({{corr}},n=npair) %>%
    select({{gene1}},{{gene2}},{{corr}})
  
  print(paste0("Top",npair,"Most correlated Pairs"))
  print(CorrPair)
  
  AntiCorrPair<-df%>%
    dplyr:: slice_min({{corr}},n=npair) %>%
    select({{gene1}},{{gene2}},{{corr}})
  
  print(paste0("Top",npair,"Most Antocorrelated Pairs"))
  print( AntiCorrPair)
  
}

## Sample function Call  
GetTopCorr(df=corr_fin,
           gene1=g_1,
           gene2=g_2,
           corr=pearson_corr,
           npair=10)