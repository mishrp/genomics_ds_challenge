#=======================#
# 
#  Genomics Exercise
#  Part 3
#=======================#

## Part 3: Find top 10 correlated gene pairs & other highly correlted or anticorrelated gene pairs.
##         Higly correlated = correlation coefficient above 7 and p-value less than 0.05. 
##        Add multiple hypothesis correction to the p-value



#libraries

library(tidyverse) # data manipulation and plottin
library(data.table) # to save data frame as data.table


options(scipen=50)


# 3a. highly correlated genes
eval_gene<-unique(as.vector(as.matrix(corr10[,1:2])))

high_corr_all<-map_dfr(1:length(eval_gene),
                       function(x) data.frame(id=eval_gene[x],
                                              filter(corr_fin, (g_1==eval_gene[x]|g_2==eval_gene[x]) & abs(pearson_corr)>0.7 & pvalue<0.05)))

# 3b. emove duplicate pair [i.e., alredy identified in corr 10]
high_corr_retain<-high_corr_all%>% 
  anti_join(corr10,by=c("g_1","g_2")) %>% 
  mutate(g_3=ifelse(id==g_1,g_2,g_1)) %>% 
  select(id,g_3,pearson_corr,pvalue)

head(high_corr_retain)


# 3c.Create data frame of gene 1, gene 2 and gene 3 (that is highly correlation with genes 1 or 2)

merg_g1<-select(corr10,g_1,g_2)%>% 
  inner_join (high_corr_retain, by=c("g_1"="id")) 


high_corr_set<-select(corr10,g_1,g_2)%>%
  inner_join(high_corr_retain, by=c("g_2"="id")) %>% 
  bind_rows(merg_g1) %>%  
  arrange(g_1) %>% 
  mutate(gene_pair=paste0(g_1,",",g_2)) %>% 
  unique(.)


# 3d. get the number of tests that were performed in deriving  corr_set 
(num_test=length(eval_gene)* (nrow(nrv_cln)-length(eval_gene)))

#3e. Bonferroni correction
(bonf_pvalue<-0.05/num_test)

#3f. bonferroni-corrected  p-values

high_corr_fin<-high_corr_set %>% 
  mutate(bonf_corr_pvalue=bonf_pvalue) %>% 
  filter(pvalue<bonf_corr_pvalue) %>% 
  select(g_1,g_2,g_3,pearson_corr,pvalue,bonf_corr_pvalue)


