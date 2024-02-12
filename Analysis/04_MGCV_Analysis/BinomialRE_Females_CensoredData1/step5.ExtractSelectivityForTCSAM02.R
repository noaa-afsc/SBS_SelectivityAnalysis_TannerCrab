#--step 5: extract estimated annual selectivity for assessment model
require(dplyr);

dfrRbyY = wtsUtilities::getObj("./results_RData/dfrRbyY.RData");


mdfr = dfrRbyY %>% dplyr::select(y,z,mnR,seR) %>%
                   tidyr::pivot_longer(c("mnR","seR"),names_to="type",values_to="value") %>%
                   dplyr::arrange(type,y,z);

dfrMnRs = mdfr %>% tidyr::pivot_wider(id_cols=c("type","y"),names_from="z",values_from="value");

write.csv(dfrMnRs,file="./results_csvs/selectivityByYear.Mean.csv",row.names=FALSE);


