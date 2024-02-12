#--plot empirical selectivity results
dirPrj = rstudioapi::getActiveProject();
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);

source(file.path(dirThs,"r_PlotEmpiricalSelectivity.R"));
lst = wtsUtilities::getObj(file.path(dirThs,"rda_Step1_EmpiricalSelectivityFromData.RData"));

#--create plots
pdf(file.path(dirThs,"fig_EmpiricalSelectivityFromData.pdf"),width=8,height=6);

#--females
dfrZCsp<-dfrZCs[(dfrZCs$x=="female")&(dfrZCs$z<=125),];
dfrESsp<-dfrESs[(dfrESs$x=="female")&(dfrESs$z<=125),];
plotEmpiricalSelectivity(dfrZCsp,dfrESsp,
                          plotPoints=TRUE,
                          points=list(alpha=0.2,size=3,dodge=0),
                          plotLines=TRUE,
                          plotViolins=FALSE,
                          plotSmooths=TRUE,
                          smooths=list(method="gam",formula=y~s(x,bs="cs",k=5),knots=c(25,50,75,100,125)));
#--males
dfrZCsp<-dfrZCs[(dfrZCs$x=="male")&(dfrZCs$z<=185),];
dfrESsp<-dfrESs[(dfrESs$x=="male")&(dfrESs$z<=185),];
plotEmpiricalSelectivity(dfrZCsp,dfrESsp,
                          plotPoints=TRUE,
                          points=list(alpha=0.2,size=3,dodge=0),
                          plotLines=TRUE,
                          plotViolins=FALSE,
                          plotSmooths=TRUE,
                          smooths=list(method="gam",formula=y~s(x,bs="cs",k=7),knots=c(25,50,75,100,125,150,175)));
dev.off();

#--alernative version
require(ggplot2)
p1 = ggplot(data=dfrZCs,mapping=aes_string(x="z",y="n",colour="fleet")) +
       geom_line() + geom_point() +
       facet_grid(rows = y~x,scales = "free_y") +
       labs(x="size (mm CW)",y="crab sampled",colour="")
ggsave("fig_ComparisonNumSampled.pdf",p1,device=pdf,width=9,height=8);
p1 = ggplot(data=dfrZCs,mapping=aes_string(x="z",y="val",colour="fleet")) +
       geom_line() + geom_point() +
       facet_grid(rows = y~x,scales = "free_y") +
       labs(x="size (mm CW)",y="abundance (millions)",colour="")
ggsave("fig_ComparisonZCs.pdf",p1,device=pdf,width=9,height=8);
tmp = dfrESs;
tmp$totN = tmp$n_BSFRF+tmp$n_NMFS;
tmp$y = factor(tmp$y)
p1 = ggplot(data=tmp,mapping=aes_string(x="z",y="emp_sel",colour="y",size="totN")) +
       geom_point() + scale_size_area() + 
       geom_line(mapping=aes_string(x="z",y="emp_sel",colour="y"),inherit.aes=FALSE) +
       geom_smooth(mapping=aes_string(x="z",y="emp_sel"),inherit.aes = FALSE) +
       facet_grid(rows = x~.,scales = "free_y") +
       labs(x="size (mm CW)",y="empirical\ncatchability",colour="study\nyear",size="crab\nsampled")
ggsave("fig_ComparisonESs.pdf",p1,device=pdf,width=6.5,height=5);

#--create table of numbers of crab sampled, by year, sex, and gear
require(tables);
tmp  = dfrZCs |> dplyr::group_by(fleet,y,x) |>
                  dplyr::summarize(`number sampled`=wtsUtilities::Sum(n),
                            `total abundance`=wtsUtilities::Sum(val)) |>
                  dplyr::ungroup() |>
                  tidyr::pivot_longer(cols=c("number sampled","total abundance"),names_to="type");
tbl = tables::tabular(Heading("year")*Factor(y,name="year")~
                        Factor(type)*
                        Factor(x,name="sex")*(
                          Heading("group")*Factor(fleet)*value*wtsUtilities::Sum
                        ),data=tmp);
colLabels(tbl) = colLabels(tbl)[c(2,4,6)];
write.csv.tabular(tbl,file="SBS_Summary.csv");


