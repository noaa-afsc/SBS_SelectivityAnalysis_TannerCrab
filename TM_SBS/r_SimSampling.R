
r_a = 6; #--ratio of NMFS area swept to BSFRF area swept
i = 0;
lstRes = list();
for (l_B in seq(2,100,2)) {
  for (S_N in seq(0.1,1.0,0.1)){
    nB = rpois(100,l_B);
    nN = rpois(100,l_B*r_a*S_N)
    lstRes[[i<-i+1]] = tibble::tibble(l_B=l_B,S_N=S_N,nB=nB,nN=nN,r=nN/(r_a*nB));
  }
}
dfrRes = dplyr::bind_rows(lstRes); 
rm(lstRes);
dfr = dfrRes |> dplyr::filter(l_B==6)
ggplot(dfr,aes(x=S_N,y=r)) + 
  geom_violin(aes(group=S_N)) +
  geom_abline(slope=1,linetype=3) +
  geom_point(data=dfr |> dplyr::group_by(S_N) |> 
                         dplyr::summarize(mean=mean(r[is.finite(r)],na.rm=TRUE)) |> 
                         dplyr::ungroup(),
             mapping=aes(x=S_N,y=mean),shape=21,colour="green",inherit.aes=FALSE);


ggplot(dfrZCs_RawByStn |> dplyr::filter(gear=="BSFRF",sex=="MALE"),
       aes(x=size,y=caught)) + 
  geom_point(alpha=0.2) + 
  facet_wrap(~year,scales="free_y")

