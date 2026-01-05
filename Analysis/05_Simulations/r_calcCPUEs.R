#--get observed CPUEs by haul
  #----want to use tcsamSurveyData functions for this: they scale by area swept
  require(ggplot2);
  require(rlang);
  require(tables);
  require(wtsSizeComps);
  
  dirThs = dirname(rstudioapi::getActiveDocumentContext()$path)
  dfr = wtsUtilities::getObj(file.path(dirThs,"../04_HaulLevelCatchability/rda_Step1b_dfrPropsAllSBSYears.RData"));
  
  #--calc CPUE by year, sex, size, and haul
  dfrCPUE = dfr |> dplyr::select(YEAR,HAULJOIN,SEX,SIZE,numNMFS,aswNMFS,expfNMFS,numBSFRF,aswBSFRF,expfBSFRF) |> 
              dplyr::mutate(cpueNMFS=numNMFS*expfNMFS,cpueBSFRF=numBSFRF*expfBSFRF) |> 
              dplyr::rename_with(tolower,.cols=c(HAULJOIN,YEAR,SEX,SIZE));
  
  ggplot(dfrCPUE,aes(x=size)) + 
    geom_point(aes(y=cpueNMFS),colour="green",alpha=0.3) + 
    geom_point(aes(y=cpueBSFRF,x=size+2),colour="blue",alpha=0.3) + 
    scale_y_log10() + 
    labs(x="size (mm CW)",y="CPUE (num/sq. nmi.)") + 
    wtsPlots::getStdTheme()
    
  ggplot(dfrCPUE |> dplyr::filter(numBSFRF>4),aes(x=size)) + 
    geom_smooth(aes(y=cpueNMFS),colour="green",alpha=0.3) + 
    geom_point(aes(y=cpueNMFS),colour="green",alpha=0.3) + 
    geom_smooth(aes(y=cpueBSFRF,x=size+2),colour="blue",alpha=0.3) + 
    geom_point(aes(y=cpueBSFRF,x=size+2),colour="blue",alpha=0.3) + 
    scale_y_log10() + 
    labs(x="size (mm CW)",y="CPUE (num/sq. nmi.)") + 
    wtsPlots::getStdTheme()
    
  #--calc mean and variance of CPUE by year,sex,size
  dfrMnCPUE = dfrCPUE |> dplyr::group_by(year,sex,size) |> 
                dplyr::summarize(cpueNMFS=mean(cpueNMFS),
                                 varNMFS=var(cpueNMFS),
                                 cpueBSFRF=mean(cpueBSFRF),
                                 varBSFRF=var(cpueBSFRF))
  ggplot(dfrMnCPUE,aes(x=size)) + 
    geom_smooth(aes(y=cpueNMFS),colour="green",alpha=0.3) + 
    geom_point(aes(y=cpueNMFS),colour="green",alpha=0.3) + 
    geom_smooth(aes(y=cpueBSFRF,x=size+2),colour="blue",alpha=0.3) + 
    geom_point(aes(y=cpueBSFRF,x=size+2),colour="blue",alpha=0.3) + 
    scale_y_log10() + 
    labs(x="size (mm CW)",y="mean CPUE (num/sq. nmi.)") + 
    wtsPlots::getStdTheme()
