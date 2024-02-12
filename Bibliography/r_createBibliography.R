#--assemble bib style bibliography
pth = "~/Work/Projects/Bibliography";#--path to main bibliography folder

#--review available files
cat(paste("'",list.files(path=pth,pattern=glob2rx("*.bib")),"'",sep=''),sep=",\n");

#--select files to concatenate
files = c('bib_ADFG.bib',
          'bib_ADMB.bib',
          'bib_CPT.bib',
          'bib_GeneralRefs.bib',
          'bib_GPT.bib',
          'bib_NMFS-Surveys.bib',
          'bib_NPFMC.bib',
          'bib_PIBKC.bib',
          'bib_R-packages.bib',
          'bib_SBS_Studies.bib',
          'bib_SSC.bib',
          'bib_TannerCrab.bib',
          'bib_TannerCrabSAFEs.bib'
          )
bib = "%--Tanner crab SBS Selectivity Studies bibliography file";
for (f in files){
  fn = file.path(pth,f);
  txt = readLines(con=fn);
  bib = c(bib,txt);
}

dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
cat(bib,sep="\n",file=file.path(dirThs,"bib_SBS-TannerCrab-doc.bib"));

