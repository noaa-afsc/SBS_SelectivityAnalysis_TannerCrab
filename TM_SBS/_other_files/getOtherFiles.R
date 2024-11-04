#--copy "other files" for processing qmd files from original locations
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
bib = path.expand("~/Work/Projects/Bibliography/AllRefs.bib");
csl = system.file("files/CJFAS.csl",package="wtsQMD");
ltx = system.file("files/ltx_ExtraLatexIncludes.tex",package="wtsQMD");
fls = c(bib,csl,ltx);
for (f in fls) file.copy(f,file.path(dirThs,basename(f)),overwrite=TRUE);

