FOR /f "usebackq" %%x IN (`Rscript -e "Rcpp:::CxxFlags()"`) DO SET PKG_CXXFLAGS=%%x
R CMD SHLIB biSBMWin_R.cpp