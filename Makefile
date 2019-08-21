# Not BSD Friendly. Sorry. Pull-Requests welcome.

PACKAGE = scbursts
VERSION = 1.7

LIBS := -lgmp -lmpfr -llikelihood

RCPPINCL := 		$(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RCPPLIBS := 		$(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)

c_sources := 		$(wildcard *.c)
c_sharedlibs := 	$(patsubst %.c,%.o,$(c_sources))

cpp_sources := 		$(wildcard *.cpp)
cpp_sharedlibs := 	$(patsubst %.cpp,%.o,$(cpp_sources))

%.o : %.c
	R CMD SHLIB $<

%.o : %.cpp
	PKG_CPPFLAGS="$(RCPPFLAGS) $(RCPPINCL)" PKG_LIBS="$(LIBS) $(RLDFLAGS) $(RCPPLIBS)" R CMD SHLIB $<

all: $(c_sharedlibs) $(cpp_sharedlibs) deps docs build check

all-fast: $(c_sharedlibs) $(cpp_sharedlibs) deps docs build fastcheck

update-deps:
	Rscript -e 'install.packages("devtools",  repos="http://cran.rstudio.com")'
	Rscript -e 'install.packages("roxygen2",  repos="http://cran.rstudio.com")'
	Rscript -e 'install.packages("rmarkdown", repos="http://cran.rstudio.com")'
	Rscript -e 'install.packages("knitr",     repos="http://cran.rstudio.com")'
	Rscript -e 'install.packages("readxl",    repos="http://cran.rstudio.com")'
	Rscript -e 'install.packages("tibble",    repos="http://cran.rstudio.com")'
	Rscript -e 'install.packages("tinytex",   repos="http://cran.rstudio.com")'
	Rscript -e 'install.packages("Rcpp",      repos="http://cran.rstudio.com")'
	Rscript -e 'install.packages("RcppEigen", repos="http://cran.rstudio.com")'

deps:
	Rscript -e 'if (!require("devtools"))  install.packages("devtools",  repos="http://cran.rstudio.com")'
	Rscript -e 'if (!require("roxygen2"))  install.packages("roxygen2",  repos="http://cran.rstudio.com")'
	Rscript -e 'if (!require("rmarkdown")) install.packages("rmarkdown", repos="http://cran.rstudio.com")'
	Rscript -e 'if (!require("knitr"))     install.packages("knitr",     repos="http://cran.rstudio.com")'	
	Rscript -e 'if (!require("readxl"))    install.packages("readxl",    repos="http://cran.rstudio.com")'
	Rscript -e 'if (!require("tibble"))    install.packages("tibble",    repos="http://cran.rstudio.com")'
	Rscript -e 'if (!require("tinytex"))   install.packages("tinytex",   repos="http://cran.rstudio.com")'
	Rscript -e 'if (!require("Rcpp"))      install.packages("Rcpp",      repos="http://cran.rstudio.com")'
	Rscript -e 'if (!require("RcppEigen")) install.packages("RcppEigen", repos="http://cran.rstudio.com")'

docs:
	R -e 'devtools::document()'

build: deps NAMESPACE $(c_sharedlibs) $(cpp_sharedlibs) 
	R -e 'Rcpp::compileAttributes()'
	Rscript -e 'library(roxygen2); roxygenize()'
	PKG_CPPFLAGS="$(RCPPFLAGS) $(RCPPINCL)" PKG_LIBS="$(LIBS) $(RLDFLAGS) $(RCPPLIBS)" R CMD build .

NAMESPACE:
	$(MAKE) docs

$(PACKAGE)_$(VERSION).tar.gz:
	$(MAKE) build

install: $(PACKAGE)_$(VERSION).tar.gz
	R CMD INSTALL $(PACKAGE)_$(VERSION).tar.gz

check: $(PACKAGE)_$(VERSION).tar.gz
	R CMD check $(PACKAGE)_$(VERSION).tar.gz --as-cran

fastcheck: build
	R CMD check $(PACKAGE)_$(VERSION).tar.gz

clean:
	$(RM) -r $(PACKAGE).Rcheck/
	$(RM) $(wildcard src/*.o)
	$(RM) src/scbursts.so
	$(RM) $(PACKAGE)_$(VERSION).tar.gz

$(PACKAGE).Rcheck:
	$(MAKE) check

export: $(PACKAGE)_$(VERSION).tar.gz $(PACKAGE).Rcheck
	@echo Copying tarball and manuals to ../build/
	@mkdir -p ../build/
	@cp scbursts.Rcheck/scbursts-manual.pdf ../build/
	@cp scbursts.Rcheck/scbursts/doc/scbursts.pdf ../build/
	@cp $(PACKAGE)_$(VERSION).tar.gz ../build
	$(MAKE) clean



