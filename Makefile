PACKAGE = scbursts
VERSION = 1.2

all: deps docs build check

update-deps:
	Rscript -e 'install.packages("devtools",  repos="http://cran.rstudio.com")'
	Rscript -e 'install.packages("roxygen2",  repos="http://cran.rstudio.com")'
	Rscript -e 'install.packages("rmarkdown", repos="http://cran.rstudio.com")'
	Rscript -e 'install.packages("knitr",     repos="http://cran.rstudio.com")'
	Rscript -e 'install.packages("readxl",     repos="http://cran.rstudio.com")'
	Rscript -e 'install.packages("tinytex",   repos="http://cran.rstudio.com")'

deps:
	Rscript -e 'if (!require("devtools"))  install.packages("devtools",  repos="http://cran.rstudio.com")'
	Rscript -e 'if (!require("roxygen2"))  install.packages("roxygen2",  repos="http://cran.rstudio.com")'
	Rscript -e 'if (!require("rmarkdown")) install.packages("rmarkdown", repos="http://cran.rstudio.com")'
	Rscript -e 'if (!require("knitr"))     install.packages("knitr",     repos="http://cran.rstudio.com")'	
	Rscript -e 'if (!require("readxl"))     install.packages("readxl",     repos="http://cran.rstudio.com")'
	Rscript -e 'if (!require("tinytex"))   install.packages("tinytex",   repos="http://cran.rstudio.com")'

docs:
	R -e 'devtools::document()'

build: deps NAMESPACE
	R CMD build .

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



