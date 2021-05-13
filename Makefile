all: initc docs test check

initc:
	R --slave -e "Rcpp::compileAttributes()"
	R --slave -e "tools::package_native_routine_registration_skeleton('.', 'src/init.c', character_only = FALSE)"

docs: man readme vigns

data:
	Rscript --slave inst/extdata/simulate_data.R

readme:
	R --slave -e "rmarkdown::render('README.Rmd')"

man:
	R --slave -e "devtools::document()"

vigns:
	rm -f doc/*.html
	rm -f vignettes/*.html
	rm -f inst/doc/*.html
	R --slave -e "options(rmarkdown.html_vignette.check_title = FALSE);devtools::build_vignettes(install = FALSE)"
	rm -f vignettes/*.utf8.md
	rm -f vignettes/*.md
	cp -Rf doc inst/

purl_vigns:
	R --slave -e "lapply(dir('vignettes', '^.*\\\\.Rmd$$'), function(x) knitr::purl(file.path('vignettes', x), gsub('.Rmd', '.R', x, fixed = TRUE)))"

quicksite:
	R --slave -e "options(rmarkdown.html_vignette.check_title = FALSE);pkgdown::build_site(run_dont_run = TRUE, lazy = TRUE)"
	cp -Rf doc inst/

site:
	R --slave -e "pkgdown::clean_site()"
	R --slave -e "options(rmarkdown.html_vignette.check_title = FALSE);pkgdown::build_site(run_dont_run = TRUE, lazy = TRUE)"
	cp -Rf doc inst/

test:
	R --slave -e "devtools::test()" > test.log 2>&1
	rm -f tests/testthat/Rplots.pdf

quickcheck:
	echo "\n===== R CMD CHECK =====\n" > check.log 2>&1
	R --slave -e "devtools::check(build_args = '--no-build-vignettes', args = '--no-build-vignettes', run_dont_test = TRUE, vignettes = FALSE)" >> check.log 2>&1se
	cp -Rf doc inst/

check:
	echo "\n===== R CMD CHECK =====\n" > check.log 2>&1
	R --slave -e "devtools::check(build_args = '--no-build-vignettes', args = '--no-build-vignettes', run_dont_test = TRUE, vignettes = FALSE)" >> check.log 2>&1
	cp -Rf doc inst/

install:
	R --slave -e "devtools::install_local('../surveyvoi', force = TRUE, upgrade = 'never')"

examples:
	rm -f examples.log
	R --slave -e "devtools::run_examples(test = TRUE);warnings()"  >> examples.log
	rm -f Rplots.pdf

wbcheck:
	R --slave -e "devtools::check_win_devel()"
	cp -Rf doc inst/

solarischeck:
	R --slave -e "rhub::check(platform = 'solaris-x86-patched', email = 'jeffrey.hanson@uqconnect.edu.au', show_status = FALSE)"

sancheck:
	R --slave -e "rhub::check(platform = 'linux-x86_64-rocker-gcc-san', email = 'jeffrey.hanson@uqconnect.edu.au', show_status = FALSE)"

vgcheck:
	R --slave -e "rhub::check(platform = 'ubuntu-gcc-devel', valgrind = TRUE, email = 'jeffrey.hanson@uqconnect.edu.au', show_status = FALSE)"

.PHONY: initc docs data site test check checkwb build install man readme vigns site quicksite benchmark examples solarischeck wbcheck sancheck vgcheck
