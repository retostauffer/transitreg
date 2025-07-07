


.PHONY: pdf
pdf:
	R CMD Sweave stauffer.Rnw
	pdflatex stauffer.tex
	pdflatex stauffer.tex
	pdflatex stauffer.tex
	#make clean

slides:
	R CMD Sweave slides.Rnw
	pdflatex slides.tex
	pdflatex slides.tex
	make clean

figure1:
	Rscript make_figure1.R

clean:
	-rm *.tex
	-rm *.nav
	-rm *.out
	-rm *.snm
	-rm *.toc
	-rm *.vrb
