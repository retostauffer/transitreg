


.PHONY: pdf
pdf:
	R CMD Sweave stauffer.Rnw
	pdflatex stauffer.tex
	pdflatex stauffer.tex
	pdflatex stauffer.tex
	make clean


figure1:
	Rscript make_figure1.R

clean:
	-rm stauffer.aux
	-rm stauffer.log
	-rm stauffer.tex



