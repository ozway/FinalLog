MONTH=february15

default:
	pdflatex $(MONTH).tex

clean:
	@rm -fv *.log
	@rm -fv *.aux
	@rm -fv *.toc

cleanpdf: clean
	rm -fv *.pdf


tar: 
	tar czvf ../$(MONTH)Bio.tar.gz ../$(MONTH)/*.tex ../$(MONTH)/*.pdf ../$(MONTH)/data

backup: tar
