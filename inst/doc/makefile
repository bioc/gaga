gagamanual: gagamanual.Rnw
	echo 'Sweave("gagamanual.Rnw")' | R --vanilla --quiet
	pdflatex gagamanual.tex
	bibtex gagamanual
	pdflatex gagamanual.tex
	pdflatex gagamanual.tex
	rm gagamanual.tex
	rm gagamanual.log
	rm gagamanual.aux
	rm gagamanual.out
	rm gagamanual.blg
	rm gagamanual.bbl
	rm gagamanual-fig*.pdf
	rm gagamanual-fig*.eps
