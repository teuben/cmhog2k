#

all:	hydro.ps  hydro_2.ps

hydro.ps:	hydro.dvi
	dvips hydro -o

hydro.dvi:	hydro.tex
	latex hydro


hydro_2.ps: hydro.ps
	psnup -2 hydro.ps  > hydro_2.ps
