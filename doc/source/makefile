NAME_SOURCE=pip
NAME=pip

all:
	latex $(NAME_SOURCE).tex
	latex $(NAME_SOURCE).tex
	dvips -o $(NAME).ps $(NAME_SOURCE).dvi
clean:
	rm -f $(NAME).aux $(NAME).log $(NAME).dvi

html:
	latex2html -t PipLib -white -address cedric.bastoul@inria.fr \
	-show_section_numbers \
	-antialias -split 4 -info 0 -html_version 3.2,math pip.tex

cleanall: clean
	rm -f $(NAME).ps

tgz: cleanall
	(cd .. && \
	tar cvf piplib-doc.tar piplib-doc && \
	gzip piplib-doc.tar)

