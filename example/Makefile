# Please enter here the locations for PipLib include and libraries if they
# aren't the default values (/usr/lib and /usr/include).
PIPLIB_INC = $(HOME)/progs/linux/include
PIPLIB_LIB = $(HOME)/progs/linux/lib

CC = gcc
LDLIBS= -lpiplib32
CFLAGS= -DLINEAR_VALUE_IS_INT -I $(PIPLIB_INC) -L $(PIPLIB_LIB)

ROOT = ..

PIPTEST = \
	big \
	brisebarre \
	cg1 \
	esced \
	ex \
	ex2 \
	fimmel \
	max \
	negative \
	small \
	sor1d \
	square \
	sven

BITS = 64

example: example.c 
	@echo "             /*-----------------------------------------------*"
	@echo "              *                Making example                 *"
	@echo "              *-----------------------------------------------*/"
	$(CC) example.c -o example $(CFLAGS) $(LDLIBS)

pip:
	$(MAKE) total -C ../

clean:
	@echo "             /*-----------------------------------------------*"
	@echo "              *               Cleaning example                *"
	@echo "              *-----------------------------------------------*/"
	-rm -f example example.exe core

test:
	failedtest=0; \
	for x in $(PIPTEST) ; do \
		echo "Verify file $$x" ; \
		$(ROOT)/obj$(BITS)_$(TARGET)/example$(BITS)$(EXEC_SUFFIX) < $$x.pip > xyz ;\
		diff -w xyz $$x.ll ; \
		result=$$?; \
		if [ "$$result" -eq "1" ]; then \
			echo "Error: $$x is not the same"; \
			failedtest=`expr $$failedtest + 1`; \
		else \
			echo "$$x passed"; \
	        fi; \
	done ; \
	rm xyz ; \
	if [ $$failedtest != 0 ]; then \
		echo "$$failedtest tests failed"; \
	else \
	        echo "PIP works correctly :-)"; \
        fi

# Include the shared compilation rules
include $(ROOT)/Makefile
