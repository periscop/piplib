# Please enter here the locations for PipLib include and libraries if they
# aren't the default values (/usr/lib and /usr/include).
PIPLIB_INC = /usr/local/include #$(HOME)/progs/linux/include
PIPLIB_LIB = /usr/local/lib     #$(HOME)/progs/linux/lib

CC = gcc
LDLIBS= -lpiplib32
CFLAGS= -I $(PIPLIB_INC) -L $(PIPLIB_LIB)


example: example.c
	@echo "             /*-----------------------------------------------*"
	@echo "              *                Making example                 *"
	@echo "              *-----------------------------------------------*/"
	$(CC) example.c -o example $(CFLAGS) $(LDLIBS)

clean:
	@echo "             /*-----------------------------------------------*"
	@echo "              *               Cleaning example                *"
	@echo "              *-----------------------------------------------*/"
	-rm -f example example.exe core
