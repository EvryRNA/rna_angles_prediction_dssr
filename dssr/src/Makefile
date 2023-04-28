CFLAGS		= -ansi -pedantic -W -Wall -Wextra -Wunused -Wshadow -O3
CC		= gcc
LIBS		= -lm

OBJS		= \
		ana_fncs.o  \
		app_fncs.o  \
		cmn_fncs.o  \
		fncs_slre.o  \
		nrutil.o  \
		reb_fncs.o 

HEADERS		= \
		x3dna.h  \
		x3dna_fncs.h 

PROGRAMS	=  analyze  

all:		$(PROGRAMS)

# [02]: analyze
analyze:	$(OBJS) analyze.o
			$(CC) $(OBJS) analyze.o -o analyze $(LIBS)
			mv analyze ../bin
analyze.o:	$(HEADERS)

# Dependency of common objs on the headers
ana_fncs.o:	$(HEADERS)
app_fncs.o:	$(HEADERS)
cmn_fncs.o:	$(HEADERS)
fncs_slre.o:	$(HEADERS)
nrutil.o:	$(HEADERS)
reb_fncs.o:	$(HEADERS)

etags:
		etags *.[ch]

.PHONY : clean
clean:
		@rm -f *.o
