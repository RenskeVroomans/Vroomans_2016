#
# Simple C++ Makefile by Anton Crombach, A.B.M.Crombach@bio.uu.nl
#

OBJPATH = obj/
BINPATH = bin/

# Targets
all: 
	@cd $(OBJPATH); \
	make all

ana:
	@cd $(OBJPATH); \
	make all

slib:
	@cd $(OBJPATH); \
	make slib

fluq:
	@cd $(OBJPATH); \
	make fluq

.PHONY: clean realclean distclean
clean:
	@cd $(OBJPATH); make clean

distclean: 
	@cd $(OBJPATH); make distclean; \
	cd ../$(BINPATH); rm -f *

realclean:
	@cd $(OBJPATH); make realclean

