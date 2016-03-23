
EXT = htslib angsd
TOOLS = ngsSim ngsPopGen ngsUtils ngsDist ngsF


all: $(EXT) $(TOOLS)

.PHONY: $(EXT) $(TOOLS)
$(EXT) $(TOOLS):
	@$(MAKE) -C $@ HTSSRC=../htslib/;

test:
	@for i in $(TOOLS); do $(MAKE) -C $$i test; done

clean:
	@for i in $(EXT) $(TOOLS); do $(MAKE) -C $$i clean; done
