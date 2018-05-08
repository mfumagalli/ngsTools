
export HTSSRC=../htslib/

EXT = htslib angsd
TOOLS = ngsSim ngsPopGen ngsUtils ngsDist ngsLD ngsF ngsF-HMM

all: $(EXT) $(TOOLS)

.PHONY: $(EXT) $(TOOLS)
$(EXT) $(TOOLS):
	@$(MAKE) -C $@

test:
	@for i in $(TOOLS); do $(MAKE) -C $$i test; done

clean:
	@for i in $(EXT) $(TOOLS); do $(MAKE) -C $$i clean; done
