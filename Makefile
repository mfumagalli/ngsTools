
TOOLS = ngsSim angsd ngsF ngsPopGen ngsUtils

all: $(TOOLS)

.PHONY: $(TOOLS)
$(TOOLS):
	$(MAKE) -C $@;

test:
	@for i in $(TOOLS); do $(MAKE) -C $$i test; done

clean:
	@for i in $(TOOLS); do $(MAKE) -C $$i clean; done
