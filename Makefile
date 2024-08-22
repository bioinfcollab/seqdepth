TMP_DATADIR := $(shell mktemp -d)
NF_TMPDIR := $(shell mktemp -d)
SINGULARITY_TMPDIR := $(shell mktemp -d)
container := ${CONTAINER}
containersdir := 'containers'
seqdepth_scripts := ${seqdepth_scripts}
ifndef seqdepth_scripts
override seqdepth_scripts = ./scripts
endif
ifndef container
override container=$(shell readlink -f ./containers/seqdepth.sif)
endif
bold := $(shell tput bold)
sgr0 := $(shell tput sgr0)

REQBIN    := singularity nextflow
REQBIN_download   := wget

.PHONY: default
default: info;

info:
	@echo Seqdepth pipeline requires nextflow and singularity
	@echo the following targets are available:
	@echo - $(bold)build_container:$(sgr0)
	@echo "\t"+ will build containers/seqdepth.sif with required tools and R
	@echo "\t"+ see def/seqdepth.def_tmpl for details
	@echo "\t"+ requires configured singularity \(adjusted /etc/sub{uid,gid}\) or sudo privileges
	@echo - $(bold)download_container:$(sgr0)
	@echo "\t"+ download seqdepth container avoiding the need to build it
	@echo - $(bold)test:$(sgr0)
	@echo "\t"+ run a test by downloading small dataset and running scripts on it

testscripts:
	cd scripts && prove
gendata:
	$(info $(bold)Downloading test data$(sgr0))
	cp scripts/gentestfiles.sh $(TMP_DATADIR)
	cd $(TMP_DATADIR) && singularity exec -B $(PWD) $(container) ./gentestfiles.sh

build_container:
	$(foreach bin,$(REQBIN),\
		$(if $(shell command -v $(bin) 2> /dev/null),,$(error `$(bin)` is missing, aborting...)))
ifndef seqdepth_scripts
	$(info example: $(bold)make build_containers seqdepth_scripts=/path/to/seqdepth$(sgr0))
	$(error $(bold)seqdepth_scripts$(sgr0) env var is not set)
endif
	perl -pe 's{<seqdepth_scripts>}{$(seqdepth_scripts)}' def/seqdepth.def_tmpl > $(SINGULARITY_TMPDIR)/seqdepth.def
	singularity build --fakeroot $(containersdir)/seqdepth.sif $(SINGULARITY_TMPDIR)/seqdepth.def
	@test -d $(SINGULARITY_TMPDIR) && rm -rf $(SINGULARITY_TMPDIR)

download_container:
	$(foreach bin,$(REQBIN_download),\
		$(if $(shell command -v $(bin) 2> /dev/null),,$(error `$(bin)` is missing, aborting...)))
	wget -O containers/seqdepth.sif https://sharing.biotec.tu-dresden.de/index.php/s/bxPbyoEctgsfMDc/download

runtest:
ifndef container
	$(error $(bold)CONTAINER$(sgr0) env var is not set)
endif
	$(info $(bold)Running nextflow pipeline$(sgr0))
	nextflow run main.nf -profile local -w $(NF_TMPDIR) --datadir $(TMP_DATADIR) --outdir $(NF_TMPDIR)/result --container $(container) --temp $(NF_TMPDIR) --hpcproject test --containerOptions --task hypergeometric --split 10

.PHONY: clean
clean:
	$(info $(bold)Cleanup temporary data$(sgr0))
	@echo "Cleaning ${TMP_DATADIR} and ${NF_TMPDIR}"
	@test -d $(TMP_DATADIR) && rm -rf $(TMP_DATADIR)
	@test -d $(NF_TMPDIR) && rm -rf $(NF_TMPDIR)

.PHONY: test
test:
ifeq ("$(wildcard $(container))","")
	$(info example: $(bold)make test CONTAINER=/full/path/to/singularity_container$(sgr0))
	$(error $(bold)CONTAINER$(sgr0) env var is not set)
endif
	$(MAKE) gendata runtest clean
