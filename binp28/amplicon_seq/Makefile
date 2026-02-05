PYTHON ?= python3
DATA_DIR := data/

# Apptainer images
IMG_DIR := bin/images/
FASTQC_IMG := fastqc.sif
TRIMMOMATIC_IMG := trimmomatic.sif
PANDASEQ_IMG := pandaseq.sif
VSEARCH_IMG := vsearch.sif

# Apptainer URLS
FASTQC_URL := oras://community.wave.seqera.io/library/fastqc:0.12.1--104d26ddd9519960
TRIMMOMATIC_URL := docker://staphb/trimmomatic:latest
VSEARCH_URL := oras://community.wave.seqera.io/library/vsearch:2.30.3--db74bd085ba128f3
PANDASEQ_URL := docker://dromero93/pandaseq:latest

# VSEARCH BINARIES
VSEARCH_BIN_URL := https://github.com/torognes/vsearch/releases/download/v2.30.4/vsearch-2.30.4-linux-x86_64.tar.gz

# Reference Datatbase for chimera checking
CHIMERA_DATABASE := https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.gold.bacteria.zip

# RDP Classifier Database
RDP_CLASSIFIER := https://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/rdp_classifier_2.14.zip

.DEFAULT_GOAL := all
SHELL := bash
.PHONY := all download clean install
.SUFFIXES:
.DELETE_ON_ERROR:



download: # Download the required files from the course servver
	@if [ ! -d data/ ]; then \
		echo "data folder does not exists.."; \
		echo "Create the data folder."; \
		mkdir data; \
	else \
		echo "data folder already exists.."; \
	fi
	@if [ ! -f data/AmpliconData.tar.gz ]; then \
		echo "Copying data from the resource folder.." ; \
		cp -r /resources/binp28/Data/AmpliconData.tar.gz data/; \
	fi
	@echo "Un-packing the contents of the file.."
	@cd data/ && tar -xvzf AmpliconData.tar.gz && cd AmpliconData/ && unzip fastq_files.zip
	@echo "Clearing files.."
	@cd data/AmpliconData/ && rm -rf fastq_files.zip XML*
	@echo "[download] ok"
 
# Install all the softwares required
apptainer-install: $(FASTQC_IMG) $(TRIMMOMATIC_IMG) $(VSEARCH_IMG) $(PANDASEQ_IMG) 
	@echo "[install] All singularity images pulled."


# Rules to pull the images from Sequera and docker
# $@ -> means that it checks for the target file, in this case it would be the fastqc.sif
$(FASTQC_IMG):

	@if [ ! -d "$(IMG_DIR)" ]; then \
		echo "Image directory does not exists. Creating one..";\
		mkdir -p "$(IMG_DIR)";\
	else \
		echo "Directory already exists.";\
	fi

	@if [ ! -f bin/images/$@ ]; then \
		echo "Pulling FASTQC image.."; \
		apptainer pull bin/images/$@ $(FASTQC_URL); \
	else \
		echo "FASTQC image already exists."; \
	fi

$(TRIMMOMATIC_IMG):
	@if [ ! -f bin/images/$@ ]; then \
	    echo "Pulling Trimmomatic image..."; \
	    apptainer pull bin/images/$@ $(TRIMMOMATIC_URL); \
	else \
	    echo "Trimmomatic image already exists."; \
	fi

$(PANDASEQ_IMG):
	@if [ ! -f bin/images/$@ ]; then \
	    echo "Pulling Pandaseq image..."; \
	    apptainer pull bin/images/$@ $(PANDASEQ_URL); \
	else \
	    echo "Pandaseq image already exists."; \
	fi

$(VSEARCH_IMG):
	@if [ ! -f bin/images/$@ ]; then \
	    echo "Pulling Vsearch image..."; \
	    apptainer pull bin/images/$@ $(VSEARCH_URL); \
	else \
	    echo "Vsearch image already exists."; \
	fi

install-vsearch-binaries: # Install vsearch according to the exercies
	@if [ ! -d "$(IMG_DIR)" ]; then \
	    mkdir -p "$(IMG_DIR)"; \
	fi
	@cd bin/ && if [ ! -d vsearch_github_bin ]; then \
	    echo "Downloading VSEARCH binaries..."; \
	    wget "$(VSEARCH_BIN_URL)"; \
	    echo "Extracting..."; \
	    tar -xvzf vsearch-2.30.4-linux-x86_64.tar.gz; \
	    mv vsearch-2.30.4-linux-x86_64 vsearch_github_bin; \
	    rm vsearch-2.30.4-linux-x86_64.tar.gz; \
	    echo "VSEARCH installed."; \
	else \
	    echo "VSEARCH binaries already exist."; \
	fi

download-chimera-database: # Install the chimera database
	@if [ ! -d "$(DATA_DIR)/chimera_database" ]; then \
		echo "Chimera database not present. Downloading.."; \
		mkdir -p "$(DATA_DIR)/chimera_database" && cd "$(DATA_DIR)/chimera_database"; \
		wget "$(CHIMERA_DATABASE)"; \
		echo "Extracting";\
		unzip -p silva.gold.bacteria.zip silva.gold.align | sed -e "s/[.-]//g" > gold.fasta ;\
		echo "Database extracted". ;\
	else \
		echo "Database already exists.";\
	fi

download-rdp-database: # Download the latest version of the RDP database
	@if [ ! -d "$(DATA_DIR)/rdp_database/rdp_classifier_2.14" ]; then \
		echo "RDP Classifier database does not exist.."; \
		mkdir -p "$(DATA_DIR)/rdp_database" && cd "$(DATA_DIR)/rdp_database";\
		wget "$(RDP_CLASSIFIER)";\
		echo "Extracting..";\
		unzip rdp_classifier_2.14.zip;\
		echo "Database extracted.";\
	else \
		echo "Database already exist.";\
	fi
