ENV_NAME := logan
ENV_FILE := environment.yml
CONDA    := $(shell command -v mamba 2>/dev/null || command -v conda 2>/dev/null)

KO_SIG_URL  := https://zenodo.org/records/10045253/files/KOs_sketched_scaled_1000.sig.zip
KO_SIG_FILE := KOs_sketched_scaled_1000.sig.zip

METADATA_S3_URI   := s3://sra-pub-metadata-us-east-1/sra/metadata/
METADATA_DIR      ?= aws_sra_metadata
METADATA_DB       ?= sra_metadata.db
METADATA_BUILD_SQL := database/build_metadata_db.sql
ACCESSION_SQL      := database/wgs_metagenome.sql

.PHONY: install env create-db get-accessions data clean-env run-test

install: env data

env:
ifeq ($(CONDA),)
	$(error Neither mamba nor conda was found on PATH. Install Miniconda/Miniforge first)
endif
	@if $(CONDA) env list | awk '{print $$1}' | grep -qx "$(ENV_NAME)"; then \
		echo "Updating existing '$(ENV_NAME)' conda environment..."; \
		$(CONDA) env update -n $(ENV_NAME) -f $(ENV_FILE) --prune; \
	else \
		echo "Creating '$(ENV_NAME)' conda environment..."; \
		$(CONDA) env create -n $(ENV_NAME) -f $(ENV_FILE); \
	fi
	@echo "Done. Activate it with: conda activate $(ENV_NAME)"

create-db:
	# Synchronize the public SRA metadata snapshot, then build/update the DuckDB database.
	aws s3 sync --no-sign-request $(METADATA_S3_URI) "$(METADATA_DIR)"
	cd "$(METADATA_DIR)" && duckdb "$(abspath $(METADATA_DB))" < "$(abspath $(METADATA_BUILD_SQL))"

# One command performs metadata synchronization/database creation first, then writes
# wgs_metagenome_accessions.txt using the version-controlled selection query.
get-accessions: create-db
	duckdb "$(METADATA_DB)" < "$(ACCESSION_SQL)"

data: $(KO_SIG_FILE)

$(KO_SIG_FILE):
	wget -O $(KO_SIG_FILE) $(KO_SIG_URL)

clean-env:
ifeq ($(CONDA),)
	$(error Neither mamba nor conda was found on PATH)
endif
	$(CONDA) env remove -n $(ENV_NAME)

# Runs the Nextflow project against its small test manifest, inside the
# conda env that env: creates -- this is what makes `nextflow` itself reproducible.
run-test:
	$(CONDA) run -n $(ENV_NAME) nextflow run nextflow -profile test -resume
