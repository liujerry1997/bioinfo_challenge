# Makefile for testing, building and deploying

VERSION = v0.1.0
DOCKER_URL = hliu2023/tempus_challenge:$(VERSION)


all: format build test

build:
	docker buildx build . -t $(DOCKER_URL)

deploy:
	docker push $(DOCKER_URL)

format:
	black --include '.py$$' -l 100 bin tests

# Make sure your local environment has Nextflow and needed python dependencies installed
test:
	nextflow run main.nf -profile test

test_unit: format build
	docker run --rm --entrypoint pytest -v $(PWD):$(PWD) -w $(PWD) $(DOCKER_URL) -s tests 

# Use this make command to run the pipeline with specified profile
run: format build
	nextflow run main.nf -profile XXX

shell:
	docker run -it  -v $$PWD:$$PWD -w $$PWD $(DOCKER_URL) bash

