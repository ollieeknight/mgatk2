.PHONY: all setup check e2e clean

BIN ?= .venv/bin/
FIXTURES = tests/fixtures

all: check e2e

setup:
	python3 -m venv .venv
	.venv/bin/pip install --upgrade pip
	.venv/bin/pip install -e '.[dev]'

check:
	$(BIN)ruff format src/ tests/
	$(BIN)ruff check --fix src/ tests/
	$(BIN)vulture src tests
	$(BIN)pytest -q

# One end-to-end run of every single-cell pipeline against the committed fixtures.
e2e:
	$(BIN)mgatk2 run -i $(FIXTURES)/10x_atac/outs -o .test-work/hdf5 -f hdf5 -t 2
	$(BIN)mgatk2 run -i $(FIXTURES)/10x_atac/outs -o .test-work/txt -f txt -t 2
	$(BIN)mgatk2 tenx -i $(FIXTURES)/10x_atac/outs -o .test-work/tenx -t 2
	$(BIN)mgatk2 call -i $(FIXTURES)/10x_atac/outs -o .test-work/call -t 2
	$(BIN)mgatk2 run -i $(FIXTURES)/10x_multi/outs -o .test-work/multi -f hdf5 -t 2

clean:
	rm -rf build dist .test-work src/*.egg-info .venv .ruff_cache .pytest_cache
	find . -type d -name __pycache__ -exec rm -rf {} +
