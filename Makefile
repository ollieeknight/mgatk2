.PHONY: all setup check clean

all: check

setup:
	python3 -m venv .venv
	.venv/bin/pip install --upgrade pip
	.venv/bin/pip install -e '.[dev]'

# Format, lint, unit tests, then one end-to-end run of every pipeline
# against the committed fixtures.
check:
	.venv/bin/ruff format src/ tests/
	.venv/bin/ruff check --fix src/ tests/
	.venv/bin/pytest -q
	.venv/bin/mgatk2 run -i tests/fixtures/10x_atac/outs -o .test-work/hdf5 -f hdf5 -t 2
	.venv/bin/mgatk2 run -i tests/fixtures/10x_atac/outs -o .test-work/txt -f txt -t 2
	.venv/bin/mgatk2 tenx -i tests/fixtures/10x_atac/outs -o .test-work/tenx -t 2
	.venv/bin/mgatk2 call -i tests/fixtures/10x_atac/outs -o .test-work/call -t 2
	.venv/bin/mgatk2 run -i tests/fixtures/10x_multi/outs -o .test-work/multi -f hdf5 -t 2

clean:
	rm -rf build dist .test-work *.egg-info src/*.egg-info
	rm -rf .venv .ruff_cache .pytest_cache
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	find . -type f -name ".DS_Store" -delete
