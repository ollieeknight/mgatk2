.PHONY: env setup format lint test check build check-all fixtures run tenx call run-multi integration all benchmark benchmark-hdf5 clean

env: setup

setup:
	python3 -m venv .venv
	.venv/bin/pip install --upgrade pip
	.venv/bin/pip install -e '.[dev]'

format:
	.venv/bin/ruff format src/ tests/
	.venv/bin/ruff check --fix src/ tests/

lint:
	.venv/bin/ruff check src/ tests/
	.venv/bin/ruff format --check src/ tests/

test:
	.venv/bin/pytest -q

check:
	$(MAKE) lint
	$(MAKE) test

build:
	.venv/bin/python -m build --no-isolation

check-all: check build
	.venv/bin/mgatk2 --help
	.venv/bin/mgatk2 run --help
	.venv/bin/mgatk2 tenx --help
	.venv/bin/mgatk2 call --help
	.venv/bin/mgatk2 paired --help
	.venv/bin/mgatk2 wes --help

fixture:
	.venv/bin/python tests/create_integration_fixture.py .test-work

fixtures:
	.venv/bin/python tests/create_integration_fixture.py tests/fixtures/10x_atac

run:
	.venv/bin/mgatk2 run --help
	.venv/bin/mgatk2 run -i tests/fixtures/10x_atac/outs -o .test-work/run-hdf5 -f hdf5 -t 1
	.venv/bin/mgatk2 run -i tests/fixtures/10x_atac/outs -o .test-work/run-txt -f txt -t 1

tenx:
	.venv/bin/mgatk2 tenx --help
	.venv/bin/mgatk2 tenx -i tests/fixtures/10x_atac/outs -o .test-work/tenx -t 1

call:
	.venv/bin/mgatk2 call --help
	.venv/bin/mgatk2 call -i tests/fixtures/10x_atac/outs -o .test-work/call -t 1

run-multi:
	.venv/bin/mgatk2 run -i tests/fixtures/10x_multi/outs -o .test-work/run-multi -f hdf5 -t 1

integration: run tenx call run-multi

all: check-all integration

benchmark: fixture
	cd tests && ../.venv/bin/python benchmark.py

benchmark-hdf5: fixture
	cd tests && ../.venv/bin/python benchmark.py --format hdf5

clean:
	rm -rf build dist *.egg-info
	rm -rf .ruff_cache .venv .mypy_cache .flake8
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type d -name *.egg-info -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	find . -type f -name ".DS_Store" -delete
