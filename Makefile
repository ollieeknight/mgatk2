.PHONY: setup format lint check test clean check-all benchmark benchmark-hdf5

setup:
	python3 -m venv .venv
	.venv/bin/pip install --upgrade pip
	.venv/bin/pip install -e . ruff

format:
	.venv/bin/ruff format src/
	.venv/bin/ruff check --fix src/

lint:
	.venv/bin/ruff check src/

run:
	.venv/bin/mgatk2 run --help
	cd tests && ../.venv/bin/mgatk2 run -o run_hdf5_output -f hdf5
	cd tests && ../.venv/bin/mgatk2 run -o run_txt_output -f txt

tenx:
	.venv/bin/mgatk2 tenx --help
	cd tests && ../.venv/bin/mgatk2 tenx -o tenx_output
	cd tests && ../.venv/bin/mgatk2 run -o run_output

benchmark:
	cd tests && python benchmark.py

benchmark-hdf5:
	cd tests && python benchmark.py --format hdf5

clean:
	rm -rf build dist *.egg-info
	rm -rf .ruff_cache .venv .mypy_cache .flake8
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type d -name *.egg-info -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	find . -type f -name ".DS_Store" -delete
