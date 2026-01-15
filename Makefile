.PHONY: setup format lint check test clean

install:
	python3 -m venv .venv
	.venv/bin/pip install --upgrade pip
	.venv/bin/pip install -e . ruff mypy isort black flake8

format:
	.venv/bin/isort src/
	.venv/bin/black src/

lint:
	.venv/bin/ruff check --fix src/

check: format lint

run:
	cd tests && ../.venv/bin/mgatk2 run -o run_output

tenx:
	cd tests && ../.venv/bin/mgatk2 tenx -o tenx_output

clean:
	rm -rf build dist *.egg-info
	rm -rf .ruff_cache .venv .mypy_cache .flake8
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type d -name *.egg-info -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	find . -type f -name ".DS_Store" -delete
