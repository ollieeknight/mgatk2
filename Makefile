.PHONY: setup format lint check test clean

setup:
	python3 -m venv .venv
	.venv/bin/pip install --upgrade pip
	.venv/bin/pip install -e .
	.venv/bin/pip install black isort ruff flake8 mypy

format:
	.venv/bin/black src/
	.venv/bin/isort src/

lint:
	.venv/bin/ruff check src/
	.venv/bin/flake8 src/

check: format lint

test:
	cd tests && conda run -n mgatk2 mgatk2 run -o run_output

clean:
	rm -rf build dist *.egg-info
	rm -rf .flake8 .ruff_cache .venv
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
