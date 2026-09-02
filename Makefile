.PHONY: install design report test lint

install:
	python3 -m pip install -e ".[dev,figures]"

design:
	cosmepep design

report:
	cosmepep report

test:
	python3 -m pytest tests -q

lint:
	python3 -m ruff check src tests
	python3 -m mypy src
