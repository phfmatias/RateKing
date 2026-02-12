.PHONY: clean install dev-install

clean:
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	rm -rf build/ dist/ *.egg-info .pytest_cache .coverage htmlcov/

install: clean
	pip install .

dev-install: clean
	pip install -e .

reinstall: clean
	pip uninstall -y RateKing
	pip install -e .
