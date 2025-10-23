# CountMut Makefile
# Common development tasks

.PHONY: help install install-dev test lint format clean build publish docs

help: ## Show this help message
	@echo "CountMut Development Commands"
	@echo "============================="
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

install: ## Install package in production mode
	uv sync

install-dev: ## Install package in development mode with all dependencies
	uv sync --extra dev

test: ## Run tests
	uv run python -m pytest tests/ -v

test-cov: ## Run tests with coverage
	uv run python -m pytest tests/ -v --cov=countmut --cov-report=html --cov-report=term

test-thread: ## Run thread safety tests
	uv run python -m pytest tests/test_thread_safety.py -v

test-perf: ## Run performance tests
	uv run python -m pytest tests/test_performance.py -v

lint: ## Run linting
	uv run ruff check countmut/
	uv run ruff format --check countmut/

format: ## Format code
	uv run ruff format countmut/
	uv run ruff check --fix countmut/

format-check: ## Check code formatting
	uv run ruff format --check countmut/

clean: ## Clean build artifacts
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .pytest_cache/
	rm -rf .coverage
	rm -rf htmlcov/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

build: ## Build package
	uv build

build-check: ## Check build
	uv run twine check dist/*

publish-test: ## Publish to TestPyPI
	uv run twine upload --repository testpypi dist/*

publish: ## Publish to PyPI
	uv run twine upload dist/*

docs: ## Generate documentation
	@echo "Documentation is available in README.md"

example: ## Run example script
	uv run python examples/basic_usage.py

benchmark: ## Run performance benchmark
	uv run python -c "import time; from countmut.core import count_mutations; print('Running benchmark...'); start=time.time(); count_mutations('test.bam', 'test.fa', 'output.tsv'); print(f'Time: {time.time()-start:.2f}s')"

check: lint format-check test ## Run all checks

ci: ## Run CI pipeline locally
	uv run ruff check countmut/
	uv run ruff format --check countmut/
	uv run python -m pytest tests/ -v --cov=countmut --cov-report=xml

all: clean install-dev test lint format build ## Run full development pipeline

# Development helpers
dev-setup: ## Set up development environment
	pip install -e ".[dev]"
	pre-commit install || echo "pre-commit not available"

# Docker helpers (if needed)
docker-build: ## Build Docker image
	docker build -t countmut .

docker-test: ## Run tests in Docker
	docker run --rm countmut pytest tests/ -v

# Release helpers
version-patch: ## Bump patch version
	bump2version patch

version-minor: ## Bump minor version
	bump2version minor

version-major: ## Bump major version
	bump2version major

release: clean build publish ## Create and publish release

# Help for specific targets
test-help: ## Show test options
	@echo "Test Options:"
	@echo "  test        - Run all tests"
	@echo "  test-cov    - Run tests with coverage"
	@echo "  test-thread - Run thread safety tests"
	@echo "  test-perf   - Run performance tests"

lint-help: ## Show linting options
	@echo "Linting Options:"
	@echo "  lint         - Run linting checks"
	@echo "  format       - Format code"
	@echo "  format-check - Check code formatting"