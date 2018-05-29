import pytest

def pytest_addoption(parser):
    parser.addoption("--make_interactive", action="store_true",
        help="in interactive mode, some plots are shown and their "
             "correctness is queried via input()")

@pytest.fixture
def make_interactive(request):
    return request.config.getoption("--make_interactive")