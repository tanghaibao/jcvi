import pytest


@pytest.mark.parametrize(
    "input,output", [("AT5G54690.2", "AT5G54690"), ("evm.test.1", "evm.test.1")],
)
def test_gene_name(input, output):
    from jcvi.utils.cbook import gene_name

    assert gene_name(input) == output
