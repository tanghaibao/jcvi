from jcvi.formats.sizes import Sizes

FASTA_CONTENTS = """
>chr1
ACGT
>chr2
ACGTTCGA
>chr3
ACGTTCGAGAGA
"""


def test_sizes():
    from jcvi.apps.base import cleanup
    from jcvi.formats.base import write_file

    fastafile = "test.fa"
    write_file(fastafile, FASTA_CONTENTS, skipcheck=True)
    sizes = Sizes(fastafile)
    assert len(sizes) == 3
    assert sizes.sizes_mapping["chr1"] == 4
    assert sizes.sizes_mapping["chr2"] == 8
    assert sizes.sizes_mapping["chr3"] == 12
    assert sizes.cumsizes_mapping["chr1"] == 0
    assert sizes.cumsizes_mapping["chr2"] == 4
    assert sizes.cumsizes_mapping["chr3"] == 12
    cleanup(fastafile, sizes.filename)
