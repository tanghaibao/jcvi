import pytest


@pytest.mark.parametrize(
    "blastfile,pctid,hitlen,inverse,expected",
    [
        ("blastfile", 90, 100, False, "blastfile.P90L100"),
        ("blastfile", 98.22, 0, True, "blastfile.P98_2L0.inverse"),
        ("blastfile", 99.0, 33, True, "blastfile.P99L33.inverse"),
        ("a.b_dd", 100, 5, False, "a.b_dd.P100L5"),
    ],
)
def test_filtered_blastfile_name(
    blastfile: str,
    pctid: float,
    hitlen: int,
    inverse: bool,
    expected: str,
):
    from jcvi.formats.blast import filtered_blastfile_name

    assert filtered_blastfile_name(blastfile, pctid, hitlen, inverse) == expected
