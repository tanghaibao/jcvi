def test_tabulate():
    from jcvi.utils.table import tabulate

    data = {(1, "a"): 3, (1, "b"): 4, (2, "a"): 5, (2, "b"): 0}
    assert (
        tabulate(data)
        == """===========
o    a    b
-----------
1    3    4
2    5    0
-----------"""
    )
    assert (
        tabulate(data, transpose=True)
        == """===========
o    1    2
-----------
a    3    5
b    4    0
-----------"""
    )
