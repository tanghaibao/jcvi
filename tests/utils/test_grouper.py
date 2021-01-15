def test_grouper():
    from jcvi.utils.grouper import Grouper

    g = Grouper()
    g.join("a", "b")
    g.join("b", "c")
    g.join("d", "e")
    assert list(g) == [["a", "b", "c"], ["d", "e"]]
    assert g.joined("a", "b")
    assert g.joined("a", "c")
    assert "f" not in g
    assert not g.joined("a", "d")
    del g["b"]
    assert list(g) == [["a", "c"], ["d", "e"]]
