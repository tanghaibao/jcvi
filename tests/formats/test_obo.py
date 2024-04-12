def test_oboreader():
    import os
    from jcvi.formats.obo import GODag_from_GO

    go, obo_file = GODag_from_GO()
    r1, r2, r3 = [
        rec
        for i, rec in enumerate(sorted(set(go.values()), key=lambda x: x.item_id))
        if i < 3
    ]
    r4 = go["GO:0000010"]
    assert r1.item_id == "GO:0000001"
    assert r1.name == "mitochondrion inheritance"
    assert r2.item_id == "GO:0000002"
    assert r2.namespace == "biological_process"
    assert r3.item_id == "GO:0000006"
    assert tuple(sorted(r4.alt_ids)) == ("GO:0036422",)

    if os.path.exists(obo_file):
        os.remove(obo_file)
