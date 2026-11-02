#!/usr/bin/env python
# -*- coding: UTF-8 -*-


def test_blastfilter_main(tmp_path, monkeypatch):
    from jcvi.compara.blastfilter import main

    monkeypatch.chdir(tmp_path)
    for name in ("a", "b"):
        with open(f"{name}.bed", "w") as fw:
            for i in range(4):
                print(f"chr1\t{i * 1000}\t{i * 1000 + 500}\t{name}{i}", file=fw)
    hit = "{}\t{}\t90\t100\t1\t0\t1\t100\t1\t100\t1e-50\t{}\n"
    with open("a.b.last", "w") as fw:
        fw.write(hit.format("a0", "b0", 50))
        fw.write(hit.format("a0", "b0", 100))  # best HSP of the pair is kept
        fw.write(hit.format("a0", "b3", 60))  # cscore 0.6 < 0.7
        fw.write(hit.format("a3", "b3", 80))
        fw.write(hit.format("a2", "b2", 90))  # excluded below
    with open("excl.anchors", "w") as fw:
        fw.write("###\nb2\ta2\t90\n")

    main(["a.b.last", "--qbed=a.bed", "--sbed=b.bed", "--exclude=excl.anchors"])
    with open("a.b.last.filtered") as fp:
        rows = [x.split()[:2] + [x.split()[-1]] for x in fp]
    assert rows == [["a0", "b0", "100"], ["a3", "b3", "80"]]
