from jcvi import cli


def test_normalize_module_name():
    assert cli._normalize_module_name("formats.fasta") == "jcvi.formats.fasta"
    assert cli._normalize_module_name("jcvi.formats.fasta") == "jcvi.formats.fasta"


def test_main_dispatches_module(monkeypatch):
    called = {}

    def fake_run_module(module_name, run_name=None, alter_sys=False):
        called["module_name"] = module_name
        called["run_name"] = run_name
        called["alter_sys"] = alter_sys
        called["argv"] = cli.sys.argv[:]
        return {}

    monkeypatch.setattr(cli.runpy, "run_module", fake_run_module)
    rc = cli.main(["formats.fasta", "extract", "in.fa"])

    assert rc == 0
    assert called["module_name"] == "jcvi.formats.fasta"
    assert called["run_name"] == "__main__"
    assert called["alter_sys"] is True
    assert called["argv"] == ["jcvi.formats.fasta", "extract", "in.fa"]


def test_main_version_exits():
    try:
        cli.main(["--version"])
    except SystemExit as exc:
        assert exc.code == 0
    else:
        raise AssertionError("Expected --version to exit")


def test_main_no_module_prints_help(capsys):
    rc = cli.main([])
    assert rc == 0
    captured = capsys.readouterr()
    assert "Run jcvi modules" in captured.out


def test_main_missing_module_raises_import_error():
    try:
        cli.main(["does.not.exist"])
    except ImportError:
        pass
    else:
        raise AssertionError("Expected ImportError for missing module")
