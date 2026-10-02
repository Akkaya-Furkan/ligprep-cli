import pytest

from ligprep.errors import LigprepError
from ligprep.inputs import parse_line, read_list


def test_bare_text_is_a_name():
    e = parse_line("aspirin")
    assert (e.kind, e.value) == ("name", "aspirin")


def test_bare_integer_is_a_cid():
    assert parse_line("2244").kind == "cid"


@pytest.mark.parametrize("line,kind,value", [
    ("cid:3672", "cid", "3672"),
    ("name:acetic acid", "name", "acetic acid"),
    ("smiles:CC(=O)O", "smiles", "CC(=O)O"),
    ("SMILES: CC(=O)O ", "smiles", "CC(=O)O"),
])
def test_prefixes(line, kind, value):
    e = parse_line(line)
    assert (e.kind, e.value) == (kind, value)


@pytest.mark.parametrize("line", ["", "   ", "# only a comment"])
def test_blank_and_comment_lines_are_skipped(line):
    assert parse_line(line) is None


def test_trailing_comment_is_stripped():
    assert parse_line("aspirin  # painkiller").value == "aspirin"


def test_cid_prefix_requires_an_integer():
    with pytest.raises(LigprepError):
        parse_line("cid:abc")


def test_label_is_filesystem_safe():
    assert parse_line("name:acetic acid / USP").label == "acetic_acid_USP"
    assert parse_line("smiles:CC(=O)O").label == "smiles"


def test_read_list_reports_the_line_number(tmp_path):
    f = tmp_path / "bad.txt"
    f.write_text("aspirin\ncid:xyz\n")
    with pytest.raises(LigprepError, match=r"bad\.txt:2"):
        read_list(f)


def test_read_list_keeps_order(tmp_path):
    f = tmp_path / "l.txt"
    f.write_text("# header\naspirin\n\n2519\n")
    assert [e.value for e in read_list(f)] == ["aspirin", "2519"]
