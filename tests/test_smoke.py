"""Smoke test on trastuzumab (Herceptin) VH/VL."""

from abagent.developability import analyse_chain
from abagent.numbering import cdrs
from abagent.sequence import Antibody, SequenceError, clean
from abagent.tools import analyse_antibody

# Trastuzumab variable regions (public sequences)
TRAST_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
TRAST_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR"
    "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)


def test_clean_rejects_bad():
    try:
        clean("EVQZZ")
    except SequenceError:
        return
    raise AssertionError("expected SequenceError")


def test_antibody_both_chains():
    ab = Antibody(heavy=TRAST_VH, light=TRAST_VL, name="trastuzumab")
    assert set(ab.chains()) == {"H", "L"}


def test_developability_runs():
    rep = analyse_chain(TRAST_VH, "H")
    assert rep.length == len(TRAST_VH)
    assert 4.0 < rep.pI < 10.0
    # trastuzumab VH has a known N-glyc-free but liability-bearing CDR set
    assert any("Deamidation" in liab.kind for liab in rep.liabilities)


def test_cdr3_found():
    loops = cdrs(TRAST_VH, "H")
    assert loops.cdr3  # should find the WGGDGFYAMDY-region H3


def test_tool_text():
    out = analyse_antibody(heavy=TRAST_VH, light=TRAST_VL, name="trastuzumab")
    assert "Chain H" in out and "Chain L" in out
    assert "pI" in out


def test_dotenv_strips_inline_comments():
    """Regression: `KEY=value  # comment` must not keep the comment in the value."""
    import os
    import tempfile

    from abagent.agent import _load_dotenv

    d = tempfile.mkdtemp()
    cwd = os.getcwd()
    try:
        os.chdir(d)
        with open(".env", "w") as fh:
            fh.write("ABAGENT_PROVIDER=gemini          # or: anthropic\n")
        os.environ.pop("ABAGENT_PROVIDER", None)
        _load_dotenv()
        assert os.environ["ABAGENT_PROVIDER"] == "gemini"
    finally:
        os.chdir(cwd)
        os.environ.pop("ABAGENT_PROVIDER", None)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
    print("\n--- sample report ---")
    print(analyse_antibody(heavy=TRAST_VH, light=TRAST_VL, name="trastuzumab"))
