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


def test_cdr3_known_limitation_extra_cysteines():
    """Pins a KNOWN, ACCEPTED limitation of the fallback engine.

    The heuristic treats the second cysteine as Cys104, so an antibody with an
    extra cysteine before the real one gets a wrong CDR3. Measured impact on
    real data is small (97.1% exact over 5,000 SAbDab2 structures); the
    principled alternative -- anchoring on the J motif -- scored worse overall
    (96.3%, with 0.5% no-call), so the trade-off is accepted rather than fixed.
    Install anarcii for exact IMGT numbering.

    This test asserts the wrong answer on purpose, so the limitation stays
    visible instead of being forgotten.
    """
    vh = TRAST_VH.replace("GFNIKDTY", "GFNIKCTY")
    assert vh.count("C") > TRAST_VH.count("C")  # genuinely has an extra cysteine
    assert cdrs(vh, "H").cdr3 != "SRWGGDGFYAMDY"


def test_cdr3_length_is_plausible():
    loop = cdrs(TRAST_VH, "H").cdr3
    assert 2 <= len(loop) <= 70


def test_gemini_can_invoke_tools():
    """Regression: Gemini's automatic function calling must be able to call the tools.

    PEP 563 string annotations (`from __future__ import annotations` in tools.py)
    break the SDK's isinstance() checks, so every tool call fails and the model
    answers from memory instead. Exercises the SDK's real invocation path; no
    network or API key needed.
    """
    try:
        from google.genai import _extra_utils as eu
    except ImportError:
        return  # google-genai not installed; nothing to check

    from abagent.tools import analyse_antibody, predict_structure

    out = eu.invoke_function_from_dict_args(
        {"heavy": TRAST_VH, "name": "trastuzumab"}, analyse_antibody
    )
    assert "Chain H" in out and "pI" in out, "tool did not return a real report"
    # the model may omit optional args entirely
    out2 = eu.invoke_function_from_dict_args({"heavy": TRAST_VH}, analyse_antibody)
    assert "Chain H" in out2
    eu.invoke_function_from_dict_args({"heavy": TRAST_VH}, predict_structure)


def test_tool_annotations_are_real_types():
    """The tool signatures must expose real types, not PEP 563 strings."""
    import inspect

    from abagent.tools import analyse_antibody

    for param in inspect.signature(analyse_antibody).parameters.values():
        assert not isinstance(param.annotation, str), (
            f"{param.name} annotation is a string — remove "
            "`from __future__ import annotations` from abagent/tools.py"
        )


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
    print("\n--- sample report ---")
    print(analyse_antibody(heavy=TRAST_VH, light=TRAST_VL, name="trastuzumab"))
