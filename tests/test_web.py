"""Tests for the local web UI (no server framework, no network)."""

from abagent.web import analyse_payload

VH = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
      "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")


def test_payload_is_json_serialisable_and_complete():
    import json

    data = analyse_payload(heavy=VH, name="trastuzumab")
    json.dumps(data)  # must not raise
    chain = data["chains"][0]
    assert chain["chain"] == "H"
    assert chain["pI"] > 0 and chain["length"] == len(VH)
    assert chain["cdr"]["cdr3"] == "SRWGGDGFYAMDY"
    assert chain["liabilities"]


def test_cdr3_liabilities_are_flagged():
    """The UI colours by consequence, so in_cdr3 must be set correctly."""
    data = analyse_payload(heavy=VH)
    chain = data["chains"][0]
    assert any(liability["in_cdr3"] for liability in chain["liabilities"])
    # positions marked in_cdr3 must actually fall inside the called loop
    seq, loop = chain["sequence"], chain["cdr"]["cdr3"]
    start = seq.find(loop) + 1
    for liability in chain["liabilities"]:
        if liability["in_cdr3"]:
            assert start <= liability["position"] <= start + len(loop) - 1


def test_bad_sequence_returns_error_not_exception():
    assert "error" in analyse_payload(heavy="EVQZZ")
    assert "error" in analyse_payload()


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
