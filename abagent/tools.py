"""Tool definitions the Claude agent can call, plus their implementations.

Each tool is a small, deterministic Python function. The agent decides which
to call; the functions do the actual antibody analysis and return plain text
so results stay inspectable.
"""

from __future__ import annotations

from .developability import analyse_chain
from .numbering import cdrs
from .report import format_chain
from .sequence import Antibody, SequenceError
from .structure import predict_fv

# --- JSON schemas advertised to the model ------------------------------------
TOOL_SCHEMAS = [
    {
        "name": "analyse_antibody",
        "description": (
            "Analyse an antibody variable region: number it, pull out the CDR "
            "loops, and produce a developability report (sequence liabilities, "
            "isoelectric point, net charge, hydrophobicity). Give heavy and/or "
            "light chain amino-acid sequences."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "heavy": {"type": "string", "description": "Heavy-chain VH amino-acid sequence (optional)."},
                "light": {"type": "string", "description": "Light-chain VL amino-acid sequence (optional)."},
                "name": {"type": "string", "description": "Label for this antibody."},
            },
        },
    },
    {
        "name": "predict_structure",
        "description": (
            "Predict the 3D Fv structure of the antibody (ABodyBuilder2) and "
            "write a PDB file. Optional heavy engine; reports clearly if not installed."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "heavy": {"type": "string"},
                "light": {"type": "string"},
            },
        },
    },
]


def analyse_antibody(heavy: str = "", light: str = "", name: str = "query") -> str:
    try:
        ab = Antibody(heavy=heavy or None, light=light or None, name=name)
    except SequenceError as exc:
        return f"Error: {exc}"
    blocks = [f"Antibody: {ab.name}"]
    for chain, seq in ab.chains().items():
        rep = analyse_chain(seq, chain)
        loops = cdrs(seq, chain)
        blocks.append(format_chain(rep, loops))
    return "\n\n".join(blocks)


def predict_structure(heavy: str = "", light: str = "") -> str:
    res = predict_fv(heavy or None, light or None)
    return res.message


DISPATCH = {
    "analyse_antibody": analyse_antibody,
    "predict_structure": predict_structure,
}


def run_tool(name: str, tool_input: dict) -> str:
    fn = DISPATCH.get(name)
    if fn is None:
        return f"Unknown tool: {name}"
    try:
        return fn(**tool_input)
    except Exception as exc:  # keep the agent loop alive on bad input
        return f"Tool {name} failed: {exc}"
