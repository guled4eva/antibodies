"""abagent: an AI agent for antibody developability analysis.

Built for the AI x Antibodies Hackathon (Oslo, 14 Nov 2026).
"""

from .developability import analyse_chain
from .numbering import cdrs
from .sequence import Antibody, SequenceError
from .tools import analyse_antibody, predict_structure

__version__ = "0.1.0"
__all__ = [
    "Antibody", "SequenceError", "analyse_chain", "cdrs",
    "analyse_antibody", "predict_structure",
]
