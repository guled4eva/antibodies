"""Optional Fv structure prediction.

Wraps ImmuneBuilder/ABodyBuilder2 when installed. That model is heavy, so
this is a thin, optional wrapper: if the package is missing it returns a
clear "not available" result instead of raising, and hackathon day can
`pip install ImmuneBuilder` to turn it on.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class StructureResult:
    available: bool
    pdb_path: str | None
    message: str


def predict_fv(heavy: str | None, light: str | None,
               out_path: str = "fv_model.pdb") -> StructureResult:
    try:
        from ImmuneBuilder import ABodyBuilder2
    except Exception:
        return StructureResult(
            available=False, pdb_path=None,
            message="Structure engine not installed. Run `pip install ImmuneBuilder` "
                    "to enable Fv prediction.",
        )
    try:
        predictor = ABodyBuilder2()
        seqs = {}
        if heavy:
            seqs["H"] = heavy
        if light:
            seqs["L"] = light
        antibody = predictor.predict(seqs)
        antibody.save(out_path)
        return StructureResult(available=True, pdb_path=out_path,
                               message=f"Fv model written to {out_path}")
    except Exception as exc:  # pragma: no cover - depends on optional engine
        return StructureResult(available=False, pdb_path=None,
                               message=f"Structure prediction failed: {exc}")
