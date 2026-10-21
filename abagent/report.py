"""Turn analysis dataclasses into compact, readable text for the agent."""

from __future__ import annotations

from .developability import ChainReport
from .numbering import CDRSet


def format_chain(report: ChainReport, cdr: CDRSet | None) -> str:
    lines = [
        f"Chain {report.chain}  ({report.length} aa)",
        f"  pI {report.pI} | net charge @pH7 {report.net_charge_pH7:+.2f} | GRAVY {report.gravy}",
    ]
    if cdr:
        lines.append(f"  CDRs ({cdr.engine}): H1/L1={cdr.cdr1 or '-'}  "
                     f"H2/L2={cdr.cdr2 or '-'}  H3/L3={cdr.cdr3 or '-'}")
    if report.liabilities:
        lines.append(f"  Sequence liabilities ({len(report.liabilities)}):")
        for liab in report.liabilities:
            in_cdr3 = cdr and liab.motif and cdr.cdr3 and _overlaps(report, liab, cdr)
            flag = "  <-- in CDR3" if in_cdr3 else ""
            lines.append(f"    pos {liab.position:>3}  {liab.kind:<20} ({liab.note}){flag}")
    else:
        lines.append("  No sequence liabilities flagged.")
    return "\n".join(lines)


def _overlaps(report: ChainReport, liab, cdr: CDRSet) -> bool:
    # cheap check: is the liability motif a substring of the CDR3 region
    return bool(cdr.cdr3) and liab.motif[:2] in cdr.cdr3
