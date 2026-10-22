# abagent — antibody developability agent

Built for the **AI × Antibodies Hackathon**, Oslo, 14 November 2026.

You paste an antibody's heavy and/or light chain into a chat agent and ask
questions. The agent calls deterministic analysis tools and explains the
results the way a drug-discovery team needs them — *which* sequence liabilities
actually matter, what the charge and hydrophobicity mean for formulation, and
what to check or change next.

The point of the hackathon is "results a drug-discovery team can actually
trust", so the design keeps the science in small, inspectable Python functions
and uses Claude only to orchestrate and explain — never to invent numbers.

## What it does now

- **Numbering + CDRs** — IMGT numbering and CDR extraction (`anarcii` when
  installed; a conserved-residue heuristic as a fallback so the skeleton runs).
- **Developability report** — sequence liabilities (N-glycosylation sequons,
  Asn deamidation, Asp isomerisation, Asp-Pro fragmentation, Met/Trp oxidation,
  free cysteines), isoelectric point, net charge at pH 7, GRAVY hydrophobicity.
- **Structure (optional)** — Fv prediction via ABodyBuilder2 when
  `ImmuneBuilder` is installed.

## Install

```bash
pip install -e .                 # core: anthropic + google-genai (pure Python)
pip install -e ".[numbering]"    # add anarcii for real IMGT numbering (heavy)
pip install -e ".[structure]"    # add ImmuneBuilder for Fv structure (heavy)
export ANTHROPIC_API_KEY=...      # hackathon provides Claude credits
```

## Data

Validation / benchmarking uses the **SAbDab2 ML Dataset** (Oxford OPIG) —
16,415 curated antibody and antibody–antigen structures with train/test splits.
It is ~958 MB, so it is not committed; fetch it on demand:

```bash
python scripts/get_data.py     # downloads + verifies + extracts into ./data (gitignored)
```

The split CSVs carry ground-truth IMGT annotations (sequences, CDR loops,
numbering) for every structure, so they can validate the tools without parsing
a single `.cif`. Measure the CDR finder against them:

```bash
python scripts/benchmark_cdrs.py --limit 5000      # exact / partial / wrong
```

Licensed CC-BY-4.0. If you use it, cite:
> Capel, H. L., Vavourakis, O., Williams, B. H., Taylor, C., & Deane, C. M.
> (2026). *SAbDab2 Machine Learning Dataset* (v0.3.0) [Dataset]. Zenodo.
> https://doi.org/10.5281/zenodo.23343734

## Web UI

```bash
python -m abagent.web          # opens http://localhost:8000
```

Paste a variable region, get a report with the liabilities marked on the
sequence itself, ranked by consequence (a liability inside CDR3 is treated
differently from the same chemistry in framework). Standard library only — no
Flask, FastAPI or uvicorn — so a demo cannot fail because a web framework
would not install.

## Run

```bash
abagent          # start the chat agent
```

Then paste a sequence and ask, e.g. *"any liabilities in the CDRs I should worry
about?"*

### Choosing the model (Claude or Gemini)

The same tools run under either model; pick with `ABAGENT_PROVIDER`:

```bash
ABAGENT_PROVIDER=anthropic abagent     # default; needs ANTHROPIC_API_KEY
ABAGENT_PROVIDER=gemini    abagent     # needs GEMINI_API_KEY
```

Gemini has a free tier — get a key at https://aistudio.google.com/apikey and put
`GEMINI_API_KEY=...` in your `.env`. Override the model name with `ABAGENT_MODEL`
(defaults: `claude-opus-4-8`, `gemini-3.8-flash`).

To confirm the whole agent path works **without an API key or network** (runs
the real tool-dispatch loop against a scripted fake client):

```bash
abagent --selftest
```

## Building a standalone executable (optional)

Bundles Python + the code into one double-clickable file for teammates who
don't have Python. It still needs an `ANTHROPIC_API_KEY` and internet to chat.
PyInstaller cannot cross-compile, so **build on the OS you want the binary for**
— run this from Windows PowerShell to get `abagent.exe`:

```bash
pip install pyinstaller
pyinstaller packaging/abagent.spec      # output in dist/  (abagent.exe on Windows)
dist/abagent --selftest                 # verify the frozen build
```

The heavy optional engines (torch/anarcii/ImmuneBuilder) are excluded from the
build to keep it small (~50 MB); the frozen binary uses the fallback engines.

## Layout

| file | role |
|------|------|
| `abagent/sequence.py` | parse / validate / clean sequences |
| `abagent/numbering.py` | numbering + CDR loops (anarcii or fallback) |
| `abagent/developability.py` | liabilities, pI, charge, hydropathy (pure Python) |
| `abagent/structure.py` | optional Fv structure prediction |
| `abagent/report.py` | format analysis into readable text |
| `abagent/tools.py` | tool schemas + implementations the agent calls |
| `abagent/agent.py` | the Claude tool-use loop and REPL |

## For hackathon day

The skeleton is deliberately small so the day goes to the real challenge. Likely
next steps once the organisers announce the specific problem:

- Swap the heuristic CDR finder for `anarcii` (one `pip install`).
- Add the specific property the challenge targets (affinity, stability, …) as a
  new tool in `tools.py` — the agent picks it up automatically from `TOOL_SCHEMAS`.
- Add a `suggest_mutations` tool that proposes fixes for flagged liabilities and
  re-scores them.
- Validate against a known therapeutic (trastuzumab is used in the smoke test).

## Measured accuracy

The fallback CDR engine is scored against SAbDab2's ground-truth IMGT
annotations rather than assumed to work (`scripts/benchmark_cdrs.py`).
Headline figure is on the **held-out test split**, which SAbDab2 builds from
antibodies dissimilar to the training set:

| heavy-chain CDR3, held-out test split (n = 2,479) | |
|---|---|
| exact match | **99.0%** |
| partial overlap | 0.8% |
| wrong | 0.0% |
| no call | 0.2% |

Reproduce with:

```bash
python scripts/benchmark_cdrs.py --split test --limit 5000 --show 0
```

Notes on method, because they change how the number should be read:

- **Parameters were tuned on train, scored once on test.** Train scored 99.00%
  and test 98.99%, so the rule is not fitted to the data it is judged on.
- **Held-out matters.** Mixing train and test flatters the score by about a
  point, so the test split is what we quote.
- **Two alternatives were tried and rejected on evidence** (train split):
  anchoring on the last J motif and walking back to the nearest cysteine
  (94.5% — gives up when no WGxG/FGxG is present), and always starting at the
  last cysteine before the loop end (94.3% — CDR3s often contain their own
  cysteine, which truncates the loop).

The failure this fixed: when an antibody carries an extra cysteine before
Cys104, the primary anchor lands too early and the loop comes out ~60 residues
too long. Detecting the implausible length and re-anchoring on the conserved
`[YF]C` of framework 3 took the held-out score from 96.4% to 99.0%.

Still approximate: `pip install anarcii` for exact IMGT numbering, and CDR1/CDR2
are not called by the fallback engine at all.

## Status

Skeleton. The liability motifs and biophysical formulas are standard and
referenced in the source; numbering falls back to an approximate engine until
`anarcii` is installed. Not for clinical use.
