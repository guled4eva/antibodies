"""The chat agent: a Claude tool-use loop over the antibody tools.

You talk to it in plain language and paste sequences; it decides when to call
`analyse_antibody` / `predict_structure`, reads the results, and explains them.
"""

from __future__ import annotations

import os
import sys

from anthropic import Anthropic

from .tools import TOOL_SCHEMAS, run_tool

MODEL = os.environ.get("ABAGENT_MODEL", "claude-opus-4-8")

SYSTEM = (
    "You are an antibody engineering assistant for a drug-discovery team. "
    "When the user gives antibody sequences, call the analysis tools rather "
    "than guessing. Read the tool output and explain, in plain terms, what it "
    "means for developability: which liabilities actually matter (a liability "
    "inside CDR3 or at a paratope residue is far more serious than one buried "
    "in framework), what the pI and charge imply for formulation and "
    "aggregation, and what you would check or change next. Be concrete and "
    "honest about uncertainty, and say when a result comes from the rough "
    "fallback engine rather than a validated one. Do not fabricate numbers."
)


def _run_tools(blocks) -> list[dict]:
    results = []
    for block in blocks:
        if block.type == "tool_use":
            output = run_tool(block.name, block.input)
            results.append({
                "type": "tool_result",
                "tool_use_id": block.id,
                "content": output,
            })
    return results


def chat_once(client: Anthropic, messages: list[dict]) -> str:
    """Advance the conversation one user turn, running tools until the model stops."""
    while True:
        resp = client.messages.create(
            model=MODEL, max_tokens=2048, system=SYSTEM,
            tools=TOOL_SCHEMAS, messages=messages,
        )
        messages.append({"role": "assistant", "content": resp.content})
        if resp.stop_reason != "tool_use":
            return "".join(b.text for b in resp.content if b.type == "text")
        messages.append({"role": "user", "content": _run_tools(resp.content)})


def _load_dotenv() -> None:
    """Load KEY=value lines from a .env in the cwd or repo root, if present."""
    for path in (".env", os.path.join(os.path.dirname(__file__), "..", ".env")):
        if os.path.isfile(path):
            for line in open(path):
                line = line.strip()
                if line and not line.startswith("#") and "=" in line:
                    k, v = line.split("=", 1)
                    os.environ.setdefault(k.strip(), v.strip().strip("'\""))
            return


def repl() -> None:
    _load_dotenv()
    if not os.environ.get("ANTHROPIC_API_KEY"):
        raise SystemExit("Set ANTHROPIC_API_KEY (export it, or put it in a .env file).")
    client = Anthropic()
    messages: list[dict] = []
    print("Antibody agent. Paste sequences or ask a question. Ctrl-D to quit.\n")
    try:
        while True:
            user = input("you> ").strip()
            if not user:
                continue
            messages.append({"role": "user", "content": user})
            print("\nagent>", chat_once(client, messages), "\n")
    except (EOFError, KeyboardInterrupt):
        print("\nbye")


def selftest() -> int:
    """Exercise the whole agent path offline: no API key, no network.

    Uses a scripted fake client that asks for the `analyse_antibody` tool and
    then replies with text, so the real tool-dispatch loop in `chat_once` runs
    end to end against the real tools.
    """
    from types import SimpleNamespace

    TRAST_VH = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
                "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")

    class FakeClient:
        """Mimics Anthropic().messages.create for one tool round-trip."""

        def __init__(self) -> None:
            self.calls = 0
            self.tool_output: str | None = None
            self.messages = SimpleNamespace(create=self._create)

        def _create(self, **kw):
            self.calls += 1
            if self.calls == 1:
                block = SimpleNamespace(type="tool_use", name="analyse_antibody",
                                        id="t1", input={"heavy": TRAST_VH, "name": "trastuzumab"})
                return SimpleNamespace(content=[block], stop_reason="tool_use")
            # capture what the tool returned (fed back in the last user message)
            last = kw["messages"][-1]["content"][0]
            self.tool_output = last["content"]
            text = SimpleNamespace(type="text", text="Analysis complete; see tool output above.")
            return SimpleNamespace(content=[text], stop_reason="end_turn")

    client = FakeClient()
    messages = [{"role": "user", "content": "analyse this VH: " + TRAST_VH}]
    print("[selftest] running agent loop with a fake client (no API, no network)...")
    reply = chat_once(client, messages)

    checks = {
        "tool was called": client.calls == 2,
        "tool produced a report": client.tool_output and "Chain H" in client.tool_output,
        "report has biophysics": client.tool_output and "pI" in client.tool_output,
        "loop returned final text": bool(reply),
    }
    ok = all(checks.values())
    for name, passed in checks.items():
        print(f"  [{'PASS' if passed else 'FAIL'}] {name}")
    print("\n--- tool output the agent saw ---")
    print(client.tool_output)
    print("\n[selftest]", "ALL PASS" if ok else "FAILURES — see above")
    return 0 if ok else 1


def main() -> int:
    if "--selftest" in sys.argv:
        return selftest()
    repl()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
