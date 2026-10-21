"""The chat agent: a Claude tool-use loop over the antibody tools.

You talk to it in plain language and paste sequences; it decides when to call
`analyse_antibody` / `predict_structure`, reads the results, and explains them.
"""

from __future__ import annotations

import os

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


def repl() -> None:
    if not os.environ.get("ANTHROPIC_API_KEY"):
        raise SystemExit("Set ANTHROPIC_API_KEY first.")
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


if __name__ == "__main__":
    repl()
