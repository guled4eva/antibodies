"""LLM backends behind one interface, so the same antibody tools can be driven
by either Claude or Gemini.

Pick the provider with the ABAGENT_PROVIDER env var ("anthropic" or "gemini").
Each backend keeps its own conversation history and exposes `.send(text) -> str`.
The tools themselves (abagent.tools) are provider-agnostic.
"""

from __future__ import annotations

import os

from .tools import TOOL_SCHEMAS, analyse_antibody, predict_structure, run_tool

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

# Python callables the models may call, by tool name.
TOOL_FUNCS = [analyse_antibody, predict_structure]

DEFAULT_MODELS = {
    "anthropic": "claude-opus-4-8",
    "gemini": "gemini-2.0-flash",
}


# --------------------------------------------------------------------------- #
# Anthropic
# --------------------------------------------------------------------------- #
def _run_tools(blocks) -> list[dict]:
    results = []
    for block in blocks:
        if block.type == "tool_use":
            results.append({
                "type": "tool_result",
                "tool_use_id": block.id,
                "content": run_tool(block.name, block.input),
            })
    return results


def chat_once(client, messages: list[dict], model: str) -> str:
    """One user turn against Claude, running tools until the model stops."""
    while True:
        resp = client.messages.create(
            model=model, max_tokens=2048, system=SYSTEM,
            tools=TOOL_SCHEMAS, messages=messages,
        )
        messages.append({"role": "assistant", "content": resp.content})
        if resp.stop_reason != "tool_use":
            return "".join(b.text for b in resp.content if b.type == "text")
        messages.append({"role": "user", "content": _run_tools(resp.content)})


class AnthropicBackend:
    name = "anthropic"

    def __init__(self, model: str | None = None) -> None:
        from anthropic import Anthropic

        if not os.environ.get("ANTHROPIC_API_KEY"):
            raise RuntimeError("Set ANTHROPIC_API_KEY (export it or put it in .env).")
        self.client = Anthropic()
        self.model = model or os.environ.get("ABAGENT_MODEL") or DEFAULT_MODELS["anthropic"]
        self.messages: list[dict] = []

    def send(self, text: str) -> str:
        self.messages.append({"role": "user", "content": text})
        return chat_once(self.client, self.messages, self.model)


# --------------------------------------------------------------------------- #
# Gemini (automatic function calling: the SDK introspects the Python tools)
# --------------------------------------------------------------------------- #
class GeminiBackend:
    name = "gemini"

    def __init__(self, model: str | None = None) -> None:
        from google import genai
        from google.genai import types

        key = os.environ.get("GEMINI_API_KEY") or os.environ.get("GOOGLE_API_KEY")
        if not key:
            raise RuntimeError("Set GEMINI_API_KEY (export it or put it in .env).")
        self.client = genai.Client(api_key=key)
        self.model = model or os.environ.get("ABAGENT_MODEL") or DEFAULT_MODELS["gemini"]
        config = types.GenerateContentConfig(
            system_instruction=SYSTEM,
            tools=TOOL_FUNCS,  # SDK builds declarations + auto-executes these
        )
        self.chat = self.client.chats.create(model=self.model, config=config)

    def send(self, text: str) -> str:
        resp = self.chat.send_message(text)
        return resp.text or "(no text returned)"


# --------------------------------------------------------------------------- #
def get_backend(provider: str | None = None):
    """Return a backend instance for the chosen provider (default: anthropic)."""
    provider = (provider or os.environ.get("ABAGENT_PROVIDER") or "anthropic").lower()
    if provider == "anthropic":
        return AnthropicBackend()
    if provider == "gemini":
        return GeminiBackend()
    raise RuntimeError(f"Unknown provider {provider!r}; use 'anthropic' or 'gemini'.")
