"""The chat agent REPL and an offline selftest.

You talk to it in plain language and paste sequences; the chosen backend
(Claude or Gemini) decides when to call the antibody tools, reads the results,
and explains them. Pick the model with ABAGENT_PROVIDER=anthropic|gemini.
"""

from __future__ import annotations

import os
import sys

from .backends import chat_once, get_backend


def _load_dotenv() -> None:
    """Load KEY=value lines from a .env in the cwd or repo root, if present."""
    for path in (".env", os.path.join(os.path.dirname(__file__), "..", ".env")):
        if os.path.isfile(path):
            for line in open(path):
                line = line.strip()
                if line and not line.startswith("#") and "=" in line:
                    k, v = line.split("=", 1)
                    # drop inline comments:  VALUE   # note  ->  VALUE
                    for sep in (" #", "\t#"):
                        v = v.split(sep, 1)[0]
                    os.environ.setdefault(k.strip(), v.strip().strip("'\""))
            return


def _friendly_error(exc: Exception) -> str:
    """Turn a provider SDK exception into one actionable sentence."""
    msg = str(exc)
    low = msg.lower()
    if any(s in low for s in ("api key not valid", "api_key_invalid", "invalid api key",
                              "unauthenticated", "unauthorized", "401", "403",
                              "permission_denied")):
        return ("The API key was rejected. Check the key in .env matches the provider: "
                "Gemini keys come from https://aistudio.google.com/apikey, "
                "Claude keys from https://console.anthropic.com.")
    if any(s in low for s in ("quota", "resource_exhausted", "rate limit", "429")):
        return "Rate limit or quota reached. Wait a minute, then try again."
    if "model" in low and ("not found" in low or "404" in low):
        return ("That model name isn't available (Google retires old ones). The error "
                "usually names the replacement — set it, e.g. "
                "ABAGENT_MODEL=gemini-3.8-flash")
    if any(s in low for s in ("503", "unavailable", "high demand", "overloaded")):
        return ("The model is busy on the provider's side (not your setup) — this is "
                "temporary. Try again in a minute, or switch model/provider, e.g. "
                "ABAGENT_PROVIDER=anthropic.")
    if any(s in low for s in ("connection", "timeout", "network", "getaddrinfo")):
        return "Could not reach the API — looks like a network problem. Check your connection."
    return f"{type(exc).__name__}: {msg[:300]}"


def repl() -> None:
    _load_dotenv()
    backend = get_backend()  # raises a clear error if the key is missing
    print(f"Antibody agent [{backend.name}:{backend.model}]. "
          "Paste sequences or ask a question. Ctrl-D to quit.\n")
    try:
        while True:
            user = input("you> ").strip()
            if not user:
                continue
            try:
                print("\nagent>", backend.send(user), "\n")
            except Exception as exc:
                print("\n[error]", _friendly_error(exc), "\n")
    except (EOFError, KeyboardInterrupt):
        print("\nbye")


def selftest() -> int:
    """Exercise the agent path offline: no API key, no network.

    1. Runs the real Anthropic tool-dispatch loop against a scripted fake client.
    2. Checks the Gemini backend wires the tools up (if google-genai is installed),
       without making any network call.
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
            self.tool_output = kw["messages"][-1]["content"][0]["content"]
            text = SimpleNamespace(type="text", text="Analysis complete; see tool output above.")
            return SimpleNamespace(content=[text], stop_reason="end_turn")

    checks: dict[str, bool] = {}

    # 1. Anthropic loop
    client = FakeClient()
    messages = [{"role": "user", "content": "analyse this VH: " + TRAST_VH}]
    print("[selftest] anthropic loop with fake client (no API, no network)...")
    reply = chat_once(client, messages, model="fake")
    checks["anthropic: tool was called"] = client.calls == 2
    checks["anthropic: tool produced a report"] = bool(client.tool_output) and "Chain H" in client.tool_output
    checks["anthropic: report has biophysics"] = bool(client.tool_output) and "pI" in client.tool_output
    checks["anthropic: loop returned text"] = bool(reply)

    # 2. Gemini wiring (construction only, dummy key, no send -> no network)
    print("[selftest] gemini backend wiring (dummy key, no network)...")
    try:
        from .backends import GeminiBackend
        os.environ.setdefault("GEMINI_API_KEY", "selftest-dummy")
        gb = GeminiBackend()
        checks["gemini: backend constructs"] = gb.name == "gemini"
        checks["gemini: chat session created"] = gb.chat is not None

        # Actually invoke a tool the way the SDK does during automatic function
        # calling. Construction succeeding does NOT mean tool calls work.
        from google.genai import _extra_utils as _eu

        from .tools import analyse_antibody as _aa

        out = _eu.invoke_function_from_dict_args(
            {"heavy": TRAST_VH, "name": "trastuzumab"}, _aa
        )
        checks["gemini: tool actually invokes"] = "Chain H" in out and "pI" in out
    except ImportError:
        print("  (google-genai not installed; skipping gemini checks)")
    except Exception as exc:  # construction should not hit the network
        checks["gemini: backend constructs"] = False
        print(f"  gemini wiring error: {exc}")

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
