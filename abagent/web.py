"""A small local web UI for the antibody quality checker.

Deliberately standard-library only: no Flask, FastAPI or uvicorn. The analysis
is the valuable part and it already works; a demo should not be able to fail
because a web framework would not install.

Run:
    python -m abagent.web            # http://localhost:8000
    python -m abagent.web --port 8080 --no-browser
"""

import argparse
import json
import socket
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from .developability import analyse_chain
from .numbering import cdrs
from .sequence import Antibody, SequenceError

STATIC_DIR = Path(__file__).parent / "static"


def analyse_payload(heavy="", light="", name="query"):
    """Run the analysis and return a JSON-serialisable report."""
    try:
        ab = Antibody(heavy=heavy or None, light=light or None, name=name)
    except SequenceError as exc:
        return {"error": str(exc)}

    chains = []
    for chain, seq in ab.chains().items():
        rep = analyse_chain(seq, chain)
        loops = cdrs(seq, chain)
        in_cdr3 = set()
        if loops.cdr3:
            start = seq.find(loops.cdr3)
            if start >= 0:
                in_cdr3 = set(range(start + 1, start + 1 + len(loops.cdr3)))
        chains.append({
            "chain": chain,
            "sequence": seq,
            "length": rep.length,
            "pI": rep.pI,
            "net_charge": rep.net_charge_pH7,
            "gravy": rep.gravy,
            "cdr": {"engine": loops.engine, "cdr1": loops.cdr1,
                    "cdr2": loops.cdr2, "cdr3": loops.cdr3},
            "liabilities": [
                {"kind": liability.kind, "position": liability.position,
                 "motif": liability.motif, "note": liability.note,
                 "in_cdr3": liability.position in in_cdr3}
                for liability in rep.liabilities
            ],
        })
    return {"name": ab.name, "chains": chains}


class Handler(BaseHTTPRequestHandler):
    server_version = "abagent"

    def log_message(self, fmt, *args):  # quieter console
        return

    def _send(self, code, body, content_type="application/json"):
        if isinstance(body, (dict, list)):
            body = json.dumps(body).encode()
        elif isinstance(body, str):
            body = body.encode()
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        path = self.path.split("?", 1)[0]
        if path in ("/", "/index.html"):
            page = STATIC_DIR / "index.html"
            if not page.exists():
                return self._send(500, "index.html missing", "text/plain")
            return self._send(200, page.read_text(encoding="utf-8"),
                              "text/html; charset=utf-8")
        if path == "/api/example":
            return self._send(200, {
                "name": "trastuzumab",
                "heavy": ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPT"
                          "NGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDY"
                          "WGQGTLVTVSS"),
                "light": ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLY"
                          "SGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"),
            })
        return self._send(404, {"error": "not found"})

    def do_POST(self):
        if self.path.split("?", 1)[0] != "/api/analyse":
            return self._send(404, {"error": "not found"})
        try:
            length = int(self.headers.get("Content-Length") or 0)
            data = json.loads(self.rfile.read(length) or b"{}")
        except (ValueError, TypeError):
            return self._send(400, {"error": "invalid JSON"})
        result = analyse_payload(
            heavy=(data.get("heavy") or "").strip(),
            light=(data.get("light") or "").strip(),
            name=(data.get("name") or "query").strip() or "query",
        )
        return self._send(400 if "error" in result else 200, result)


def _free_port(preferred):
    with socket.socket() as sock:
        try:
            sock.bind(("127.0.0.1", preferred))
            return preferred
        except OSError:
            sock.bind(("127.0.0.1", 0))
            return sock.getsockname()[1]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--no-browser", action="store_true")
    args = ap.parse_args(argv)

    port = _free_port(args.port)
    if port != args.port:
        print(f"port {args.port} busy, using {port}")
    url = f"http://{args.host}:{port}"
    server = ThreadingHTTPServer((args.host, port), Handler)
    print(f"Antibody quality checker running at {url}\nCtrl-C to stop.")
    if not args.no_browser:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
