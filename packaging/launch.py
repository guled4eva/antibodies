"""Entry point for the frozen (PyInstaller) build of abagent.

Kept tiny and import-light so the packaged binary starts fast. It just hands
off to the normal CLI dispatcher.
"""

from abagent.agent import main

if __name__ == "__main__":
    raise SystemExit(main())
