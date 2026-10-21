# PyInstaller spec for abagent. Build with:  pyinstaller packaging/abagent.spec
# Run this on the OS you want the binary for (Windows -> abagent.exe).
from PyInstaller.utils.hooks import collect_submodules, copy_metadata

datas = copy_metadata("anthropic")
hiddenimports = collect_submodules("anthropic") + [
    "abagent", "abagent.agent", "abagent.tools", "abagent.developability",
    "abagent.numbering", "abagent.sequence", "abagent.report", "abagent.structure",
]

a = Analysis(
    ["launch.py"],
    pathex=[".."],          # repo root, so the abagent package is importable
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    runtime_hooks=[],
    excludes=["torch", "ImmuneBuilder", "anarcii"],  # heavy optional engines
    noarchive=False,
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz, a.scripts, a.binaries, a.datas, [],
    name="abagent",
    console=True,
    onefile=True,
    upx=False,
)
