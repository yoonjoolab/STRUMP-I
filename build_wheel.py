"""
Build a wheel for STRUMP-I into ./dist/
"""
import subprocess, sys, shutil, os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DIST = ROOT / "dist"
DIST.mkdir(exist_ok=True)

def have_module(name: str) -> bool:
    try:
        __import__(name)
        return True
    except Exception:
        return False

if have_module("build"):
    subprocess.check_call([sys.executable, "-m", "build", "-w", str(ROOT)])
else:
    # Fallback to legacy setuptools path
    setup_py = ROOT / "setup.py"
    if not setup_py.exists():
        setup_py.write_text("from setuptools import setup; setup()\n")
    subprocess.check_call([sys.executable, str(setup_py), "bdist_wheel", "-d", str(DIST)])
print(f"Wheel(s) written to: {DIST}")
