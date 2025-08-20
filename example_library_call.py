
from pathlib import Path
import sys
from strump_i.__main__ import main as strumpi_main

def run_modeler(allele: str, pepfile: str, outdir: str):
    argv = ["Modeler", "--allele", allele, "--pepfile", pepfile, "--output_dir", outdir]
    saved = sys.argv
    try:
        sys.argv = ["strumpi"] + argv
        strumpi_main()
    finally:
        sys.argv = saved

def run_predictor(input_dir: str, output_dir: str):
    argv = ["BindingPredictor", "--input_dir", input_dir, "--output_dir", output_dir]
    saved = sys.argv
    try:
        sys.argv = ["strumpi"] + argv
        strumpi_main()
    finally:
        sys.argv = saved

if __name__ == "__main__":
    run_modeler("HLA-A*02:01", "/path/to/peptides.txt", "/abs/path/to/out")
    run_predictor("/abs/path/to/out", "/abs/path/to/out/predictions")
