# STRUMP-I

Structure-based peptide/MHC-I prediction pipeline packaged as a Python module with a CLI.

## Install from wheel

1. Build the wheel (see below) or obtain one from `dist/`.
2. Install into your environment:
   ```bash
   pip install dist/strump_i-*.whl
   ```

## Build the wheel

From the repository root (where `pyproject.toml` lives):

```bash
python build_wheel.py
# or, if you have the build tool:
python -m build -w
```

The wheel will be created in `dist/`.

## Command-line usage

After installation, the `strumpi` command is available.

### 1) Modeler

Build pMHC models for a list of peptides against a single allele.

```bash
strumpi Modeler       --allele HLA-A*02:01       --pepfile /path/to/peptides.txt       --output_dir /abs/path/to/out       --tinker_dir /abs/path/to/tinker       --foldx /abs/path/to/foldx
```

**Inputs**
- `--allele`: MHC-I allele (e.g., `HLA-A*02:01`)
- `--pepfile`: Text file with one 8–11mer peptide per line
- `--output_dir`: Output directory (created if missing)
- Optional: `--tinker_dir`, `--foldx` (paths to local installs)
- Advanced: replace defaults to scoring matrices and templates with
  `--mhc_seq_mat_fn`, `--pdb_seq_mat_fn`, `--templateinfo_fn`, `--blosum_fn`, `--pam_fn`, `--keyfile_fn`

### 2) BindingPredictor

Use the pre-trained LightGBM classifier to score modeled structures.

```bash
strumpi BindingPredictor       --input_dir /abs/path/to/out       --output_dir /abs/path/to/out/predictions
```

- Optional: `--model_dir` (defaults to the packaged reference model in `strump_i/reference/model/binding`)

### 3) SASACalculator (optional)

Compute solvent-accessible surface area (SASA) for built PDBs. Requires `freesasa` (install with `pip install freesasa` or via your package manager).

```bash
strumpi SASACalculator       --input_dir /abs/path/to/out       --threads 8
```

## Library usage

You can also call STRUMP-I from Python.

```python
from pathlib import Path
from strump_i.__main__ import main as strumpi_main

# Run the "Modeler" subcommand programmatically
import sys
argv = [
    "Modeler",
    "--allele", "HLA-A*02:01",
    "--pepfile", "/path/to/peptides.txt",
    "--output_dir", "/abs/path/to/out"
]
# Temporarily inject args
saved = sys.argv
try:
    sys.argv = ["strumpi"] + argv
    strumpi_main()
finally:
    sys.argv = saved
```

Alternatively, you may use the lower-level APIs (subject to change) in
`strump_i/functions/*` if you need fine-grained control over modeling,
energy calculation, and prediction.

## Data files

Reference data (templates, matrices, models) are packaged under
`strump_i/reference/` and installed with the wheel. Use absolute paths
when specifying output directories.

## Reproducibility tips

- Use Python >= 3.9
- Create a fresh virtualenv/conda env
- Some steps depend on external tools (`FoldX`, `Tinker`, `freesasa`).
  Install them locally and pass their paths via CLI flags.
