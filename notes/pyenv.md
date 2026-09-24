# Python environment (VERIFIED 2026-09-21)

Which interpreter runs the framework's compute scripts. Checked directly, not assumed.

## NEW MACHINE (2026-09-21): use `python3.11`

This is a Linux box, NOT the old macOS/miniconda machine. The canonical interpreter is:

```
python3.11
```

- **python3.11** (Python 3.11.15) has the full stack: `numpy 2.3.5`, `scipy 1.17.1`,
  `torch 2.12.0+cpu`, `sympy 1.14.0`. scipy IS present now (unlike the old Mac base env).
- system `python3` is **3.10.12** with `numpy 2.2.6` but a BROKEN scipy 1.8.0
  (NumPy-version mismatch warning) and NO torch. Do not use it.
- No conda on this machine.

## Consequence for the old scipy-workaround

The old Mac env had NO scipy, so `qh3_trefoil_solver_jam.py` was rewritten to a numpy
brute-force nearest-point. That still works here; scipy KDTree (`qh3_trefoil_solver_3d.py`)
also works again now. torch is CPU-only.

## Quick check command
```
python3.11 -c "import numpy,scipy,torch,sympy;print(numpy.__version__,scipy.__version__,torch.__version__,sympy.__version__)"
```

## HISTORY (old Mac machine, pre-2026-09-21, no longer valid)
Was `/usr/local/Caskroom/miniconda/base/bin/python` (torch+numpy, NO scipy); `torch_intel`
env was empty. That machine is gone; ignore those paths.
