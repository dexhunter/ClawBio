#!/usr/bin/env python3
"""Time dnasp.compute_k on the rp49 sliding-window slices.

The windows are the ingroup of the committed DnaSP fixture
``skills/dnasp/tests/fixtures/inputs/rp49_ing34_outGUA.fas`` (outgroup
``rp49.gua`` removed), cut with the same 100 bp / 25 bp loop as
``run_analysis``. Each slice is passed through ``complete_deletion`` before
the timer, which is what ``analyse_region`` does before it calls
``compute_k``. Empty slices never reach ``compute_k`` and are not timed.

Prints one JSON object: seconds, call count, and a checksum of the returned
floats. Run from the repository root:

    python osspr-intake/harness/compute_k_rp49_windows.py
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import struct
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WINDOW = 100
STEP = 25
FIXTURE = (
    ROOT
    / "skills"
    / "dnasp"
    / "tests"
    / "fixtures"
    / "inputs"
    / "rp49_ing34_outGUA.fas"
)
OUTGROUP = "rp49.gua"


def load_dnasp():
    name = "osspr_dnasp_compute_k"
    spec = importlib.util.spec_from_file_location(name, ROOT / "skills" / "dnasp" / "dnasp.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def window_slices(seqs: list[str]) -> list[tuple[int, int]]:
    """(start, end) pairs, matching run_analysis's window loop."""
    spans: list[tuple[int, int]] = []
    length = len(seqs[0])
    pos = 0
    while pos < length:
        end = min(pos + WINDOW, length)
        spans.append((pos, end))
        if end == length:
            break
        pos = min(pos + STEP, length - 1)
    return spans


def cleaned_windows(dn, seqs: list[str]) -> list[list[str]]:
    prepared: list[list[str]] = []
    for start, end in window_slices(seqs):
        raw = [seq[start:end] for seq in seqs]
        clean, net = dn.complete_deletion(raw)
        if len(clean) >= 2 and net > 0:
            prepared.append(clean)
    return prepared


def checksum(values: list[float]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(struct.pack("<d", value))
    return digest.hexdigest()


def main() -> None:
    dn = load_dnasp()
    alignment = dn.load_alignment(FIXTURE)
    kept = [(name, seq) for name, seq in zip(alignment.names, alignment.seqs) if name != OUTGROUP]
    seqs = [seq for _, seq in kept]
    windows = cleaned_windows(dn, seqs)
    started = time.perf_counter()
    values = [dn.compute_k(window) for window in windows]
    elapsed = time.perf_counter() - started
    json.dump(
        {
            "seconds": elapsed,
            "calls": len(values),
            "checksum": checksum(values),
        },
        sys.stdout,
    )
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
