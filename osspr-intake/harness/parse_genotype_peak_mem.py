#!/usr/bin/env python3
"""Peak Python allocations while parsing George Church's 23andMe file.

Uses ``tracemalloc`` around one ``genome_compare._parse_genotype_file`` call
on ``skills/genome-compare/data/george_church_23andme.txt.gz``. Tracing starts
after imports. The checksum is computed after tracing stops, so it is not
part of the peak.

    python osspr-intake/harness/parse_genotype_peak_mem.py
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import tracemalloc
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def load_genome_compare():
    name = "osspr_genome_compare_peak"
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "skills" / "genome-compare" / "genome_compare.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    sys.path.insert(0, str(ROOT))
    spec.loader.exec_module(module)
    return module


def checksum(genotypes: dict, positions: dict) -> str:
    digest = hashlib.sha256()
    for rsid in sorted(genotypes):
        pos = positions[rsid]
        digest.update(rsid.encode())
        digest.update(b"\0")
        digest.update(genotypes[rsid].encode())
        digest.update(b"\0")
        digest.update(str(pos["chrom"]).encode())
        digest.update(b"\0")
        digest.update(str(pos["pos"]).encode())
        digest.update(b"\n")
    return digest.hexdigest()


def main() -> None:
    module = load_genome_compare()
    path = module.REFERENCE_FILE
    tracemalloc.start()
    genotypes, positions = module._parse_genotype_file(path)
    _current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    json.dump(
        {
            "peak_memory_bytes": peak,
            "n_snps": len(genotypes),
            "checksum": checksum(genotypes, positions),
        },
        sys.stdout,
    )
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
