#!/usr/bin/env python3
"""Time genome_compare._parse_genotype_file on George Church's 23andMe file.

The file is the reference shipped with the skill
(``skills/genome-compare/data/george_church_23andme.txt.gz``). The timer
covers one call, which is what ``run_comparison`` does for each genome.
The checksum is taken after the timer and does not depend on dict order.

    python osspr-intake/harness/parse_genotype_wall.py
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def load_genome_compare():
    name = "osspr_genome_compare_wall"
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
    started = time.perf_counter()
    genotypes, positions = module._parse_genotype_file(path)
    elapsed = time.perf_counter() - started
    json.dump(
        {
            "seconds": elapsed,
            "n_snps": len(genotypes),
            "checksum": checksum(genotypes, positions),
        },
        sys.stdout,
    )
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
