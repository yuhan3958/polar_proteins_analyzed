"""Run the entire reproducible GH2 research pipeline from one command.

Usage:
    python run_all.py
    python run_all.py --skip-structure
    python run_all.py --skip-phylogeny

Outputs a machine-readable manifest and fails loudly on missing dependencies,
download failures, or incorrect outputs. All steps run in dependency order.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STEPS = [
    ("dataset", "run_analysis.py", "dataset_results/summary.json"),
    ("genus_audit", "scripts/audit_phylogeny.py", "dataset_results/audit/summary.json"),
    ("extended", "scripts/extended_research.py", "dataset_results/extended/summary.json"),
    ("phylogeny", "scripts/build_phylogeny.py", "dataset_results/phylogeny/summary.json"),
    ("tree_support", "scripts/tree_support_audit.py", "dataset_results/phylogeny/support_audit.json"),
    ("structure", "scripts/structure_distances.py", "dataset_results/structures/summary.json"),
    ("paper_validation", "scripts/validate_manuscript.py", "dataset_results/manuscript_validation.json"),
    ("manuscript", "scripts/generate_manuscript.py", "paper/gh2_exploratory_manuscript_ko.md"),
]


def selected_steps(skip_phylogeny=False, skip_structure=False):
    excluded = set()
    if skip_phylogeny:
        excluded.update(("phylogeny", "tree_support"))
    if skip_structure:
        excluded.add("structure")
    return [step for step in STEPS if step[0] not in excluded]


def execute(skip_phylogeny=False, skip_structure=False):
    results = []
    failures = set()
    for name, script, expected in selected_steps(skip_phylogeny, skip_structure):
        output = ROOT / expected
        if name == "paper_validation" and (failures or skip_phylogeny or skip_structure):
            results.append({"step": name, "command": script, "output": expected,
                            "status": "blocked", "success": False,
                            "reason": "Full upstream dataset, phylogeny and structure required"})
            print("[SKIPPED] paper_validation: incomplete upstream analysis", flush=True)
            failures.add(name)
            continue
        if name == "manuscript" and any(step in failures for step in ("dataset", "genus_audit", "extended", "phylogeny", "tree_support", "structure", "paper_validation")):
            results.append({"step": name, "command": script, "output": expected,
                            "status": "blocked", "success": False,
                            "reason": "Required analysis output missing or failed"})
            print("[SKIPPED] manuscript: an upstream analysis was incomplete", flush=True)
            failures.add(name)
            continue
        if name == "tree_support" and "phylogeny" in failures:
            print("[SKIPPED] tree_support: current phylogeny failed", flush=True)
            results.append({"step": name, "command": script, "output": expected,
                            "status": "blocked", "success": False,
                            "reason": "phylogeny failed in this run"})
            continue
        if name == "phylogeny" and not (
            shutil.which("mafft") and (shutil.which("FastTree") or shutil.which("fasttree"))
        ):
            print("[SKIPPED] phylogeny: MAFFT and FastTree must be installed and on PATH.", flush=True)
            print("Windows: install both under WSL Ubuntu and run Python inside WSL; "
                  "or install native binaries and add them to PATH.", flush=True)
            failures.add(name)
            results.append({"step": name, "command": script, "output": expected,
                            "status": "missing_dependency", "success": False,
                            "reason": "MAFFT and/or FastTree not found on PATH"})
            continue
        print(f"\\n=== {name}: {script} ===", flush=True)
        start = time.monotonic()
        # Do not accept a leftover output from an earlier successful run.
        try:
            output.unlink(missing_ok=True)
        except OSError as exc:
            raise RuntimeError(f"Cannot reset stale output {output}") from exc
        result = subprocess.run([sys.executable, str(ROOT / script)],
                                cwd=ROOT, check=False)
        ok = result.returncode == 0 and output.is_file() and output.stat().st_size > 0
        results.append({"step": name, "command": script,
                        "returncode": result.returncode, "output": expected,
                        "output_exists": output.is_file(),
                        "status": "success" if ok else "failed",
                        "success": ok, "elapsed_seconds": round(time.monotonic()-start, 2)})
        if not ok:
            failures.add(name)
            print(f"[FAILED] {name} (exit={result.returncode}), output={expected}", flush=True)
            # Independent downstream analyses must still be attempted.
    folder = ROOT / "dataset_results"
    folder.mkdir(parents=True, exist_ok=True)
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version.split()[0],
        "all_successful": not failures,
        "failed_steps": sorted(failures),
        "steps": results,
        "note": "Missing dependencies block only related stages. A partial run exits nonzero. "
                "Generated outputs do not establish causality or molecular flexibility."
    }
    (folder / "full_run_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8")
    print(f"\\nManifest: {folder / 'full_run_manifest.json'}")
    return 1 if failures else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-phylogeny", action="store_true",
                        help="Skip external MAFFT/FastTree stage and tree-support audit")
    parser.add_argument("--skip-structure", action="store_true",
                        help="Skip PDB downloads and structural distance stage")
    args = parser.parse_args()
    return execute(args.skip_phylogeny, args.skip_structure)


if __name__ == "__main__":
    raise SystemExit(main())
