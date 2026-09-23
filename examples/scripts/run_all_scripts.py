#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Run all FinancePy example scripts.")
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Enable interactive plotting.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=120,
        help="Maximum runtime per example in seconds (default: 120).",
    )

    parser.add_argument(
        "--fail-fast",
        action="store_true",
        help="Stop after the first failed or timed-out example.",
    )

    args = parser.parse_args()

    # This file is:
    #
    #   examples/scripts/run_all_scripts.py
    #
    # Therefore:
    #
    #   scripts_dir  = examples/scripts
    #   examples_dir = examples
    #
    scripts_dir = Path(__file__).resolve().parent
    examples_dir = scripts_dir.parent

    # Find all example*.py files below examples/scripts.
    scripts = sorted(scripts_dir.rglob("example_*.py"))

    if not scripts:
        print(f"No example*.py files found in {scripts_dir}")
        return 1

    print(f"Examples directory : {examples_dir}")
    print(f"Scripts directory  : {scripts_dir}")
    print(f"Python             : {sys.executable}")
    print(f"Plotting           : {'enabled' if args.plot else 'disabled'}")
    print(f"Timeout            : {args.timeout} seconds")
    print(f"Examples found     : {len(scripts)}")
    print()

    passed = []
    failed = []
    timed_out = []

    start_all = time.perf_counter()

    for i, script in enumerate(scripts, 1):
        relative = script.relative_to(examples_dir)

        # Convert:
        #
        #   scripts/bonds/example_bond_frn.py
        #
        # into:
        #
        #   scripts.bonds.example_bond_frn
        #
        module = ".".join(relative.with_suffix("").parts)

        print("=" * 80)
        print(f"[{i}/{len(scripts)}] {relative}")
        print(f"Module: {module}")
        print("=" * 80)

        env = os.environ.copy()

        # Prevent matplotlib windows from opening unless --plot is supplied.
        if not args.plot:
            env["MPLBACKEND"] = "Agg"

        start = time.perf_counter()

        try:
            result = subprocess.run(
                [sys.executable, "-m", module],
                cwd=examples_dir,
                env=env,
                timeout=args.timeout,
                capture_output=True,
                text=True,
            )

            elapsed = time.perf_counter() - start

            if result.returncode == 0:
                passed.append((relative, elapsed))
                print(f"\nPASS ({elapsed:.2f}s): {relative}")

            else:
                failed.append((relative, result.returncode, elapsed))

                print(
                    f"\nFAIL ({elapsed:.2f}s, "
                    f"exit={result.returncode}): {relative}"
                )

                if result.stdout:
                    print("\nSTDOUT:")
                    print(result.stdout)

                if result.stderr:
                    print("\nSTDERR:")
                    print(result.stderr)

                if args.fail_fast:
                    break

        except subprocess.TimeoutExpired:
            elapsed = time.perf_counter() - start
            timed_out.append((relative, elapsed))

            print(f"\nTIMEOUT ({elapsed:.2f}s): {relative}")

        except KeyboardInterrupt:
            print("\n\nInterrupted by user.")
            return 130

        print()

    total_time = time.perf_counter() - start_all

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)

    print(f"Total       : {len(scripts)}")
    print(f"Passed      : {len(passed)}")
    print(f"Failed      : {len(failed)}")
    print(f"Timed out   : {len(timed_out)}")
    print(f"Total time  : {total_time:.2f}s")

    if failed:
        print("\nFailed examples:")

        for script, returncode, elapsed in failed:
            print(f"  {script} " f"(exit={returncode}, {elapsed:.2f}s)")

    if timed_out:
        print("\nTimed-out examples:")

        for script, elapsed in timed_out:
            print(f"  {script} " f"({elapsed:.2f}s)")

    print()

    if failed or timed_out:
        return 1

    print("All examples passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
