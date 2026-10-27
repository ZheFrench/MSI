#!/usr/bin/env python3
"""
make_msi_locker.py

Fills the MSI locker from ADIVAR's: every run ADIVAR has finished (status 3 in
--locker-adivar) is added to --locker-msi (one run name per line, created if
missing; the runs already there are kept). msi_detection_pre_p.py (cron mode)
only takes the runs finished by ADIVAR and not in the MSI locker: after this
script, it never takes these runs, only the ones ADIVAR finishes later. Used
to start the cron without it analysing every past run.

Same rules as the cron (functions of msi_detection_pre_p.py): unfinished runs
(status other than 3) are left out, so the cron still takes them when ADIVAR
finishes them.

    make_msi_locker.py --locker-adivar adivar.json --locker-msi msi-sensor.locker.txt [--dry-run]

--dry-run: prints the runs that would be added, writes nothing.
"""

import argparse
import logging
import sys
from pathlib import Path

# msi_detection_pre_p.py sits in main-scripts/, next to this script's folder.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "main-scripts"))
from msi_detection_pre_p import list_runs_done, list_runs_msi, lock  # noqa: E402

logger = logging.getLogger("make-msi-locker")


def main():
    parser = argparse.ArgumentParser(description="Adds every run finished by ADIVAR to the MSI locker.")
    parser.add_argument("--locker-adivar", type=Path, required=True, help="ADIVAR's locker JSON (read only).")
    parser.add_argument("--locker-msi", type=Path, required=True,
                        help="MSI locker, one run per line (created if missing).")
    parser.add_argument("--dry-run", action="store_true", help="Print the runs to add, write nothing.")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(name)s: %(message)s")

    if not args.locker_adivar.is_file():
        sys.exit(f"no ADIVAR locker: {args.locker_adivar}")
    done = list_runs_done(args.locker_adivar)
    already = set(list_runs_msi(args.locker_msi))
    to_add = [run for run in done if run not in already]
    logger.info(f"{len(done)} run(s) finished by ADIVAR, {len(done) - len(to_add)} already in the MSI locker, "
                f"{len(to_add)} to add")

    for run in to_add:
        if args.dry_run:
            print(run)
        else:
            lock(args.locker_msi, run)


if __name__ == "__main__":
    main()
