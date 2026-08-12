#!/bin/bash
# Overnight coherent-background ladder. Strictly sequential (single
# process); every solve checkpoints to checkpoint.jsonl on completion,
# so the ladder can be interrupted and rerun without losing work.
# Quadrant solver (validated against full-domain at h=0.2 to 0.02%).
cd "$(dirname "$0")"
LOG=ladder.log
{
  echo "=== ladder start $(date) ==="
  nice -n 10 python3.11 coherent_pair.py --quadrant \
      --h 0.1 0.05 0.025 --eta 1e-2 3e-3 1e-3 --two-channel
  echo "=== gradf ladder done $(date) ==="
  nice -n 10 python3.11 coherent_pair.py --quadrant --prescription pullback \
      --h 0.1 0.05 --eta 1e-2 3e-3
  echo "=== pullback reference done $(date) ==="
} >> "$LOG" 2>&1
