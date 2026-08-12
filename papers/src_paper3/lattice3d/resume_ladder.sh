#!/bin/bash
# Cron-driven auto-resume for the lattice3d ladder. Fires no LLM
# requests: relaunches the checkpointed driver if it is not running,
# and disables itself (via the LADDER_DONE marker) once the full
# ladder (19 physics cases) is checkpointed and the driver exits
# cleanly. Safe against double-launch via flock.
DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR" || exit 1
[ -f LADDER_DONE ] && exit 0
exec 9>"$DIR/.resume.lock"
flock -n 9 || exit 0
if pgrep -f "python3.11 driver.py" >/dev/null; then exit 0; fi
# already complete?
N=$(grep -c '"grid"' checkpoint.jsonl 2>/dev/null || echo 0)
if [ "$N" -ge 19 ]; then touch LADDER_DONE; exit 0; fi
echo "[resume_ladder] relaunch at $(date), $N/19 cases done" >> driver.log
nice -n 10 python3.11 driver.py >> driver.log 2>&1
if [ $? -eq 0 ]; then
  N=$(grep -c '"grid"' checkpoint.jsonl 2>/dev/null || echo 0)
  [ "$N" -ge 19 ] && touch LADDER_DONE
fi
