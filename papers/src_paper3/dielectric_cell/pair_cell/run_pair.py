"""Physics-run entry point. Executes the checkpointed run plan:
validation gates (validate2d), the coarse-grid phase (phase1), the
symmetry-reduced fine-grid ladder (phase2), then assembles results.md
(make_results). All solves append to checkpoint.jsonl on completion and
are skipped on re-run, so this script is safely resumable."""
import validate2d
import phase1
import phase2
import make_results

if __name__ == "__main__":
    validate2d.main()
    phase1.main()
    phase2.main()
    make_results.main()
