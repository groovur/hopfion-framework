import json
import os
import time

CHECKPOINT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "checkpoint.jsonl")


def append(record):
    """Append one result record to checkpoint.jsonl, flushed to disk
    immediately so no completed solve is ever lost."""
    record = dict(record)
    record["timestamp"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    with open(CHECKPOINT_PATH, "a") as f:
        f.write(json.dumps(record) + "\n")
        f.flush()
        os.fsync(f.fileno())


def load():
    """Return all checkpointed records as a list of dicts."""
    if not os.path.exists(CHECKPOINT_PATH):
        return []
    out = []
    with open(CHECKPOINT_PATH) as f:
        for line in f:
            line = line.strip()
            if line:
                out.append(json.loads(line))
    return out


def have(records, **match):
    """True if a record matching all key=value pairs already exists."""
    for r in records:
        if all(r.get(k) == v for k, v in match.items()):
            return True
    return False
