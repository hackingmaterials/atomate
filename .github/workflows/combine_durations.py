import json
from pathlib import Path

split_prefix = "split-"
durations_path = Path(".test_durations")

split_paths = Path(".").glob(f"{split_prefix}*/{durations_path.name}")

try:
    previous_durations = json.loads(durations_path.read_text())
except FileNotFoundError:
    previous_durations = {}

new_durations = previous_durations.copy()

for path in split_paths:
    durations = json.loads(path.read_text())
    new_durations.update(durations)

durations_path.parent.mkdir(parents=True, exist_ok=True)
durations_path.write_text(json.dumps(new_durations))
