"""Output paths for reproducible runs without changing recorded evidence."""
from datetime import datetime, timezone
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def output_directory(path):
    path = Path(path).resolve()
    if path.is_relative_to(ROOT) and not path.is_relative_to(ROOT / 'outputs'):
        raise ValueError('Write reruns under outputs/ or outside the repository')
    path.mkdir(parents=True, exist_ok=True)
    return path


def output_file(path):
    path = Path(path).resolve()
    output_directory(path.parent)
    return path


def report_output(filename):
    parser = argparse.ArgumentParser(description=__import__('__main__').__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'outputs/validation' / filename)
    args = parser.parse_args()
    return output_file(args.output)


def save_report(report, path):
    path.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
