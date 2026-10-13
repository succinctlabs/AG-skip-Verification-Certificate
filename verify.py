#!/usr/bin/env python3
"""Run the paper certificate with SageMath; logs go to the ignored .audit/ directory."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--paper', type=Path, help='paper directory or explicit_instance.tex to compare with the data')
    parser.add_argument('--smoke', action='store_true', help='installation/smoke check only; NOT a certificate')
    parser.add_argument('--sage', default=os.environ.get('SAGE', 'sage'), help='SageMath executable (default: sage)')
    args = parser.parse_args()
    executable = shutil.which(args.sage)
    if executable is None:
        parser.error('SageMath was not found. Install SageMath and put sage on PATH, or pass --sage /path/to/sage.')
    root = Path(__file__).resolve().parent
    artifacts = root / '.audit'
    artifacts.mkdir(exist_ok=True)
    stem = 'smoke' if args.smoke else 'verification'
    result = artifacts / (stem + '.json')
    # Remove a previous success report before starting a new attempt.
    result.unlink(missing_ok=True)
    command = [executable, '-python', str(root / 'F2_human_audit/certify.py'), '--output', str(result)]
    if args.paper:
        command += ['--paper', str(args.paper.resolve())]
    if args.smoke:
        command += ['--smoke']
    env = dict(os.environ, PYTHONUNBUFFERED='1')
    with (artifacts / (stem + '.log')).open('w') as log:
        process = subprocess.Popen(command, cwd=root, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        try:
            for line in process.stdout:
                print(line, end='', flush=True)
                log.write(line)
                log.flush()
            code = process.wait()
        except KeyboardInterrupt:
            process.terminate()
            process.wait()
            return 130
    if code:
        print(f'FAILED (exit {code}); see {artifacts / (stem + ".log")}', file=sys.stderr)
        return code
    print(f'Report: {result}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
