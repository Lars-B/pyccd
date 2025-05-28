import subprocess
from pathlib import Path


def generate_transcope_help():
    output_file = Path(__file__).parent / "transcope_help.txt"
    with open(output_file, "w") as f:
        subprocess.run(['transcope', '--help'], stdout=f, check=True)
