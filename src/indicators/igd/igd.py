import os
import sys
import numpy as np
from pathlib import Path
from pymoo.indicators.igd import IGD

def read_reference_set(filepath):
    points = []
    with open(filepath, "r") as f:
        for line in f:
            vals = line.strip().split()
            if len(vals) == 2:
                points.append([float(vals[0]), float(vals[1])])
    return np.array(points)

def read_solutions(filepath):
    executions, current = [], []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                if current:
                    executions.append(np.array(current))
                    current = []
                continue
            vals = line.split()
            if len(vals) == 2:
                current.append([float(vals[0]), float(vals[1])])
        if current:
            executions.append(np.array(current))
    return executions

def compute_igd(reference_set, executions):
    igd = IGD(reference_set)
    return [igd(exe) for exe in executions]

def main(data_file_path: str, ref_file_path: str, output_file_path: str):
    data_file = Path(data_file_path)
    ref_file = Path(ref_file_path)
    output_file = Path(output_file_path)

    if not data_file.exists():
        print(f"[ERRO] Data file not found: {data_file}")
        sys.exit(1)
    if not ref_file.exists():
        print(f"[ERRO] Reference file not found: {ref_file}")
        sys.exit(1)

    reference_set = read_reference_set(str(ref_file))
    execs = read_solutions(str(data_file))
    igd_values = compute_igd(reference_set, execs)

    with open(output_file, "w") as f:
        for v in igd_values:
            f.write(f"{v}\n")
        f.write("\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Uso: python3 igd.py <data_file> <reference_set> <output_file>")
        sys.exit(1)
    
    data_file_arg = sys.argv[1]
    ref_file_arg = sys.argv[2]
    output_file_arg = sys.argv[3]
    
    main(data_file_arg, ref_file_arg, output_file_arg)