# Statistical Tests Suite

A comprehensive suite for generating statistical indicators for executions of different optimization algorithms. This project provides automated analysis and comparison tools for multi-objective optimization algorithm results.

## Overview

This statistical testing suite analyzes the performance of optimization algorithms by computing various quality indicators and performing statistical tests. It's designed to work with results from multiple algorithms across different problem instances and runs.

## Requirements

- **Python 3** (for generating comparative tables)
- **C/C++ compiler** (for building statistical indicators)
- **Make** utility
- **algorithm_results directory** containing the execution results

## Directory Structure Expected

The suite expects a specific directory structure for algorithm results:

```
algorithm_results/
├── <algorithm_name>/           # e.g., COMOLSD, MOEAD, NSGA2
│   ├── <instance_name>/        # Problem instance name
│   │   ├── <run_1>/           # Individual run results
│   │   ├── <run_2>/
│   │   └── ...
│   └── ...
└── ...
```

### Example Structure
```
algorithm_results/
├── COMOLSD/
│   ├── instance_01/
│   │   ├── run_1/
│   │   ├── run_2/
│   │   └── run_3/
│   └── instance_02/
│       ├── run_1/
│       └── run_2/
├── MOEAD/
│   └── instance_01/
│       ├── run_1/
│       └── run_2/
└── NSGA2/
    └── instance_01/
        ├── run_1/
        └── run_2/
```

## Usage

### Quick Start

1. Ensure you have the `algorithm_results/` directory with your execution results
2. Run the complete analysis:

```bash
./run.sh
```

### What the Script Does

The `run.sh` script performs the following steps:

1. **Builds the project**: Compiles all C/C++ statistical indicators
2. **Generates metrics**: Runs analysis on all algorithm results
3. **Creates comparative table**: Generates `comparative_results.csv`
4. **Cleanup**: Removes temporary files

### Output

- **Log file**: `log.txt` - Contains detailed execution logs
- **Comparative table**: `comparative_results.csv` - Final comparison results
- **Analysis directory**: `analysis/` - Intermediate analysis files

## Configuration

### Instance Selection

- **Default behavior**: Processes all instances found in `algorithm_results/`
- **Custom selection**: Add instance names to `src/instances.txt` (one per line)

### Statistical Indicators

The suite includes several quality indicators:

- **Hypervolume** (`src/indicators/hypervolume/`)
- **Additive Epsilon** (`src/indicators/additive_epsilon/`)
- **Inverted Generational Distance (IGD)** (`src/indicators/igd/`)

### Statistical Tests

- **Kruskal-Wallis test** (`src/indicators/kruskal/`)
- **Mann-Whitney U test** (`src/indicators/mann_whitney/`)
- **Wilcoxon signed-rank test** (`src/indicators/wilcoxon/`)

### Utilities

- **Normalization** (`src/utils/normalize/`)
- **Filtering** (`src/utils/filter/`)
- **Boundary calculation** (`src/utils/bound/`)

## Important Configuration Notes

⚠️ **Before running the analysis, configure the parameters for:**
- **Bound calculation** (`src/utils/bound/bound_param.txt`)
- **Normalization** (`src/utils/normalize/normalize_param.txt`)
- **Hypervolume** (`src/indicators/hypervolume/hyp_ind_param_NORM.txt`)
- **Filtering** (`src/utils/filter/filter_param.txt`)

**Objective Configuration**: 
- First objective: **cost** (minimization)
- Second objective: **power** (maximization)

## Project Structure

```
├── run.sh                          # Main execution script
├── src/
│   ├── build_comparative_table.py  # Generates final comparison table
│   ├── instances.txt               # Optional: specific instances to process
│   ├── run_analysis.sh            # Core analysis script
│   ├── Makefile                    # Build configuration
│   ├── indicators/                 # Quality indicators
│   │   ├── additive_epsilon/
│   │   ├── hypervolume/
│   │   ├── igd/
│   │   ├── kruskal/
│   │   ├── mann_whitney/
│   │   └── wilcoxon/
│   └── utils/                      # Utility tools
│       ├── bound/
│       ├── conta_media_menor/
│       ├── dcdflib/
│       ├── filter/
│       └── normalize/
├── analysis/                       # Generated analysis files
├── pareto_union/                   # Pareto front unions
├── logs/                          # Execution logs
└── comparative_results.csv        # Final results table
```

## Troubleshooting

1. **Build errors**: Check that you have a C/C++ compiler and make installed
2. **Missing algorithm_results**: Ensure the directory exists and follows the expected structure
3. **Permission errors**: Make sure `run.sh` is executable (`chmod +x run.sh`)
4. **Python errors**: Verify Python 3 is installed and accessible

## Log Files

Check `log.txt` for detailed execution information and error messages. The log contains output from all build and analysis steps.

## Contributing

When adding new indicators or statistical tests, follow the existing directory structure and update the build system accordingly.