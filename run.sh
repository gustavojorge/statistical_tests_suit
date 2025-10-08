#!/bin/bash

BLUE="\033[0;34m"
GREEN="\033[0;32m"
YELLOW="\033[0;33m"
RED="\033[0;31m"
PURPLE="\033[0;35m"
NC="\033[0m" 

PROJECT_NAME="==================== Statistical Tests Suit ===================="

# Define o diretório raiz (onde está este script)
ROOT_DIR=$(pwd)
SRC_DIR="$ROOT_DIR/src"
ANALYSIS_DIR="$ROOT_DIR/analysis"
TABLE_SCRIPT="$SRC_DIR/build_comparative_table.py"
PROCESSED_INSTANCES_FILE="$ANALYSIS_DIR/processed_instances.txt"
PYTHON_CMD="/usr/bin/python3"

# Corrigido: criar diretório de logs no local certo
mkdir -p "$ROOT_DIR/logs"
LOG_FILE="$ROOT_DIR/log.txt"

echo -e "${YELLOW}--- ${PROJECT_NAME} ---${NC}"
echo -e ""

echo -e "${YELLOW}
#############################################################################
# DON'T FORGET TO SET THE PARAMETERS OF BOUND, NORMALIZE, HV, FILTER        #
# The first objective is cost (min), and the second objective is power (max)#
#############################################################################"
echo -e ""


# Redirect all output from run_analysis.sh to the log file
touch "$LOG_FILE"

# --- Step 1: Clean and build the project ---
echo -e "${BLUE}[BUILDING PROJECT]${NC}"
pushd "$SRC_DIR" > /dev/null
make clean >> "$LOG_FILE" 2>&1
make >> "$LOG_FILE" 2>&1
popd > /dev/null
echo -e "${GREEN}  -> Build completed.${NC}"

# --- Step 2: Run the analysis script ---
echo -e "${BLUE}[GENERATING METRICS]${NC}"
pushd "$SRC_DIR" > /dev/null
./run_analysis.sh >> "$LOG_FILE" 2>&1
popd > /dev/null
echo -e "${GREEN}  -> Generating metrics completed.${NC}"

# --- Step 3: Generate the comparison table ---
echo -e "${BLUE}[GENERATING COMPARATIVE TABLE]${NC}"
"$PYTHON_CMD" "$TABLE_SCRIPT" "$ANALYSIS_DIR" >> "$LOG_FILE" 2>&1
echo -e "${GREEN}  -> Comparative table generated in $ROOT_DIR/comparative_results.csv${NC}"

# --- Step 4: Clean up temporary files ---
echo -e "${BLUE}[CLEANING TEMPORARY FILES]${NC}"
rm -f "$PROCESSED_INSTANCES_FILE"
echo -e "${GREEN}  -> File $PROCESSED_INSTANCES_FILE removed.${NC}"

echo "---"
echo -e "${YELLOW}Process completed successfully. Check 'log.txt' for details.${NC}"
