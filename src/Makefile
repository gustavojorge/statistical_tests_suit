CC=gcc
CXX=g++ 
CFLAGS=-Wall -pedantic -O2
LDFLAGS=-lm

BIN_DIR=bin
UTILS_DIR=utils
INDICATORS_DIR=indicators

DCDFLIB_OBJ=$(UTILS_DIR)/dcdflib/dcdflib.o

UTILS_EXEC=$(BIN_DIR)/bound $(BIN_DIR)/normalize $(BIN_DIR)/filter
IND_EXEC=$(BIN_DIR)/eps_ind $(BIN_DIR)/hyp_ind $(BIN_DIR)/mann-whit $(BIN_DIR)/kruskal-wallis $(BIN_DIR)/wilcoxon-sign
EXECUTABLES=$(UTILS_EXEC) $(IND_EXEC)

all: $(BIN_DIR) $(EXECUTABLES)
	@echo "Compilation complete."
	@echo " "

$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

#########################
# Utils
#########################

	@echo "[BUILDING FILES]"

$(BIN_DIR)/bound: $(UTILS_DIR)/bound/bound.cc
	@echo "--> Compiling bound"
	@$(CXX) $(CFLAGS) $< -o $@ $(LDFLAGS) >/dev/null 2>&1

$(BIN_DIR)/normalize: $(UTILS_DIR)/normalize/normalize.cc
	@echo "--> Compiling normalize"
	@$(CXX) $(CFLAGS) $< -o $@ $(LDFLAGS) >/dev/null 2>&1

$(BIN_DIR)/filter: $(UTILS_DIR)/filter/filter.cc
	@echo "--> Compiling filter"
	@$(CXX) $(CFLAGS) $< -o $@ $(LDFLAGS) >/dev/null 2>&1

#########################
# Indicators
#########################

$(BIN_DIR)/eps_ind: $(INDICATORS_DIR)/additive_epsilon/eps_ind.c
	@echo "--> Compiling eps_ind"
	@$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) >/dev/null 2>&1

$(BIN_DIR)/hyp_ind: $(INDICATORS_DIR)/hypervolume/hyp_ind.c
	@echo "--> Compiling hyp_ind"
	@$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) >/dev/null 2>&1

$(BIN_DIR)/mann-whit: $(INDICATORS_DIR)/mann_whitney/mann-whit.cc $(DCDFLIB_OBJ)
	@echo "--> Compiling mann-whit"
	@$(CXX) $(CFLAGS) -I$(UTILS_DIR)/dcdflib $^ -o $@ $(LDFLAGS) >/dev/null 2>&1

$(BIN_DIR)/kruskal-wallis: $(INDICATORS_DIR)/kruskal/kruskal-wallis.cc $(DCDFLIB_OBJ)
	@echo "--> Compiling kruskal-wallis"
	@$(CXX) $(CFLAGS) -I$(UTILS_DIR)/dcdflib $^ -o $@ $(LDFLAGS) >/dev/null 2>&1

$(BIN_DIR)/wilcoxon-sign: $(INDICATORS_DIR)/wilcoxon/wilcoxon-sign.cc $(DCDFLIB_OBJ)
	@echo "--> Compiling wilcoxon-sign"
	@$(CXX) $(CFLAGS) -I$(UTILS_DIR)/dcdflib $^ -o $@ $(LDFLAGS) >/dev/null 2>&1

#########################
# Dcdflib
#########################
$(DCDFLIB_OBJ): $(UTILS_DIR)/dcdflib/dcdflib.cc
	@echo "--> Compiling dcdflib"
	$(CXX) $(CFLAGS) -c $< -o $@ >/dev/null 2>&1

clean:
	@rm -rf ../ParetoUnion/ >/dev/null 2>&1
	@rm -rf ../logs/ >/dev/null 2>&1
	@rm -rf ../analysis/ >/dev/null 2>&1
	@rm -rf ../comparative_results.csv >/dev/null 2>&1
	@rm -rf ../log.txt >/dev/null 2>&1
	@rm -rf $(BIN_DIR) >/dev/null 2>&1
	@rm -f $(UTILS_DIR)/*/*.o $(INDICATORS_DIR)/*/*.o >/dev/null 2>&1
	@rm -f *~ >/dev/null 2>&1
	@echo "Cleaning completed."
