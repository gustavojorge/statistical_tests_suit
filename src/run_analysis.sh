#!/bin/bash
# Ordem:
# MOEAD
# COMOLSD
# NSGA2

##### NÃO ESQUECER DE CONFIGURAR OS PARAMETROS DO bound, normalize, HV, filter
## Primeiro objetivo é o custo (min), e o segundo objetivo é a potência (max)

ROOT_DIR=$(dirname "$PWD")
PYTHON_CMD="python3"

# Define o caminho para o arquivo de instâncias
INSTANCES_FILE="$ROOT_DIR/src/instances.txt"

# Cria a lista de instâncias a serem processadas
INSTANCES=()
if [ -s "$INSTANCES_FILE" ]; then
  echo "[READING INSTANCES]"
  while IFS= read -r line; do
    line="${line%$'\r'}"
    if [ -n "$line" ]; then
      INSTANCES+=("$line")
    fi
  done < "$INSTANCES_FILE"
else
  echo "[Empty instance file. Processing all instances]."
  for algorithm_dir in "$ROOT_DIR"/algorithm_results/*; do
    if [ -d "$algorithm_dir" ]; then
      for instance_dir in "$algorithm_dir"/*; do
        if [ -d "$instance_dir" ]; then
          instance_name=$(basename "$instance_dir")
          if [[ ! " ${INSTANCES[*]} " =~ " ${instance_name} " ]]; then
            INSTANCES+=("$instance_name")
          fi
        fi
      done
    fi
  done
fi

mkdir -p "$ROOT_DIR"/pareto_union/MOEAD
mkdir -p "$ROOT_DIR"/pareto_union/COMOLSD
mkdir -p "$ROOT_DIR"/pareto_union/NSGA2
mkdir -p "$ROOT_DIR"/logs
mkdir -p "$ROOT_DIR"/analysis

for p in "${INSTANCES[@]}"
do
  echo "--- [Processing instance: $p]"
  mkdir -p "$ROOT_DIR"/analysis/"$p"/utils
  mkdir -p "$ROOT_DIR"/analysis/"$p"/epsilon_additive
  mkdir -p "$ROOT_DIR"/analysis/"$p"/hypervolume
  mkdir -p "$ROOT_DIR"/analysis/"$p"/igd
  mkdir -p "$ROOT_DIR"/analysis/"$p"/mann_whitney
  mkdir -p "$ROOT_DIR"/analysis/"$p"/kruskal

  # Une todas as execuções da instância p em um único arquivo
  for i in {1..20}
  do
    cat "$ROOT_DIR"/algorithm_results/MOEAD/"$p"/"$i"/"$p"_moead_1000000.txt >> "$ROOT_DIR"/pareto_union/MOEAD/"$p"_union_pareto_file.out
    cat "$ROOT_DIR"/algorithm_results/COMOLSD/"$p"/"$i"/"$p"_comolsd_1000000.txt >> "$ROOT_DIR"/pareto_union/COMOLSD/"$p"_union_pareto_file.out
    cat "$ROOT_DIR"/algorithm_results/NSGA2/"$p"/"$i"/"$p"_nsga2_1000000.txt >> "$ROOT_DIR"/pareto_union/NSGA2/"$p"_union_pareto_file.out
    echo "" >> "$ROOT_DIR"/pareto_union/MOEAD/"$p"_union_pareto_file.out
    echo "" >> "$ROOT_DIR"/pareto_union/COMOLSD/"$p"_union_pareto_file.out
    echo "" >> "$ROOT_DIR"/pareto_union/NSGA2/"$p"_union_pareto_file.out
  done

  # Cria inputBound.in juntando os três arquivos de união
  cat "$ROOT_DIR"/pareto_union/MOEAD/"$p"_union_pareto_file.out \
      "$ROOT_DIR"/pareto_union/COMOLSD/"$p"_union_pareto_file.out \
      "$ROOT_DIR"/pareto_union/NSGA2/"$p"_union_pareto_file.out \
      > "$ROOT_DIR"/analysis/"$p"/utils/inputBound.in

  # Bound
  echo "- [RUNNING] bound for instance $p"
  "$ROOT_DIR"/src/bin/bound "$ROOT_DIR"/src/utils/bound/bound_param.txt \
    "$ROOT_DIR"/analysis/"$p"/utils/inputBound.in \
    "$ROOT_DIR"/analysis/"$p"/utils/bound.out

  # Normalize
  echo "- [RUNNING] normalize for instance $p"
  "$ROOT_DIR"/src/bin/normalize "$ROOT_DIR"/src/utils/normalize/normalize_param.txt "$ROOT_DIR"/analysis/"$p"/utils/bound.out "$ROOT_DIR"/pareto_union/MOEAD/"$p"_union_pareto_file.out "$ROOT_DIR"/analysis/"$p"/utils/moead_normalizado.out
  "$ROOT_DIR"/src/bin/normalize "$ROOT_DIR"/src/utils/normalize/normalize_param.txt "$ROOT_DIR"/analysis/"$p"/utils/bound.out "$ROOT_DIR"/pareto_union/COMOLSD/"$p"_union_pareto_file.out "$ROOT_DIR"/analysis/"$p"/utils/comolsd_normalizado.out
  "$ROOT_DIR"/src/bin/normalize "$ROOT_DIR"/src/utils/normalize/normalize_param.txt "$ROOT_DIR"/analysis/"$p"/utils/bound.out "$ROOT_DIR"/pareto_union/NSGA2/"$p"_union_pareto_file.out "$ROOT_DIR"/analysis/"$p"/utils/nsga2_normalizado.out

  # Filter
  echo "- [RUNNING] filter for instance $p"
  cat "$ROOT_DIR"/analysis/"$p"/utils/moead_normalizado.out "$ROOT_DIR"/analysis/"$p"/utils/comolsd_normalizado.out "$ROOT_DIR"/analysis/"$p"/utils/nsga2_normalizado.out >> "$ROOT_DIR"/analysis/"$p"/utils/inputFilter.in
  "$ROOT_DIR"/src/bin/filter "$ROOT_DIR"/src/utils/filter/filter_param.txt "$ROOT_DIR"/analysis/"$p"/utils/inputFilter.in "$ROOT_DIR"/analysis/"$p"/reference_set.out

  # Hypervolume
  echo "- [RUNNING] hyp_ind for instance $p"
  "$ROOT_DIR"/src/bin/hyp_ind "$ROOT_DIR"/src/indicators/hypervolume/hyp_ind_param_NORM.txt "$ROOT_DIR"/analysis/"$p"/utils/moead_normalizado.out "$ROOT_DIR"/analysis/"$p"/reference_set.out "$ROOT_DIR"/analysis/"$p"/hypervolume/HV_moead.out
  "$ROOT_DIR"/src/bin/hyp_ind "$ROOT_DIR"/src/indicators/hypervolume/hyp_ind_param_NORM.txt "$ROOT_DIR"/analysis/"$p"/utils/comolsd_normalizado.out "$ROOT_DIR"/analysis/"$p"/reference_set.out "$ROOT_DIR"/analysis/"$p"/hypervolume/HV_comolsd.out
  "$ROOT_DIR"/src/bin/hyp_ind "$ROOT_DIR"/src/indicators/hypervolume/hyp_ind_param_NORM.txt "$ROOT_DIR"/analysis/"$p"/utils/nsga2_normalizado.out "$ROOT_DIR"/analysis/"$p"/reference_set.out "$ROOT_DIR"/analysis/"$p"/hypervolume/HV_nsga2.out

  echo "" >> "$ROOT_DIR"/analysis/"$p"/hypervolume/HV_moead.out
  echo "" >> "$ROOT_DIR"/analysis/"$p"/hypervolume/HV_comolsd.out
  echo "" >> "$ROOT_DIR"/analysis/"$p"/hypervolume/HV_nsga2.out

  # Epsilon
  echo "- [RUNNING] eps_ind for instance $p"
  "$ROOT_DIR"/src/bin/eps_ind "$ROOT_DIR"/src/indicators/additive_epsilon/eps_ind_param.txt "$ROOT_DIR"/analysis/"$p"/utils/moead_normalizado.out "$ROOT_DIR"/analysis/"$p"/reference_set.out "$ROOT_DIR"/analysis/"$p"/epsilon_additive/esp_ad_moead.out
  "$ROOT_DIR"/src/bin/eps_ind "$ROOT_DIR"/src/indicators/additive_epsilon/eps_ind_param.txt "$ROOT_DIR"/analysis/"$p"/utils/comolsd_normalizado.out "$ROOT_DIR"/analysis/"$p"/reference_set.out "$ROOT_DIR"/analysis/"$p"/epsilon_additive/esp_ad_comolsd.out
  "$ROOT_DIR"/src/bin/eps_ind "$ROOT_DIR"/src/indicators/additive_epsilon/eps_ind_param.txt "$ROOT_DIR"/analysis/"$p"/utils/nsga2_normalizado.out "$ROOT_DIR"/analysis/"$p"/reference_set.out "$ROOT_DIR"/analysis/"$p"/epsilon_additive/esp_ad_nsga2.out

  echo "" >> "$ROOT_DIR"/analysis/"$p"/epsilon_additive/esp_ad_moead.out
  echo "" >> "$ROOT_DIR"/analysis/"$p"/epsilon_additive/esp_ad_comolsd.out
  echo "" >> "$ROOT_DIR"/analysis/"$p"/epsilon_additive/esp_ad_nsga2.out

  # IGD
  echo "- [RUNNING] igd.py for instance $p"
  "$PYTHON_CMD" "$ROOT_DIR"/src/indicators/igd/igd.py "$ROOT_DIR"/analysis/"$p"/utils/moead_normalizado.out "$ROOT_DIR"/analysis/"$p"/reference_set.out "$ROOT_DIR"/analysis/"$p"/igd/IGD_moead.out
  "$PYTHON_CMD" "$ROOT_DIR"/src/indicators/igd/igd.py "$ROOT_DIR"/analysis/"$p"/utils/comolsd_normalizado.out "$ROOT_DIR"/analysis/"$p"/reference_set.out "$ROOT_DIR"/analysis/"$p"/igd/IGD_comolsd.out
  "$PYTHON_CMD" "$ROOT_DIR"/src/indicators/igd/igd.py "$ROOT_DIR"/analysis/"$p"/utils/nsga2_normalizado.out "$ROOT_DIR"/analysis/"$p"/reference_set.out "$ROOT_DIR"/analysis/"$p"/igd/IGD_nsga2.out

  # Kruskal-Wallis
  echo "- [RUNNING] kruskal-wallis for instance $p"
  cat "$ROOT_DIR"/analysis/"$p"/hypervolume/HV_moead.out "$ROOT_DIR"/analysis/"$p"/hypervolume/HV_comolsd.out "$ROOT_DIR"/analysis/"$p"/hypervolume/HV_nsga2.out >> "$ROOT_DIR"/analysis/"$p"/kruskal/hv_kruskal.in
  "$ROOT_DIR"/src/bin/kruskal-wallis "$ROOT_DIR"/analysis/"$p"/kruskal/hv_kruskal.in "$ROOT_DIR"/src/indicators/kruskal/kruskalparam.txt "$ROOT_DIR"/analysis/"$p"/kruskal/hv_saidakruskal.out >> "$ROOT_DIR"/logs/log_hv_kruskal.txt 2>&1
  echo "" >> "$ROOT_DIR"/logs/log_hv_kruskal.txt

  cat "$ROOT_DIR"/analysis/"$p"/epsilon_additive/esp_ad_moead.out "$ROOT_DIR"/analysis/"$p"/epsilon_additive/esp_ad_comolsd.out "$ROOT_DIR"/analysis/"$p"/epsilon_additive/esp_ad_nsga2.out >> "$ROOT_DIR"/analysis/"$p"/kruskal/eps_kruskal.in
  "$ROOT_DIR"/src/bin/kruskal-wallis "$ROOT_DIR"/analysis/"$p"/kruskal/eps_kruskal.in "$ROOT_DIR"/src/indicators/kruskal/kruskalparam.txt "$ROOT_DIR"/analysis/"$p"/kruskal/eps_saidakruskal.out >> "$ROOT_DIR"/logs/log_eps_kruskal.txt 2>&1
  echo "" >> "$ROOT_DIR"/logs/log_eps_kruskal.txt

  cat "$ROOT_DIR"/analysis/"$p"/igd/IGD_moead.out "$ROOT_DIR"/analysis/"$p"/igd/IGD_comolsd.out "$ROOT_DIR"/analysis/"$p"/igd/IGD_nsga2.out >> "$ROOT_DIR"/analysis/"$p"/kruskal/igd_kruskal.in
  "$ROOT_DIR"/src/bin/kruskal-wallis "$ROOT_DIR"/analysis/"$p"/kruskal/igd_kruskal.in "$ROOT_DIR"/src/indicators/kruskal/kruskalparam.txt "$ROOT_DIR"/analysis/"$p"/kruskal/igd_saidakruskal.out >> "$ROOT_DIR"/logs/log_igd_kruskal.txt 2>&1
  echo "" >> "$ROOT_DIR"/logs/log_igd_kruskal.txt
done

# Salva a lista de instâncias processadas em um arquivo temporário
INSTANCES_LIST_FILE="$ROOT_DIR/analysis/processed_instances.txt"
printf "%s\n" "${INSTANCES[@]}" > "$INSTANCES_LIST_FILE"
echo " "
echo "--- List of processed instances saved in $INSTANCES_LIST_FILE"

