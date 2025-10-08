import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
import re

# O script espera o caminho do diretório 'analysis' como argumento
if len(sys.argv) < 2:
    print("Erro: O caminho para o diretório de análise não foi fornecido.")
    sys.exit(1)

base_dir = Path(sys.argv[1])
results = []

def find_file(instance_path, sub_dir, expected_name):
    search_path = os.path.join(instance_path, sub_dir)
    if not os.path.exists(search_path):
        return None
    for f in os.listdir(search_path):
        if expected_name.lower() in f.strip().lower():
            return os.path.join(search_path, f)
    return None

def parse_kruskal_wallis(file_path):
    """
    Lê todas as linhas do arquivo Kruskal-Wallis e retorna apenas aquelas
    com p-value <= 0.05. Se nenhuma linha atender ao critério, retorna 'H0'.
    """
    if not file_path or not os.path.exists(file_path):
        return "N/A"

    significant_lines = []
    pattern = re.compile(r"p-value of ([0-9.eE+-]+)")

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            match = pattern.search(line)
            if match:
                try:
                    p_value = float(match.group(1))
                    if p_value <= 0.05:
                        significant_lines.append(line)
                except ValueError:
                    continue

    if significant_lines:
        return " | ".join(significant_lines)
    else:
        return "H0"

def read_values_and_calc_mean(file_path):
    """Lê os valores de um arquivo e calcula a média."""
    if not file_path or not os.path.exists(file_path):
        return np.nan
    with open(file_path, "r") as f:
        values = [float(line.strip()) for line in f if line.strip()]
    return np.mean(values) if values else np.nan

# Lendo a lista de instâncias processadas
processed_instances_file = base_dir / "processed_instances.txt"
if not processed_instances_file.exists():
    print("Erro: O arquivo de instâncias processadas não foi encontrado.")
    sys.exit(1)

with open(processed_instances_file, 'r') as f:
    instances = [line.strip() for line in f if line.strip()]

# Itera apenas sobre as instâncias da lista
for instance in instances:
    instance_path = os.path.join(base_dir, instance)
    if not os.path.isdir(instance_path):
        print(f"Aviso: Diretório da instância {instance} não encontrado.")
        continue

    try:
        # Encontra todos os arquivos de métricas
        hv_moead_file = find_file(instance_path, "hypervolume", "HV_moead")
        hv_comolsd_file = find_file(instance_path, "hypervolume", "HV_comolsd")
        hv_nsga2_file = find_file(instance_path, "hypervolume", "HV_nsga2")
        
        eps_moead_file = find_file(instance_path, "epsilon_additive", "esp_ad_moead")
        eps_comolsd_file = find_file(instance_path, "epsilon_additive", "esp_ad_comolsd")
        eps_nsga2_file = find_file(instance_path, "epsilon_additive", "esp_ad_nsga2")
        
        igd_moead_file = find_file(instance_path, "igd", "IGD_moead") 
        igd_comolsd_file = find_file(instance_path, "igd", "IGD_comolsd")
        igd_nsga2_file = find_file(instance_path, "igd", "IGD_nsga2")
        
        # Encontra os arquivos do Kruskal-Wallis
        kruskal_hv_file = find_file(instance_path, "kruskal", "hv_saidakruskal")
        kruskal_eps_file = find_file(instance_path, "kruskal", "eps_saidakruskal")
        kruskal_igd_file = find_file(instance_path, "kruskal", "igd_saidakruskal")

        # Calcula as médias
        hv_moead_mean = read_values_and_calc_mean(hv_moead_file)
        hv_comolsd_mean = read_values_and_calc_mean(hv_comolsd_file)
        hv_nsga2_mean = read_values_and_calc_mean(hv_nsga2_file)
        
        eps_moead_mean = read_values_and_calc_mean(eps_moead_file)
        eps_comolsd_mean = read_values_and_calc_mean(eps_comolsd_file)
        eps_nsga2_mean = read_values_and_calc_mean(eps_nsga2_file)
        
        igd_moead_mean = read_values_and_calc_mean(igd_moead_file) 
        igd_comolsd_mean = read_values_and_calc_mean(igd_comolsd_file)
        igd_nsga2_mean = read_values_and_calc_mean(igd_nsga2_file)

        # Processa resultados do Kruskal-Wallis separadamente
        kruskal_hv_result = parse_kruskal_wallis(kruskal_hv_file)
        kruskal_eps_result = parse_kruskal_wallis(kruskal_eps_file)
        kruskal_igd_result = parse_kruskal_wallis(kruskal_igd_file)

        results.append({
            "Instance": instance,
            "HV_MOEA_D": hv_moead_mean,
            "HV_COMOLS_D": hv_comolsd_mean,
            "HV_NSGA2": hv_nsga2_mean,
            "EPS_MOEA_D": eps_moead_mean,
            "EPS_COMOLS_D": eps_comolsd_mean,
            "EPS_NSGA2": eps_nsga2_mean,
            "IGD_MOAE_D": igd_moead_mean,
            "IGD_COMOLS_D": igd_comolsd_mean,
            "IGD_NSGA2": igd_nsga2_mean,
            "Kruskal Wallis Test (HV)": kruskal_hv_result,
            "Kruskal Wallis Test (EPS)": kruskal_eps_result,
            "Kruskal Wallis Test (IGD)": kruskal_igd_result
        })

    except Exception as e:
        print(f"Erro ao processar a instância {instance}: {e}")

if results:
    column_order = [
        "Instance",
        "HV_MOEA_D", "HV_COMOLS_D", "HV_NSGA2",
        "EPS_MOEA_D", "EPS_COMOLS_D", "EPS_NSGA2",
        "IGD_MOAE_D", "IGD_COMOLS_D", "IGD_NSGA2",
        "Kruskal Wallis Test (HV)", "Kruskal Wallis Test (EPS)", "Kruskal Wallis Test (IGD)"
    ]
    
    df = pd.DataFrame(results, columns=column_order)
    for col in df.columns:
        if "Instance" not in col and "Kruskal" not in col:
            df[col] = df[col].map('{:.4f}'.format)
            
    df.to_csv(base_dir.parent / "comparative_results.csv", index=False)
else:
    print("Nenhum resultado processado.")

