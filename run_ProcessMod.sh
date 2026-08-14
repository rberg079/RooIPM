#!/bin/bash
#SBATCH --job-name=IPM_muFun             # Nom de la tâche
#SBATCH --output=job_output_muFun.txt    # Fichier de sortie (unique par tâche)
#SBATCH --error=job_error_muFun.txt      # Fichier d'erreur (unique par tâche)
#SBATCH --nodes=1                       # Nombre de nœuds à utiliser
#SBATCH --cpus-per-task=4               # Nombre de cœurs par tâche
#SBATCH --mem=64G                       # Mémoire allouée 
#SBATCH --time=24:00:00                 # Temps maximum

# Chargement des modules nécessaires
module load r/4.5.0                     # Version de R à charger

# Réglage des locales (optionnel)
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

# Exécution du script R avec l'identifiant de réplication
Rscript ProcessMod.R