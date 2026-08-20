#!/bin/bash

#######################################################################################
# Description : This script automates the process of running OrthoFinder and generates Core proteome
# Written by  : Rodrigo Alex Henríquez Arancibia
# Date Written: February, 2024
# Updated     : March, 2026
# Usage       : ./run_orthofinder.sh
#######################################################################################

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Real directory where script is located — works even when called via absolute path
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#######################################################################################
# Resolve OrthoFinder binary
#######################################################################################
ORTHOFINDER_BIN=""
ROOT_DIR=""

if [[ -x "$SCRIPT_DIR/OrthoFinder/orthofinder" ]]; then
    ORTHOFINDER_BIN="$SCRIPT_DIR/OrthoFinder/orthofinder"
    ROOT_DIR="$SCRIPT_DIR/OrthoFinder"
elif [[ -x "/opt/conda/envs/bio_env/bin/orthofinder" ]]; then
    ORTHOFINDER_BIN="/opt/conda/envs/bio_env/bin/orthofinder"
    ROOT_DIR="/opt/conda/envs/bio_env/bin"
elif command -v orthofinder &>/dev/null; then
    ORTHOFINDER_BIN="$(command -v orthofinder)"
    ROOT_DIR="$(dirname "$ORTHOFINDER_BIN")"
else
    echo "[ERROR] OrthoFinder not found! Install it via conda or place it in $SCRIPT_DIR/OrthoFinder/"
    exit 1
fi

echo "[INFO] Using OrthoFinder: $ORTHOFINDER_BIN"

#######################################################################################
# Resolve pangenome.pl and list.pl
#######################################################################################
PANGENOME_PL=""
LIST_PL=""

for _candidate in "$ROOT_DIR/pangenome.pl" "$SCRIPT_DIR/pangenome.pl" "/usr/local/bin/pangenome.pl"; do
    if [[ -f "$_candidate" ]]; then
        PANGENOME_PL="$_candidate"
        break
    fi
done

for _candidate in "$ROOT_DIR/list.pl" "$SCRIPT_DIR/list.pl" "/usr/local/bin/list.pl"; do
    if [[ -f "$_candidate" ]]; then
        LIST_PL="$_candidate"
        break
    fi
done

if [[ -z "$PANGENOME_PL" ]]; then
    echo "[ERROR] pangenome.pl not found! Place it alongside orthofinder or in /usr/local/bin/"
    exit 1
fi

if [[ -z "$LIST_PL" ]]; then
    echo "[ERROR] list.pl not found! Place it alongside orthofinder or in /usr/local/bin/"
    exit 1
fi

echo "[INFO] Using pangenome.pl: $PANGENOME_PL"
echo "[INFO] Using list.pl:      $LIST_PL"

#######################################################################################
# Configuration
#######################################################################################
TARGET_DIR="OrthoFinder_Results"
DATE=$(date +"%Y-%m-%d_%H-%M-%S")
seq_dir="$(pwd)/$TARGET_DIR/Sequences-$DATE"
mkdir -p "$TARGET_DIR"

# Capture the working directory (where .faa files are and where core.faa will be exported)
EXPORT_DIR="$(pwd)"

cd "$TARGET_DIR"

echo "##############################################"
echo " Pangenome Analysis with OrthoFinder"
echo "##############################################"

#######################################################################################
# Functions
#######################################################################################

setup_environment() {
    echo "############################"
    echo " Setting up environment ..."
    echo "############################"

    mkdir -p "$seq_dir"

    if mv ./*.faa "$seq_dir/" 2>/dev/null; then
        echo "✓ Files moved to $seq_dir"
    elif mv ../*.faa "$seq_dir/" 2>/dev/null; then
        echo "✓ Files moved to $seq_dir"
    else
        echo "[ERROR] No .faa files found!"
        exit 1
    fi

    echo "✓ Environment setup complete"
}

run_orthofinder() {
    echo "############################"
    echo " Running OrthoFinder ..."
    echo "############################"

    local count
    count=$(ls -1 "$seq_dir"/*.faa 2>/dev/null | wc -l)

    if [[ "$count" -lt 2 ]]; then
        echo "[ERROR] Found $count .faa file(s) in $seq_dir."
        echo "OrthoFinder requires at least 2 species to run."
        exit 1
    fi

    if ! "$ORTHOFINDER_BIN" -f "$seq_dir" -og; then
        echo "[ERROR] OrthoFinder execution failed!"
        exit 1
    fi
}

core_proteome() {
    echo "################################"
    echo " Generating Core Proteome ..."
    echo "################################"

    ORTHOGROUPS_DIR=$(find "$seq_dir" -type d -name "Orthogroups" -print -quit)

    if [[ -z "$ORTHOGROUPS_DIR" ]]; then
        echo "[ERROR] Orthogroups directory not found!"
        echo "OrthoFinder may have failed or changed its output structure"
        exit 1
    fi

    echo "[INFO] Found Orthogroups directory: $ORTHOGROUPS_DIR"

    echo "[INFO] Concatenating all protein sequences..."
    cat "$seq_dir"/*.faa > "$ORTHOGROUPS_DIR/all_proteins.faa"
    echo "✓ Created all_proteins.faa"

    cd "$ORTHOGROUPS_DIR" || exit 1

    echo "[INFO] Running pangenome.pl..."
    if ! perl "$PANGENOME_PL" Orthogroups.tsv Orthogroups_UnassignedGenes.tsv; then
        echo "[ERROR] pangenome.pl failed!"
        exit 1
    fi

    echo "[INFO] Processing core gene list..."
    cp core.txt core.list
    perl -pi -e "s/\,.+\n/\n/g" core.list
    perl -pi -e "s/\t.+\n/\n/g" core.list

    echo "[INFO] Extracting core proteome sequences..."
    if ! perl "$LIST_PL" core.list all_proteins.faa > core.faa; then
        echo "[ERROR] list.pl failed!"
        exit 1
    fi

    echo "✓ Core proteome generated successfully"
}

print_summary() {
    echo "##############################################"
    echo " Analysis Complete!"
    echo "##############################################"
    echo ""
    echo "Output files:"
    echo "  - Core proteome : $ORTHOGROUPS_DIR/core.faa"
    echo "  - Full results  : $ORTHOGROUPS_DIR"
    echo ""

    if [[ -f "$ORTHOGROUPS_DIR/core.faa" ]]; then
        local num_seqs
        num_seqs=$(grep -c "^>" "$ORTHOGROUPS_DIR/core.faa")
        echo "[INFO] Core proteome contains $num_seqs sequences"
    fi

    cp "$ORTHOGROUPS_DIR/core.faa" "$EXPORT_DIR/core.faa"
    echo "✓ core.faa exported to: $EXPORT_DIR/core.faa"

    echo ""
    echo "Finished successfully!"
}

#######################################################################################
# Main
#######################################################################################
main() {
    setup_environment
    run_orthofinder
    core_proteome
    print_summary
    echo "Analysis Finished Successfully!"
}

main "$@"