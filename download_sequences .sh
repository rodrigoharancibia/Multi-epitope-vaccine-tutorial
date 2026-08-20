#!/usr/bin/env bash
set -euo pipefail

# ✅ FIX: Resolve 'datasets' binary with 3-level fallback
if [[ -x "/usr/local/bin/datasets" ]]; then
    DATASETS="/usr/local/bin/datasets"
elif [[ -x "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/datasets" ]]; then
    DATASETS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/datasets"
elif command -v datasets &>/dev/null; then
    DATASETS="$(command -v datasets)"
else
    echo "[ERROR] Instale datasets primeiro!"
    exit 1
fi

echo "[INFO] Using datasets binary: $DATASETS"

pathogen="${1:?ERROR: Forneça o patógeno. Ex: $0 'Klebsiella pneumoniae'}"
date_analysis=$(date +"%d-%m-%Y")
export BASE_DIR="$(pwd)"
export workdir="Sequences_${date_analysis}_${pathogen// /_}"
zipfile="${BASE_DIR}/${workdir}/${pathogen// /_}_${date_analysis}_data.zip"

echo "[INFO] Patógeno: ${pathogen}"
echo "[INFO] Diretório: ${workdir}"
mkdir -p "${BASE_DIR}/${workdir}"

echo "[INFO] Baixando proteomas RefSeq..."
for attempt in {1..5}; do
    echo "[INFO] Tentativa ${attempt}/5..."
    rm -f "${zipfile}"

    if "$DATASETS" download genome taxon "${pathogen}" \
        --assembly-level complete \
        --assembly-source RefSeq \
        --include protein \
        --filename "${zipfile}" 2>/dev/null && \
        unzip -t "${zipfile}" &>/dev/null; then

        echo "[INFO] Download e ZIP válidos!"
        break
    else
        rm -f "${zipfile}"
        [[ $attempt -lt 5 ]] && echo "Retry em 5s..." && sleep 5
    fi
done

[[ ! -f "${zipfile}" ]] && { echo "[ERROR] Download falhou"; exit 1; }

echo "[INFO] Extraindo ZIP..."
unzip -o "${zipfile}" -d "${BASE_DIR}/${workdir}"

echo "[INFO] Organizando arquivos .faa..."
DATA_DIR="${BASE_DIR}/${workdir}/ncbi_dataset/data"
faa_count=0
for f in "${DATA_DIR}"/GCF_*/protein.faa; do
    gcf=$(basename "$(dirname "$f")")
    dest="${BASE_DIR}/${pathogen// /_}_${gcf}_protein.faa"
    cp "$f" "${dest}"
    echo "[INFO] $(basename ${dest})"
    faa_count=$((faa_count + 1))
done

echo ""
echo "[SUCESSO] RESUMO:"
echo "[INFO] ${faa_count} arquivos .faa baixados"
echo "[INFO] Salvo em: ${BASE_DIR}"
echo "[INFO] Primeiros 5 arquivos:"
ls -lh "${BASE_DIR}"/*_protein.faa 2>/dev/null | head -5 || echo "[INFO] Nenhum arquivo encontrado"
echo ""
echo "[INFO] Para listar: ls *_protein.faa"

# ── Relatorio de metadados de assembly (assembly_metadata.csv) ──────────────
# Equivalente ao que download_sequences.py ja fazia: le o
# assembly_data_report.jsonl que o proprio 'datasets' baixa junto, e
# extrai accession, nivel de assembly, tamanho total, N50, pais, hospedeiro
# e data de coleta para cada genoma baixado neste patogeno.
#
# Requer 'jq' (parser de JSON em linha de comando).
echo ""
echo "[INFO] Gerando relatorio de metadados de assembly (assembly_metadata.csv)..."

if ! command -v jq &>/dev/null; then
    echo "[WARN] 'jq' nao encontrado no PATH -- pulando geracao de assembly_metadata.csv"
    echo "       (instale com: apt install jq / conda install -c conda-forge jq)"
else
    REPORT_FILE="${BASE_DIR}/assembly_metadata.csv"
    JSONL_FILE="${DATA_DIR}/assembly_data_report.jsonl"
    HEADER="pathogen,accession,assembly_level,total_bases,n50,country,host,collection_date"

    if [[ -f "${JSONL_FILE}" ]]; then
        NEW_ROWS="$(mktemp)"

        jq -r --arg pathogen "${pathogen}" '
            [
                $pathogen,
                (.accession // "N/A"),
                (.assemblyInfo.assemblyLevel // "N/A"),
                (.assemblyStats.totalSequenceLength // "N/A"),
                (.assemblyStats.contigN50 // "N/A"),
                (.assemblyInfo.biosample.geoLocName // "N/A"),
                (.assemblyInfo.biosample.host // "N/A"),
                (.assemblyInfo.biosample.collectionDate // "N/A")
            ] | @csv
        ' "${JSONL_FILE}" > "${NEW_ROWS}"

        if [[ ! -s "${NEW_ROWS}" ]]; then
            echo "[WARN] Nenhum registro valido em ${JSONL_FILE} -- relatorio nao atualizado"
            rm -f "${NEW_ROWS}"
        else
            # NOTA: dedup por accession (coluna 2) via awk, assumindo que o
            # campo accession nunca contem virgula (o formato do NCBI garante
            # isso). Se ja existir um registro com o mesmo accession, o novo
            # substitui o antigo -- mesmo comportamento do save_assembly_report()
            # em download_sequences.py (keep="last").
            if [[ -f "${REPORT_FILE}" ]]; then
                OLD_KEPT="$(mktemp)"
                awk -F',' 'NR==FNR { new_acc[$2]=1; next }
                           FNR==1 { next }
                           !($2 in new_acc) { print }' \
                    "${NEW_ROWS}" "${REPORT_FILE}" > "${OLD_KEPT}"

                { echo "${HEADER}"; cat "${OLD_KEPT}"; cat "${NEW_ROWS}"; } > "${REPORT_FILE}.tmp"
                mv "${REPORT_FILE}.tmp" "${REPORT_FILE}"
                rm -f "${OLD_KEPT}"
            else
                { echo "${HEADER}"; cat "${NEW_ROWS}"; } > "${REPORT_FILE}"
            fi

            rm -f "${NEW_ROWS}"
            echo "[INFO] Relatorio atualizado: ${REPORT_FILE}"
        fi
    else
        echo "[WARN] ${JSONL_FILE} nao encontrado -- relatorio nao gerado para este patogeno"
    fi
fi
