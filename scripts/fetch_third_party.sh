#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TP_DIR="${ROOT_DIR}/third_party"
mkdir -p "${TP_DIR}"

download_file() {
    local url="$1"
    local output="$2"
    if [[ -f "${output}" ]]; then
        echo "[fetch] Using cached ${output##*/}"
    else
        echo "[fetch] Downloading ${output##*/} from ${url}"
        curl -L --fail --retry 3 "${url}" -o "${output}"
    fi
}

extract_dir() {
    local archive="$1"
    local target_dir="$2"
    local temp_dir
    temp_dir="$(mktemp -d)"
    tar -xzf "${archive}" -C "${temp_dir}"

    local candidate=""
    if [[ -d "${temp_dir}/MitoFates" ]]; then
        candidate="${temp_dir}/MitoFates"
    elif [[ -d "${temp_dir}/MitoFates_1.2/MitoFates" ]]; then
        candidate="${temp_dir}/MitoFates_1.2/MitoFates"
    elif [[ -d "${temp_dir}/mitoprotII" ]]; then
        candidate="${temp_dir}/mitoprotII"
    elif [[ -d "${temp_dir}/mitoprotII-v1.101" ]]; then
        candidate="${temp_dir}/mitoprotII-v1.101"
    fi

    if [[ -z "${candidate}" ]]; then
        echo "[fetch] Unable to detect extracted directory for ${archive}" >&2
        rm -rf "${temp_dir}"
        exit 1
    fi

    rm -rf "${target_dir}"
    mv "${candidate}" "${target_dir}"
    rm -rf "${temp_dir}"
}

# MitoFates
MITOFATES_TAR="${TP_DIR}/MitoFates_1.2.tar.gz"
download_file "https://mitf.cbrc.pj.aist.go.jp/MitoFates/program/MitoFates_1.2.tar.gz" "${MITOFATES_TAR}"
extract_dir "${MITOFATES_TAR}" "${TP_DIR}/MitoFates"

# Mitoprot II
MITOPROT_TAR="${TP_DIR}/mitoprotII-v1.101.tar.gz"
download_file "ftp://ftp.biologie.ens.fr/pub/molbio/mitoprotII-v1.101.tar.gz" "${MITOPROT_TAR}"
extract_dir "${MITOPROT_TAR}" "${TP_DIR}/mitoprotII"

echo "[fetch] Third-party assets ready under ${TP_DIR}"
