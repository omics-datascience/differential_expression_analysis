#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Utilidades para consultar la API de cBioPortal, identificar estudios con datos de RNA-seq
(mRNA expression) procesables con Limma.
"""


# https://chatgpt.com/share/68829186-32b4-8011-b6fe-76e71578da4c

import argparse
import csv
import logging
import time
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd
import requests

CBIOPORTAL_API_BASE = "https://www.cbioportal.org/api"
HEADERS = {"Accept": "application/json"}
PAGE_SIZE = 10_000_000

# Tipo de datos que necesito: CONTINUOUS
DATATYPE = {
    "CONTINUOUS"
}
# Otros tipos de datos: DISCRETE y Z-SCORE

LOG = logging.getLogger(__name__)


def get_cbioportal_studies(session: Optional[requests.Session] = None) -> pd.DataFrame:
    if session is None:
        session = requests.Session()

    url = (
        f"{CBIOPORTAL_API_BASE}/studies"
        f"?projection=SUMMARY&pageSize={PAGE_SIZE}&pageNumber=0&direction=ASC"
    )
    LOG.debug("GET %s", url)
    resp = session.get(url, headers=HEADERS, timeout=120)
    resp.raise_for_status()
    studies = resp.json()

    return pd.DataFrame(
        {
            "studyId": [s.get("studyId") for s in studies],
            "name": [s.get("name") for s in studies],
            "description": [s.get("description") for s in studies],
        }
    )


def get_molecular_profiles_rna_seq(
    study_id: str, session: Optional[requests.Session] = None
) -> pd.DataFrame:
    if session is None:
        session = requests.Session()

    url = (
        f"{CBIOPORTAL_API_BASE}/studies/{study_id}/molecular-profiles"
        f"?projection=SUMMARY&pageSize={PAGE_SIZE}&pageNumber=0"
    )
    LOG.debug("GET %s", url)
    resp = session.get(url, headers=HEADERS, timeout=120)
    resp.raise_for_status()
    profiles = resp.json()

    df = pd.DataFrame(profiles)
    if df.empty:
        return df

    # Filtro por nombres de perfiles de expresión mRNA por RNA-seq
    mask = (
        (df["molecularAlterationType"] == "MRNA_EXPRESSION") &
        (df["datatype"].isin(DATATYPE)) &
        (df["name"].str.contains("RSEM|TPM|FPKM", case=False, na=False))
    )
    # # Filtro por tipo de dato CONTINUOUS 
    # mask = (df["molecularAlterationType"] == "MRNA_EXPRESSION") & (
    #     df["datatype"].isin(DATATYPE)
    # )
    # Mostrar todo MRNA_EXPRESSION, no solo los nombres específicos
    # mask = df["molecularAlterationType"] == "MRNA_EXPRESSION"
    return df.loc[mask, [
        "molecularProfileId",
        "molecularAlterationType",
        "datatype",
        "name",
        "description",
        "studyId",
    ]].reset_index(drop=True)


def get_studies_with_rna_seq(
    output_tsv: Path,
    sleep_seconds: float = 0.05,
    session: Optional[requests.Session] = None,
) -> None:
    if session is None:
        session = requests.Session()

    studies = get_cbioportal_studies(session=session)
    study_ids = studies["studyId"].tolist()

    # Prepara TSV vacío con cabecera
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    write_header(output_tsv, header=[
        "molecularProfileId",
        "molecularAlterationType",
        "datatype",
        "name",
        "description",
        "studyId",
    ])

    for sid in study_ids:
        LOG.info("Checking study: %s", sid)
        try:
            df_profiles = get_molecular_profiles_rna_seq(sid, session=session)
            if not df_profiles.empty:
                append_df_to_tsv(df_profiles, output_tsv)
        except requests.HTTPError as e:
            LOG.error("HTTP error for study %s: %s", sid, e)
        except Exception as e:
            LOG.error("Error for study %s: %s", sid, e)

        time.sleep(sleep_seconds)


def write_header(path: Path, header: Iterable[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t", quoting=csv.QUOTE_NONE, escapechar="\\")
        writer.writerow(header)


def append_df_to_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False, header=False, mode="a")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Listar estudios de cBioPortal que contienen perfiles de expresión mRNA por RNA-seq."
    )
    parser.add_argument(
        "--output-tsv",
        type=Path,
        default=Path("cbioportal_datasets_with_rna_seq.tsv"),
        help="Ruta de salida del TSV con los perfiles RNA-seq encontrados.",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.05,
        help="Espera (en segundos) entre llamadas a la API.",
    )
    parser.add_argument(
        "-v", "--verbose", action="count", default=0, help="Aumenta la verbosidad del log."
    )
    return parser.parse_args()


def setup_logging(verbosity: int) -> None:
    level = logging.WARNING
    if verbosity == 1:
        level = logging.INFO
    elif verbosity >= 2:
        level = logging.DEBUG

    logging.basicConfig(
        level=level,
        format="[%(asctime)s] %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def main() -> None:
    args = parse_args()
    setup_logging(args.verbose)

    session = requests.Session()
    get_studies_with_rna_seq(
        output_tsv=args.output_tsv,
        sleep_seconds=args.sleep,
        session=session,
    )
    LOG.info("Finalizado. Archivo generado en: %s", args.output_tsv)


if __name__ == "__main__":
    main()
