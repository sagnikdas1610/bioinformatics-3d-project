"""External structure lookup and AI-assisted prediction API clients (RCSB, AlphaFold DB, ESMFold)."""

from __future__ import annotations

import re
from typing import Literal

import requests

from bio3d.config import MAX_ESMFOLD_RESIDUES, MIN_ESMFOLD_RESIDUES
from bio3d.sequence import validate_sequence

RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.{fmt}"
ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
ESMFOLD_API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"


class ExternalServiceError(RuntimeError):
    """Raised when an external bioinformatics service returns an error or invalid payload."""


def fetch_rcsb_structure(pdb_id: str, fmt: Literal["pdb", "cif"] = "pdb") -> tuple[str, str]:
    """Fetch an experimental PDB or mmCIF structure by 4-character RCSB PDB ID."""
    clean_id = (pdb_id or "").strip().upper()
    if not re.fullmatch(r"[A-Z0-9]{4}", clean_id):
        raise ValueError("Enter a valid 4-character RCSB PDB ID (e.g., 1UBQ, 4HHB, 6LU7).")

    url = RCSB_DOWNLOAD_URL.format(pdb_id=clean_id, fmt=fmt)
    try:
        resp = requests.get(url, timeout=20)
    except requests.RequestException as exc:
        raise ExternalServiceError(f"Network error contacting RCSB PDB: {exc}") from exc

    if resp.status_code == 404:
        raise ExternalServiceError(f"PDB ID '{clean_id}' was not found on RCSB PDB (format: {fmt}).")
    if resp.status_code != 200:
        raise ExternalServiceError(f"RCSB PDB returned HTTP {resp.status_code} for {clean_id}.{fmt}.")

    text = resp.text
    marker = "_atom_site." if fmt == "cif" else "ATOM  "
    if marker not in text:
        raise ExternalServiceError(f"RCSB response for {clean_id}.{fmt} did not contain atomic coordinates.")
    return text, f"{clean_id}.{fmt}"


def fetch_alphafold_structure(uniprot_id: str) -> dict[str, str]:
    """Fetch a precomputed AlphaFold Protein Structure Database model by UniProt accession."""
    clean_acc = (uniprot_id or "").strip().upper()
    if not re.fullmatch(r"[A-Z0-9]{6,10}", clean_acc):
        raise ValueError("Enter a valid UniProt accession (e.g., P69905, P0CG48, P04637).")

    api_url = ALPHAFOLD_API_URL.format(uniprot_id=clean_acc)
    try:
        meta_resp = requests.get(api_url, timeout=20)
    except requests.RequestException as exc:
        raise ExternalServiceError(f"Network error contacting AlphaFold DB: {exc}") from exc

    if meta_resp.status_code == 404:
        raise ExternalServiceError(f"No AlphaFold DB entry found for UniProt accession '{clean_acc}'.")
    if meta_resp.status_code != 200:
        raise ExternalServiceError(f"AlphaFold DB API returned HTTP {meta_resp.status_code}.")

    try:
        entries = meta_resp.json()
    except ValueError as exc:
        raise ExternalServiceError("AlphaFold DB returned malformed JSON metadata.") from exc

    if not entries or not isinstance(entries, list):
        raise ExternalServiceError(f"AlphaFold DB returned no models for '{clean_acc}'.")

    entry = entries[0]
    pdb_url = entry.get("pdbUrl")
    if not pdb_url:
        raise ExternalServiceError(f"AlphaFold DB entry for '{clean_acc}' does not provide a PDB URL.")

    try:
        pdb_resp = requests.get(pdb_url, timeout=25)
    except requests.RequestException as exc:
        raise ExternalServiceError(f"Failed downloading AlphaFold PDB from {pdb_url}: {exc}") from exc

    if pdb_resp.status_code != 200 or "ATOM  " not in pdb_resp.text:
        raise ExternalServiceError(f"AlphaFold PDB download failed (HTTP {pdb_resp.status_code}).")

    entry_id = str(entry.get("entryId") or f"AF-{clean_acc}-F1")
    return {
        "uniprot_id": clean_acc,
        "entry_id": entry_id,
        "gene": str(entry.get("gene") or "Unknown"),
        "organism": str(entry.get("organismScientificName") or "Unknown"),
        "latest_version": str(entry.get("latestVersion") or ""),
        "pdb_url": str(pdb_url),
        "pdb_text": pdb_resp.text,
        "filename": f"{entry_id}.pdb",
    }


def predict_with_esmfold(sequence: str) -> dict[str, object]:
    """Run live single-chain de novo folding via the public ESMFold Atlas API."""
    validation = validate_sequence(sequence)
    if not validation.is_valid:
        raise ValueError(validation.message)

    seq = validation.cleaned_sequence
    if len(seq) < MIN_ESMFOLD_RESIDUES:
        raise ValueError(
            f"Sequence is too short ({len(seq)} aa). ESMFold requires at least {MIN_ESMFOLD_RESIDUES} residues."
        )
    if len(seq) > MAX_ESMFOLD_RESIDUES:
        raise ValueError(
            f"Sequence length ({len(seq)} aa) exceeds the public ESMFold API limit of "
            f"{MAX_ESMFOLD_RESIDUES} residues. Use ColabFold for longer sequences."
        )

    try:
        resp = requests.post(ESMFOLD_API_URL, data=seq, timeout=60)
    except requests.Timeout as exc:
        raise ExternalServiceError(
            "ESMFold API request timed out after 60s. Public servers may be busy; try a shorter sequence or ColabFold."
        ) from exc
    except requests.RequestException as exc:
        raise ExternalServiceError(f"Network error contacting ESMFold API: {exc}") from exc

    if resp.status_code in {403, 429, 502, 503, 504}:
        raise ExternalServiceError(
            f"Public ESMFold API is currently rate-limited or unavailable (HTTP {resp.status_code}). "
            "No fallback structure was substituted; use AlphaFold DB or the ColabFold notebook."
        )
    if resp.status_code != 200 or "ATOM  " not in resp.text:
        raise ExternalServiceError(
            f"ESMFold API returned HTTP {resp.status_code} without valid PDB ATOM coordinates."
        )

    return {
        "sequence": seq,
        "length": len(seq),
        "pdb_text": resp.text,
        "filename": f"esmfold_prediction_{len(seq)}aa.pdb",
        "service": "ESM Atlas ESMFold v1 API",
    }
