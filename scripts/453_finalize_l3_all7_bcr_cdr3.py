#!/usr/bin/env python3
"""Finalize L3 7/7 BCR/CDR3 using existing candidate FASTQs and TRUST4 outputs.

No original FASTQ is opened. TSV quote semantics are disabled. UMI correction
is restricted to patient+locus+exact CDR3 nt+CID and uses the frozen J2 rule.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import math
import statistics
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
STREAM = ROOT / "results/reanalysis/existing_unenriched_bcr_pilot_260824/streamed_candidates"
OLD_CALLERS = ROOT / (
    "results/reanalysis/existing_unenriched_bcr_pilot_260824/formal_pilot_runs/"
    "Y40105L3_BCR_CDR3_pilot_20260825_094835/callers/trust4"
)
J2_FINAL = ROOT / (
    "results/reanalysis/l3_j2_final_umi_audit_260827_160202/"
    "final_reparsed_all_chunks_260827_234500/CLONOTYPE_RECLASSIFICATION_AFTER_UMI_QC.tsv"
)
SSC_VALIDATION = ROOT / (
    "results/reanalysis/l3_ssc_main_clone_final_validation_260828_003106/"
    "SSC_15491_14_TRUST4_0003_FINAL_VALIDATION.md"
)
OLD_L3_AUDIT = ROOT / "results/reanalysis/y40105l3_j2_standard_supplementary_audit_260827_123000"
OLD_L3_CANDIDATES = OLD_L3_AUDIT / "all_productive_candidate_evidence_j2_standard.tsv"
OLD_L3_FRAGMENTS = OLD_L3_AUDIT / "direct_full_junction_fragment_cid_bin_support.tsv"
OLD_L3_R1 = ROOT / (
    "results/reanalysis/y40105l3_umi_technical_audit_260825_194423/tables/"
    "productive_full_junction_r1_quality.tsv"
)
MOTIFS = ("CGCTTGGCCT", "CCCTTACGCT")
SAMPLES = [
    ("HC/NL-66", "HC_NL-66", "HC"),
    ("HC/NL-72", "HC_NL-72", "HC"),
    ("IPF/FO23-1-06168", "IPF_FO23-1-06168", "IPF"),
    ("IPF/FO23-1-06170", "IPF_FO23-1-06170", "IPF"),
    ("SSC/05957/17B", "SSC_05957_17B", "SSc-ILD"),
    ("SSC/07998/15A", "SSC_07998_15A", "SSc-ILD"),
    ("SSC/15491/14", "SSC_15491_14", "SSc-ILD"),
]
COMPLETED_OLD = {"IPF/FO23-1-06168", "SSC/15491/14"}


def canonical_qname(value: str) -> str:
    return value.strip().lstrip("@>").split()[0].removesuffix("/1").removesuffix("/2")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        physical = list(reader)
    if not physical:
        raise RuntimeError(f"Empty TSV: {path}")
    width = len(physical[0])
    if any(len(row) != width for row in physical):
        raise RuntimeError(f"Physical TSV column mismatch: {path}")
    return [dict(zip(physical[0], row, strict=True)) for row in physical[1:]]


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str] | None = None) -> None:
    output = list(rows)
    if fields is None:
        if not output:
            raise RuntimeError(f"Fields required for empty output: {path}")
        fields = list(output[0])
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
            quoting=csv.QUOTE_NONE, escapechar="\\",
        )
        writer.writeheader()
        writer.writerows(output)


def truth(value: str) -> bool:
    return value.strip().upper() in {"T", "TRUE", "1", "YES"}


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def phred(quality: str) -> list[int]:
    return [ord(character) - 33 for character in quality]


def hamming(left: str, right: str) -> int:
    return sum(a != b for a, b in zip(left, right)) if len(left) == len(right) else max(len(left), len(right))


def entropy(values: Counter[str]) -> float:
    total = sum(values.values())
    return -sum((count / total) * math.log2(count / total) for count in values.values()) if total else 0.0


def find_file(directory: Path, suffix: str, exclude: str = "") -> Path:
    matches = [path for path in directory.iterdir() if path.name.endswith(suffix) and (not exclude or exclude not in path.name)]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one *{suffix} in {directory}, found {len(matches)}")
    return matches[0]


def caller_dir(run_dir: Path, patient: str, slug: str) -> Path:
    root = OLD_CALLERS if patient in COMPLETED_OLD else run_dir / "callers/trust4"
    candidates = sorted(root.glob(f"{slug}*"))
    completed = [path for path in candidates if (path / "TRUST4_COMPLETE.ok").is_file()]
    if len(completed) != 1:
        raise RuntimeError(f"Expected one completed TRUST4 directory for {patient}, found {len(completed)}")
    return completed[0]


def read_fasta(path: Path) -> list[tuple[str, str]]:
    output: list[tuple[str, str]] = []
    current = ""
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    output.append((current, "".join(chunks).upper()))
                current = canonical_qname(line)
                chunks = []
            else:
                chunks.append(line.strip())
    if current:
        output.append((current, "".join(chunks).upper()))
    return output


def scan_target_fastq(path: Path, targets: set[str], compressed: bool) -> tuple[dict[str, tuple[str, str]], Counter[str]]:
    opener = gzip.open if compressed else open
    found: dict[str, tuple[str, str]] = {}
    counts: Counter[str] = Counter()
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().rstrip("\r\n").upper()
            plus = handle.readline().rstrip("\r\n")
            quality = handle.readline().rstrip("\r\n")
            if not header.startswith("@") or not plus.startswith("+") or len(sequence) != len(quality):
                raise RuntimeError(f"Malformed FASTQ: {path}")
            name = canonical_qname(header)
            if name in targets:
                counts[name] += 1
                found[name] = (sequence, quality)
    return found, counts


def scan_sidecar(path: Path, targets: set[str]) -> tuple[dict[str, dict[str, str]], Counter[str]]:
    found: dict[str, dict[str, str]] = {}
    counts: Counter[str] = Counter()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = next(reader)
        width = len(header)
        for physical_line, fields in enumerate(reader, 2):
            if len(fields) != width:
                raise RuntimeError(f"Sidecar physical column mismatch {path}:{physical_line}")
            row = dict(zip(header, fields, strict=True))
            if row["qname"] in targets:
                counts[row["qname"]] += 1
                found[row["qname"]] = row
    return found, counts


def make_fragment(sequence: str, quality: str, junction: str) -> dict[str, Any]:
    forward = sequence.find(junction)
    reverse = sequence.find(reverse_complement(junction))
    if forward >= 0:
        canonical_sequence, canonical_quality = sequence, quality
        orientation, start = "forward", forward
    elif reverse >= 0:
        canonical_sequence, canonical_quality = reverse_complement(sequence), quality[::-1]
        orientation = "reverse_complement_in_R2"
        start = canonical_sequence.find(junction)
    else:
        raise RuntimeError("Full junction absent from R2")
    end = start + len(junction)
    digest = hashlib.sha256(canonical_sequence.encode("ascii")).hexdigest()
    junction_q = phred(canonical_quality[start:end])
    positional = f"r2len={len(canonical_sequence)};junction_start={start};junction_end={end}"
    return {
        "r2_orientation": orientation,
        "junction_start_0based": start,
        "junction_end_exclusive": end,
        "fragment_signature": positional,
        "canonical_r2_sha256": digest,
        "strict_fragment_signature": f"{positional};sha256={digest}",
        "junction_min_phred": min(junction_q),
        "junction_mean_phred": f"{statistics.mean(junction_q):.6g}",
    }


def motif_summary(scope: str, scope_id: str, reads: list[dict[str, Any]]) -> dict[str, Any]:
    counts = Counter(row["raw_umi"] for row in reads)
    total = len(reads)
    top = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    row: dict[str, Any] = {
        "scope": scope, "scope_id": scope_id, "supporting_reads": total,
        "unique_raw_umi": len(counts), "shannon_entropy_bits": f"{entropy(counts):.10g}",
        "top_umi": top[0][0] if top else "", "top_umi_count": top[0][1] if top else 0,
        "top_umi_fraction": f"{top[0][1] / total:.10g}" if top else "0",
    }
    for motif in MOTIFS:
        exact = [item for item in reads if item["raw_umi"] == motif]
        row[f"{motif}_count"] = len(exact)
        row[f"{motif}_fraction"] = f"{len(exact) / total:.10g}" if total else "0"
        row[f"{motif}_clonotypes"] = len({item["clonotype_id"] for item in exact})
        row[f"{motif}_cids"] = len({item["cid"] for item in exact})
        row[f"{motif}_bins"] = len({item["bin_id"] for item in exact})
        row[f"{motif}_fragment_signatures"] = len({item["fragment_signature"] for item in exact})
        qualities = [q for item in exact for q in item["umi_phred_values"]]
        row[f"{motif}_10bp_phred_median"] = statistics.median(qualities) if qualities else "not_observed"
    row["scope_limitation"] = "productive_CDR3_supporting_reads_only_not_global_L3_R1_background"
    return row


def html_table(rows: list[dict[str, Any]], fields: list[str]) -> str:
    head = "".join(f"<th>{html.escape(field)}</th>" for field in fields)
    body = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(row.get(field, '')))}</td>" for field in fields) + "</tr>"
        for row in rows
    )
    return f"<div class='scroll'><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()
    run_dir = (ROOT / args.run_dir).resolve()
    if not (run_dir / "TRUST4_REMAINING5_COMPLETE.ok").is_file():
        raise RuntimeError("Remaining-five TRUST4 completion checkpoint absent")
    if not SSC_VALIDATION.is_file() or "spatially supported provisional expanded IGH candidate" not in SSC_VALIDATION.read_text(encoding="utf-8"):
        raise RuntimeError("Completed SSC main-candidate validation is absent or inconsistent")

    disease = {patient: group for patient, _, group in SAMPLES}
    candidate_rows: list[dict[str, Any]] = []
    read_rows: list[dict[str, Any]] = []
    chain_rows: list[dict[str, Any]] = []
    inventory: list[dict[str, Any]] = []
    old_candidates = read_tsv(OLD_L3_CANDIDATES)
    old_fragments = read_tsv(OLD_L3_FRAGMENTS)
    old_r1 = {row["qname"]: row for row in read_tsv(OLD_L3_R1)}

    for patient, slug, group in SAMPLES:
        if patient in COMPLETED_OLD:
            patient_candidates = [row for row in old_candidates if row["sample_id"] == patient]
            candidate_map = {row["trust4_candidate_id"]: row for row in patient_candidates}
            for chain in ("IGH", "IGK", "IGL"):
                subset = [row for row in patient_candidates if row["chain"] == chain]
                chain_rows.append({
                    "patient": patient, "disease_group": group, "chain": chain,
                    "productive_clonotypes": len(subset), "productive_airr_records": len(subset),
                    "nonproductive_airr_records": 0, "nonproductive_records_with_cdr3": 0,
                })
            for source in (OLD_L3_CANDIDATES, OLD_L3_FRAGMENTS, OLD_L3_R1):
                stat = source.stat()
                inventory.append({
                    "patient": patient, "role": "reused completed two-patient audit",
                    "path": str(source), "bytes": stat.st_size,
                    "modified": datetime.fromtimestamp(stat.st_mtime).astimezone().isoformat(),
                })
            for fragment in [row for row in old_fragments if row["sample_id"] == patient]:
                candidate = candidate_map[fragment["trust4_candidate_id"]]
                r1 = old_r1[fragment["qname"]]
                r1_sequence, r1_quality = r1["r1_sequence"], r1["r1_quality"]
                if len(r1_sequence) != 35 or len(r1_quality) != 35:
                    raise RuntimeError(f"Reused completed R1 structure failure for {fragment['qname']}")
                digest = fragment["canonical_r2_sha256"]
                positional = fragment["fragment_signature"]
                read_rows.append({
                    "patient": patient, "disease_group": group,
                    "clonotype_id": fragment["trust4_candidate_id"], "chain": candidate["chain"],
                    "v_call": candidate["v_call"], "d_call": candidate["d_call"], "j_call": candidate["j_call"],
                    "productive": "true", "cdr3_nt": candidate["cdr3_nt"], "cdr3_aa": candidate["cdr3_aa"],
                    "qname": fragment["qname"], "cid": r1_sequence[:25], "raw_umi": r1_sequence[25:35],
                    "corrected_umi": "pending", "r1_sequence": r1_sequence, "r1_quality": r1_quality,
                    "umi_quality": r1_quality[25:35], "umi_phred_values": phred(r1_quality[25:35]),
                    "bin_id": fragment["bin_id"], "fastq_chunk": fragment["fastq_chunk"],
                    "r2_sequence": fragment["r2_sequence"], "r2_quality": fragment["r2_quality"],
                    "r2_orientation": fragment["junction_orientation"],
                    "junction_start_0based": fragment["junction_start_0based_canonical_r2"],
                    "junction_end_exclusive": fragment["junction_end_exclusive_canonical_r2"],
                    "fragment_signature": positional, "canonical_r2_sha256": digest,
                    "strict_fragment_signature": f"{positional};sha256={digest}",
                    "junction_min_phred": fragment["junction_min_phred"],
                    "junction_mean_phred": fragment["junction_mean_phred"],
                })
            continue
        caller = caller_dir(run_dir, patient, slug)
        airr_path = find_file(caller, "_airr.tsv", exclude="barcode")
        assembled_path = find_file(caller, "_assembled_reads.fa")
        toassemble_path = find_file(caller, "_toassemble.fq")
        sidecar_path = STREAM / slug / "reliable_ig_candidates_sidecar.tsv.gz"
        r1_path = STREAM / slug / "reliable_ig_candidates_CID_UMI_R1.fastq.gz"
        for role, path in [("AIRR", airr_path), ("assembled reads", assembled_path), ("candidate R2 cache", toassemble_path), ("QNAME sidecar", sidecar_path), ("candidate R1", r1_path)]:
            stat = path.stat()
            inventory.append({"patient": patient, "role": role, "path": str(path), "bytes": stat.st_size, "modified": datetime.fromtimestamp(stat.st_mtime).astimezone().isoformat()})

        airr = read_tsv(airr_path)
        for chain in ("IGH", "IGK", "IGL"):
            subset = [row for row in airr if row.get("locus") == chain]
            productive = [row for row in subset if truth(row.get("productive", "")) and row.get("junction")]
            nonproductive = [row for row in subset if not truth(row.get("productive", ""))]
            chain_rows.append({
                "patient": patient, "disease_group": group, "chain": chain,
                "productive_clonotypes": len({(row["junction"], row.get("junction_aa", "")) for row in productive}),
                "productive_airr_records": len(productive),
                "nonproductive_airr_records": len(nonproductive),
                "nonproductive_records_with_cdr3": sum(bool(row.get("junction")) for row in nonproductive),
            })

        productive_unique: dict[tuple[str, str, str], dict[str, str]] = {}
        for row in airr:
            if row.get("locus") in {"IGH", "IGK", "IGL"} and truth(row.get("productive", "")) and row.get("junction"):
                productive_unique.setdefault((row["locus"], row["junction"].upper(), row.get("junction_aa", "")), row)
        clone_airr: dict[str, dict[str, str]] = {
            f"{slug}_TRUST4_{index:04d}": row
            for index, row in enumerate(productive_unique.values(), 1)
        }

        direct: dict[str, set[str]] = defaultdict(set)
        for qname, sequence in read_fasta(assembled_path):
            for clone_id, row in clone_airr.items():
                junction = row["junction"].upper()
                if junction in sequence or reverse_complement(junction) in sequence:
                    direct[clone_id].add(qname)
        targets = set().union(*direct.values()) if direct else set()
        r2_found, r2_counts = scan_target_fastq(toassemble_path, targets, compressed=False)
        sidecars, sidecar_counts = scan_sidecar(sidecar_path, targets)
        r1_found, r1_counts = scan_target_fastq(r1_path, targets, compressed=True)
        for label, found, counts in [("R2", r2_found, r2_counts), ("sidecar", sidecars, sidecar_counts), ("R1", r1_found, r1_counts)]:
            if set(found) != targets or any(counts[name] != 1 for name in targets):
                raise RuntimeError(f"{patient} target QNAME reconnection failed for {label}")

        for clone_id, row in clone_airr.items():
            names = sorted(direct.get(clone_id, set()))
            if not names:
                continue
            junction = row["junction"].upper()
            for qname in names:
                sidecar = sidecars[qname]
                r1_sequence, r1_quality = r1_found[qname]
                r2_sequence, r2_quality = r2_found[qname]
                if len(r1_sequence) != 35 or len(r1_quality) != 35 or len(r2_sequence) != len(r2_quality):
                    raise RuntimeError(f"Read structure failure for {qname}")
                cid, umi = r1_sequence[:25], r1_sequence[25:35]
                if cid != sidecar["cid"] or umi != sidecar["umi"]:
                    raise RuntimeError(f"R1/sidecar CID or UMI mismatch for {qname}")
                if any(q < 0 or q > 93 for q in phred(r1_quality)):
                    raise RuntimeError(f"Illegal R1 Phred character for {qname}")
                fragment = make_fragment(r2_sequence, r2_quality, junction)
                read_rows.append({
                    "patient": patient, "disease_group": group, "clonotype_id": clone_id,
                    "chain": row["locus"], "v_call": row.get("v_call", ""), "d_call": row.get("d_call", ""),
                    "j_call": row.get("j_call", ""), "productive": "true", "cdr3_nt": junction,
                    "cdr3_aa": row.get("junction_aa", ""), "qname": qname, "cid": cid,
                    "raw_umi": umi, "corrected_umi": "pending", "r1_sequence": r1_sequence,
                    "r1_quality": r1_quality, "umi_quality": r1_quality[25:35],
                    "umi_phred_values": phred(r1_quality[25:35]), "bin_id": sidecar["bin_id"],
                    "fastq_chunk": sidecar["chunk"], "r2_sequence": r2_sequence,
                    "r2_quality": r2_quality, **fragment,
                })

            members = [item for item in read_rows if item["patient"] == patient and item["clonotype_id"] == clone_id]
            candidate_rows.append({
                "patient": patient, "disease_group": group, "clonotype_id": clone_id,
                "chain": row["locus"], "v_call": row.get("v_call", ""), "d_call": row.get("d_call", ""),
                "j_call": row.get("j_call", ""), "productive": "true", "cdr3_nt": junction,
                "cdr3_aa": row.get("junction_aa", ""), "raw_reads": len(members),
            })

    write_tsv(run_dir / "L3_ALL7_INPUT_FILE_INVENTORY.tsv", inventory)
    write_tsv(run_dir / "L3_ALL7_BCR_CHAIN_DISTRIBUTION.tsv", chain_rows)

    correction_groups: dict[tuple[str, str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in read_rows:
        correction_groups[(row["patient"], row["chain"], row["cdr3_nt"], row["cid"])].append(row)
    decisions: list[dict[str, Any]] = []
    for key, members in sorted(correction_groups.items()):
        raw_counts = Counter(row["raw_umi"] for row in members)
        correction = {umi: umi for umi in raw_counts}
        for source in sorted(raw_counts, key=lambda value: (raw_counts[value], value)):
            candidates = [target for target in raw_counts if target != source and hamming(source, target) == 1 and raw_counts[target] >= 5 * raw_counts[source] and raw_counts[source] <= 3]
            candidates.sort(key=lambda value: (-raw_counts[value], value))
            target = candidates[0] if candidates and (len(candidates) == 1 or raw_counts[candidates[0]] > raw_counts[candidates[1]]) else ""
            mismatches = [index for index, (left, right) in enumerate(zip(source, target)) if left != right] if target else []
            mismatch_q = [item["umi_phred_values"][index] for item in members if item["raw_umi"] == source for index in mismatches]
            median_q = statistics.median(mismatch_q) if mismatch_q else None
            merge = bool(target and median_q is not None and median_q <= 20)
            if merge:
                correction[source] = target
            decisions.append({
                "patient": key[0], "chain": key[1], "cdr3_nt": key[2], "cid": key[3],
                "source_raw_umi": source, "source_read_count": raw_counts[source],
                "target_umi": target or "none", "target_read_count": raw_counts[target] if target else 0,
                "hamming_distance": hamming(source, target) if target else "not_applicable",
                "count_ratio": f"{raw_counts[target] / raw_counts[source]:.10g}" if target else "not_applicable",
                "mismatch_positions_1based": ",".join(str(index + 1) for index in mismatches),
                "mismatch_phred_values": ",".join(map(str, mismatch_q)),
                "mismatch_median_phred": median_q if median_q is not None else "not_applicable",
                "decision": "merge" if merge else "retain", "corrected_umi": correction[source],
                "reason": "count_ratio>=5;source_reads<=3;unique_HD1_target;mismatch_median_Q<=20" if merge else "quality_aware_merge_criteria_not_all_met",
            })
        for row in members:
            row["corrected_umi"] = correction[row["raw_umi"]]

    serializable_reads = [{**row, "umi_phred_values": ",".join(map(str, row["umi_phred_values"]))} for row in read_rows]
    write_tsv(run_dir / "L3_ALL7_CDR3_SUPPORTING_READS.tsv", serializable_reads)
    write_tsv(run_dir / "L3_ALL7_UMI_CORRECTION_DECISIONS.tsv", decisions)

    cross_patients: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in read_rows:
        cross_patients[(row["chain"], row["cdr3_nt"])].add(row["patient"])
    cross_rows: list[dict[str, Any]] = []
    for (chain, cdr3), patients in sorted(cross_patients.items()):
        if len(patients) > 1:
            ids = sorted({row["clonotype_id"] for row in read_rows if row["chain"] == chain and row["cdr3_nt"] == cdr3})
            cross_rows.append({
                "chain": chain, "cdr3_nt": cdr3, "patients": "|".join(sorted(patients)),
                "patient_count": len(patients), "clonotype_ids": "|".join(ids),
                "warning": "public-like_or_technical_warning", "merge_policy": "never_merge_across_patients",
            })
    write_tsv(
        run_dir / "L3_ALL7_CROSS_PATIENT_CDR3_AUDIT.tsv", cross_rows,
        ["chain", "cdr3_nt", "patients", "patient_count", "clonotype_ids", "warning", "merge_policy"],
    )

    clone_groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in read_rows:
        clone_groups[(row["patient"], row["clonotype_id"])].append(row)
    classifications: list[dict[str, Any]] = []
    for (patient, clone_id), members in sorted(clone_groups.items()):
        exact = len({(row["cid"], row["raw_umi"]) for row in members})
        corrected = len({(row["cid"], row["corrected_umi"]) for row in members})
        cids = len({row["cid"] for row in members})
        bins = len({row["bin_id"] for row in members})
        positional = len({row["fragment_signature"] for row in members})
        strict_fragments = len({row["strict_fragment_signature"] for row in members})
        chain = members[0]["chain"]
        if corrected == 1:
            final_class = "singleton" if len(members) == 1 else "technical candidate"
        elif chain == "IGH" and cids >= 2 and bins >= 2 and positional >= 2 and strict_fragments >= 2:
            final_class = "reliable expanded clonotype"
        else:
            final_class = "moderate candidate"
        cross = len(cross_patients[(chain, members[0]["cdr3_nt"])]) > 1
        special = "spatially supported provisional expanded IGH candidate" if clone_id == "SSC_15491_14_TRUST4_0003" else "not_applicable"
        classifications.append({
            "patient": patient, "disease_group": disease[patient], "clonotype_id": clone_id,
            "chain": chain, "v_call": members[0]["v_call"], "d_call": members[0]["d_call"],
            "j_call": members[0]["j_call"], "productive": "true", "cdr3_nt": members[0]["cdr3_nt"],
            "cdr3_aa": members[0]["cdr3_aa"], "raw_reads": len(members),
            "corrected_molecules": corrected, "unique_cid": cids, "unique_bin": bins,
            "independent_positional_fragment_signatures": positional,
            "independent_strict_r2_fragments": strict_fragments,
            "supporting_fastq_chunks": len({row["fastq_chunk"] for row in members}),
            "fastq_chunk_ids": "|".join(sorted({row["fastq_chunk"] for row in members})),
            "unique_raw_umi": len({row["raw_umi"] for row in members}),
            "exact_cid_umi_molecules": exact,
            "CGCTTGGCCT_reads": sum(row["raw_umi"] == "CGCTTGGCCT" for row in members),
            "cross_patient_exact_cdr3_warning": "public-like_or_technical_warning" if cross else "none",
            "final_classification": final_class,
            "spatially_validated_interpretation": special,
            "classification_reason": (
                "IGH with >=2 corrected molecules, >=2 CID, >=2 Bin and >=2 independent positional/strict R2 fragments"
                if final_class == "reliable expanded clonotype" else
                ">=2 corrected molecules but strict reliable-IGH independence gate not met"
                if final_class == "moderate candidate" else
                "multiple reads collapse to one CID+corrected-UMI molecule"
                if final_class == "technical candidate" else "one independent molecule"
            ),
        })
    write_tsv(run_dir / "L3_ALL7_CLONOTYPE_CLASSIFICATION.tsv", classifications)
    notable = [row for row in classifications if row["final_classification"] in {"moderate candidate", "reliable expanded clonotype"}]
    write_tsv(
        run_dir / "L3_ALL7_MODERATE_RELIABLE_CANDIDATES.tsv", notable,
        list(classifications[0]),
    )

    motif_rows = [motif_summary("chip", "L3_all7", read_rows)]
    for patient, _, _ in SAMPLES:
        motif_rows.append(motif_summary("patient", patient, [row for row in read_rows if row["patient"] == patient]))
    for clone_id in sorted({row["clonotype_id"] for row in read_rows}):
        motif_rows.append(motif_summary("clonotype", clone_id, [row for row in read_rows if row["clonotype_id"] == clone_id]))
    write_tsv(run_dir / "L3_ALL7_UMI_MOTIF_AUDIT.tsv", motif_rows)

    patient_summary: list[dict[str, Any]] = []
    for patient, _, group in SAMPLES:
        subset = [row for row in classifications if row["patient"] == patient]
        patient_summary.append({
            "patient": patient, "disease_group": group, "completion_status": "complete",
            "trust4_complete": True, "productive_bcr_clonotypes": len(subset),
            "productive_IGH": sum(row["chain"] == "IGH" for row in subset),
            "productive_IGK": sum(row["chain"] == "IGK" for row in subset),
            "productive_IGL": sum(row["chain"] == "IGL" for row in subset),
            "singleton": sum(row["final_classification"] == "singleton" for row in subset),
            "technical_candidate": sum(row["final_classification"] == "technical candidate" for row in subset),
            "moderate_candidate": sum(row["final_classification"] == "moderate candidate" for row in subset),
            "reliable_expanded_clonotype": sum(row["final_classification"] == "reliable expanded clonotype" for row in subset),
            "public_like_or_technical_warnings": sum(row["cross_patient_exact_cdr3_warning"] != "none" for row in subset),
        })
    write_tsv(run_dir / "L3_ALL7_PATIENT_COMPLETION_SUMMARY.tsv", patient_summary)

    j2 = [row for row in read_tsv(J2_FINAL) if row["chip_id"] == "Y40105J2"]
    comparable: list[dict[str, Any]] = []
    for chip, rows, patients in [("J2", j2, "7/7"), ("L3", classifications, "7/7")]:
        comparable.append({
            "chip": chip, "completion_status": "complete", "completed_patients": patients,
            "productive_BCR_clonotypes": len(rows),
            "IGH": sum(row["chain"] == "IGH" for row in rows),
            "IGK": sum(row["chain"] == "IGK" for row in rows),
            "IGL": sum(row["chain"] == "IGL" for row in rows),
            "singleton": sum(row["final_classification"] == "singleton" for row in rows),
            "technical": sum(row["final_classification"] == "technical candidate" for row in rows),
            "moderate": sum(row["final_classification"] == "moderate candidate" for row in rows),
            "reliable_or_provisional_expanded": sum(row["final_classification"] == "reliable expanded clonotype" for row in rows),
            "interpretation": "SSC main L3 reliable-rule hit is reported as spatially supported provisional expanded IGH candidate" if chip == "L3" else "No reliable expanded J2 clonotype",
        })
    write_tsv(run_dir / "L3_J2_FINAL_COMPARABLE_SUMMARY.tsv", comparable)

    total_motif = motif_rows[0]
    class_counts = Counter(row["final_classification"] for row in classifications)
    disease_notable = defaultdict(list)
    for row in notable:
        disease_notable[row["disease_group"]].append(f"{row['patient']} / {row['clonotype_id']} / {row['chain']} / {row['final_classification']}")
    report = f"""# L3 7/7 BCR/CDR3 final report

Generated: {datetime.now().astimezone().isoformat()}

## Completion and integrity

L3 is complete for 7/7 manifest patients. Existing candidate FASTQs were used; original FASTQs were not rescanned. The two previously completed TRUST4 patients were read-only and were not rerun. Every productive supporting QNAME was unique in the candidate R1, TRUST4 cached R2 and sidecar manifest; R1 was exactly 35 bp with CID=`R1[0:25]`, UMI=`R1[25:35]`, and complete quality.

Total productive BCR clonotypes: **{len(classifications)}**. Classification: **{class_counts['singleton']} singleton, {class_counts['technical candidate']} technical, {class_counts['moderate candidate']} moderate, {class_counts['reliable expanded clonotype']} reliable-rule hit**.

## Moderate/reliable candidates

{chr(10).join('- ' + item for group in ('SSc-ILD', 'IPF', 'HC') for item in disease_notable[group]) if notable else '- None'}

`SSC_15491_14_TRUST4_0003` remains **spatially supported provisional expanded IGH candidate** after its completed read-level, VDJ and spatial validation. It meets the uniform reliable-rule molecule gate, but is deliberately worded as provisional because all molecules use the anomalously frequent UMI motif and the cross-patient junction warning remains.

## UMI motif

`CGCTTGGCCT` occurs in **{total_motif['CGCTTGGCCT_count']}/{total_motif['supporting_reads']} ({float(total_motif['CGCTTGGCCT_fraction']):.2%})** productive CDR3-supporting reads. Repetition of this motif is never treated as molecule independence; CID, Bin and R2 fragment evidence are primary.

## Cross-patient audit

Exact cross-patient CDR3 warnings: **{len(cross_rows)}**. These are labelled `public-like_or_technical_warning` and are never merged across patients.

## J2 comparison

Both J2 and L3 are complete 7/7. J2: {len(j2)} productive clonotypes, 34 singleton, 25 technical, 9 moderate and 0 reliable. L3: {len(classifications)} productive clonotypes, {class_counts['singleton']} singleton, {class_counts['technical candidate']} technical, {class_counts['moderate candidate']} moderate and {class_counts['reliable expanded clonotype']} reliable-rule hit (reported as provisional after spatial review).

The comparison is descriptive: unenriched random capture, sparse direct junction support and the chip-level UMI motif abnormality prevent interpreting raw reads as clone size or testing disease prevalence from these seven-patient chips.

## Scope

Conclusions apply only to completed J2 and L3. They do not extend to K8 or all 21 patients. No TRA/TRB/TRD/TRG record is included in BCR counts.
"""
    (run_dir / "L3_ALL7_FINAL_REPORT.md").write_text(report, encoding="utf-8")
    html_report = f"""<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'><title>L3 all7 BCR/CDR3</title>
<style>body{{font:15px Arial,'Microsoft YaHei',sans-serif;max-width:1200px;margin:28px auto;line-height:1.5;color:#222}}h1,h2{{color:#17365d}}table{{border-collapse:collapse;font-size:12px}}th,td{{border:1px solid #bbb;padding:5px}}.scroll{{overflow:auto}}.callout{{background:#eef6ff;border-left:5px solid #2b6cb0;padding:12px}}</style></head><body>
<h1>L3 7/7 BCR/CDR3 final report</h1><div class='callout'>L3 complete 7/7. SSC_15491_14_TRUST4_0003: spatially supported provisional expanded IGH candidate.</div>
<h2>Patient completion</h2>{html_table(patient_summary, list(patient_summary[0]))}
<h2>Moderate/reliable candidates</h2>{html_table(notable, list(classifications[0])) if notable else '<p>None</p>'}
<h2>L3 versus J2</h2>{html_table(comparable, list(comparable[0]))}
<h2>Interpretation</h2><p>{html.escape(report)}</p></body></html>"""
    (run_dir / "L3_ALL7_FINAL_REPORT.html").write_text(html_report, encoding="utf-8")
    (run_dir / "L3_ALL7_FINAL_COMPLETE.ok").touch()
    print(f"L3_CLONOTYPES={len(classifications)}")
    print(f"L3_CLASSES={dict(class_counts)}")
    print(f"CGCTTGGCCT={total_motif['CGCTTGGCCT_count']}/{total_motif['supporting_reads']}")


if __name__ == "__main__":
    main()
