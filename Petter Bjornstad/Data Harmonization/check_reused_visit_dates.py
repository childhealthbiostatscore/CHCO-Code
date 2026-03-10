"""
For each reused-visit pair from reused_visits.docx, find common procedures
and check whether the 'date' values match between the two (record_id, visit)
rows. Writes unmatched cases to unmatched_dates.csv.
"""

import re
import pandas as pd
from docx import Document

# ── 1. Parse pairs from docx ────────────────────────────────────────────────

def parse_pairs(docx_path):
    """Return deduplicated list of ((rid1, visit1), (rid2, visit2)) tuples.

    The docx uses blank lines to separate each explicit pair (sub-group of
    2 lines) within a UUID block.  We parse each sub-group as one pair and
    deduplicate so multi-subject UUIDs don't produce duplicate rows.
    """
    doc = Document(docx_path)
    line_re = re.compile(r'record_id\s+(\S+)\s+│\s+(\S+)')

    pairs_set = set()
    current_subgroup = []

    def flush_subgroup(sg):
        for i in range(len(sg)):
            for j in range(i + 1, len(sg)):
                # canonical order so (A,B) and (B,A) are the same key
                key = tuple(sorted([sg[i], sg[j]]))
                pairs_set.add(key)

    for para in doc.paragraphs:
        text = para.text.strip()
        m = line_re.search(text)
        if m:
            current_subgroup.append((m.group(1), m.group(2)))
        else:
            # blank line or separator — flush the current sub-group
            if current_subgroup:
                flush_subgroup(current_subgroup)
                current_subgroup = []

    if current_subgroup:
        flush_subgroup(current_subgroup)

    # unpack the canonical tuple back into two (rid, visit) tuples
    return [((a, b), (c, d)) for (a, b), (c, d) in pairs_set]


docx_path = "/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/reused_visits.docx"
csv_path  = "/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv"
out_path  = "/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/unmatched_dates.csv"

pairs = parse_pairs(docx_path)
print(f"Parsed {len(pairs)} unique pairs from docx.")

# ── 2. Load only the columns we need ────────────────────────────────────────

df = pd.read_csv(
    csv_path,
    usecols=["record_id", "visit", "procedure", "date"],
    low_memory=False,
    dtype=str,
)
df["date"]      = df["date"].str.strip()
df["visit"]     = df["visit"].str.strip()
df["record_id"] = df["record_id"].str.strip()
df["procedure"] = df["procedure"].str.strip()

print(f"Loaded {len(df)} rows from CSV.")

# ── 3. For each pair, find common procedures and compare dates ──────────────

rows_out = []

for (rid1, visit1), (rid2, visit2) in pairs:
    sub1 = df[(df["record_id"] == rid1) & (df["visit"] == visit1) & df["procedure"].notna()]
    sub2 = df[(df["record_id"] == rid2) & (df["visit"] == visit2) & df["procedure"].notna()]

    procs1 = set(sub1["procedure"].unique())
    procs2 = set(sub2["procedure"].unique())
    common = procs1 & procs2

    for proc in sorted(common):
        dates1 = sub1.loc[sub1["procedure"] == proc, "date"].dropna().unique()
        dates2 = sub2.loc[sub2["procedure"] == proc, "date"].dropna().unique()

        date1 = dates1[0] if len(dates1) == 1 else ("; ".join(sorted(dates1)) if len(dates1) > 1 else None)
        date2 = dates2[0] if len(dates2) == 1 else ("; ".join(sorted(dates2)) if len(dates2) > 1 else None)

        # Both None → neither has a date, not a mismatch
        if date1 is None and date2 is None:
            continue

        if date1 != date2:
            rows_out.append({
                "record_id_1": rid1,
                "visit_1":     visit1,
                "record_id_2": rid2,
                "visit_2":     visit2,
                "procedure":   proc,
                "date_1":      date1,
                "date_2":      date2,
            })

# ── 4. Save results ──────────────────────────────────────────────────────────

out = pd.DataFrame(rows_out, columns=[
    "record_id_1", "visit_1", "record_id_2", "visit_2",
    "procedure", "date_1", "date_2",
])
out = out.sort_values(["record_id_1", "visit_1", "record_id_2", "visit_2", "procedure"]).reset_index(drop=True)
out.to_csv(out_path, index=False)

print(f"\nFound {len(out)} unmatched date(s).")
print(out.to_string(index=False))
