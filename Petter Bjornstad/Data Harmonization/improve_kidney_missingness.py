"""
IMPROVE kidney variable missingness analysis.
Checks missingness by visit × procedure for the kidney-focused variable
list used in IMPROVE analyses.

Accounts for column renaming and computed variables in improve.py.

Output: HTML report saved to ~/lab-qc-reports/reports/
"""

import os
import pandas as pd
import numpy as np
from datetime import date

BASE = os.path.expanduser(
    "~/Library/CloudStorage/OneDrive-UW/"
    "Laura Pyle's files - Biostatistics Core Shared Drive/"
    "Data Harmonization/"
)
REPORT_DIR = os.path.expanduser("~/lab-qc-reports/reports/")
os.makedirs(REPORT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
df = pd.read_csv(
    BASE + "Data Clean/harmonized_dataset.csv",
    dtype={"record_id": str},
    low_memory=False,
)
improve = df[df["study"] == "IMPROVE"].copy()
print(f"IMPROVE rows: {len(improve)}, participants: {improve['record_id'].nunique()}")
print(f"Visits: {sorted(improve['visit'].dropna().unique())}")
print(f"Procedures: {sorted(improve['procedure'].dropna().unique())}")

# ---------------------------------------------------------------------------
# Variable mapping: R field name → harmonized column name, procedure, notes
# ---------------------------------------------------------------------------
# Format: (r_field_name, harmonized_col, procedure, notes)
# notes = "computed" if derived in improve.py; "dropped" if removed (was a checkbox/flag)
FIELD_MAP = [
    # --- Screening ---
    ("screen_serum_creatinine",         "creatinine_s",           "screening",      "renamed"),
    ("screen_urine_mab",                "microalbumin_u",         "screening",      "renamed"),
    ("screen_urine_cre",                "creatinine_u",           "screening",      "renamed"),
    ("screen_urine_acr",                "acr_u",                  "screening",      "renamed"),

    # --- Clamp (raw collected) ---
    ("cystatin_c",                      "cystatin_c_s",           "clamp",          "renamed"),
    ("serum_creatinine",                "creatinine_s",           "clamp",          "renamed"),
    ("clamp_urine_mab_baseline",        "microalbumin_u",         "clamp",          "renamed"),
    ("clamp_urine_mab_250",             "urine_mab_250",          "clamp",          "renamed"),
    ("clamp_urine_cre_baseline",        "creatinine_u",           "clamp",          "renamed"),
    ("clamp_urine_cre_250",             "urine_cre_250",          "clamp",          "renamed"),
    ("clamp_acr_baseline",              "acr_u",                  "clamp",          "renamed"),
    ("clamp_acr_250",                   "acr_u_pm",               "clamp",          "renamed"),
    ("clamp_urine_sodium",              "sodium_u",               "clamp",          "renamed"),
    ("clamp_glucose_bl",                "urine_glucose_bl",       "clamp",          "renamed"),
    ("urine_glucose",                   "urine_glucose",          "clamp",          "as-is"),
    ("urine_osmol",                     "urine_osmol",            "clamp",          "as-is"),
    ("clamp_urine_vol",                 "urine_vol",              "clamp",          "renamed"),
    ("iohexol_bolus",                   "(dropped — was checkbox)","clamp",         "dropped"),
    ("iohexol_time",                    "iohexol_time",           "clamp",          "as-is"),
    ("iohexol_vol",                     "iohexol_vol",            "clamp",          "as-is"),
    ("iohexol_minus_10",                "iohexol_minus_10",       "clamp",          "as-is"),
    ("iohexol_120",                     "iohexol_120",            "clamp",          "as-is"),
    ("iohexol_150",                     "iohexol_150",            "clamp",          "as-is"),
    ("iohexol_180",                     "iohexol_180",            "clamp",          "as-is"),
    ("iohexol_210",                     "iohexol_210",            "clamp",          "as-is"),
    ("iohexol_240",                     "iohexol_240",            "clamp",          "as-is"),
    ("pah_bolus",                       "(dropped — was checkbox)","clamp",         "dropped"),
    ("pah_time",                        "pah_time",               "clamp",          "as-is"),
    ("pah_vol",                         "pah_vol",                "clamp",          "as-is"),
    ("pah_minus_10",                    "pah_minus_10",           "clamp",          "as-is"),
    ("pah_90",                          "pah_90",                 "clamp",          "as-is"),
    ("pah_120",                         "pah_120",                "clamp",          "as-is"),

    # --- Outcomes (computed & renamed; merged into clamp procedure) ---
    ("egfr",                            "(dropped — section flag)","clamp",         "dropped"),
    ("gfr",                             "gfr_raw_plasma",         "clamp",          "renamed"),
    ("gfr_bsa",                         "gfr_bsa_plasma",         "clamp",          "renamed"),
    ("abs_pah",                         "pah_clear_abs",          "clamp",          "renamed"),
    ("pah_bsa",                         "pah_clear_bsa",          "clamp",          "renamed"),
    ("rpf",                             "erpf_raw_plasma",        "clamp",          "renamed"),
    ("erpf_bsa",                        "erpf_bsa_plasma",        "clamp",          "renamed"),
    ("ecv",                             "ecv",                    "clamp",          "as-is"),
    ("gfr_ecv_percent",                 "gfr_ecv_percent",        "clamp",          "as-is"),
    ("gfr_ecv_std",                     "gfr_ecv_std",            "clamp",          "as-is"),
    ("ff",                              "ff",                     "clamp",          "computed: gfr_raw_plasma / erpf_raw_plasma"),
    ("rbf",                             "rbf",                    "clamp",          "computed: erpf_raw_plasma / (1 - hct/100)"),
    ("filtration_pressure",             "deltapf",                "clamp",          "computed: (gfr_raw_plasma/60) / kfg"),
    ("glomerular_pressure",             "glomerular_pressure",    "clamp",          "computed: pg + deltapf + 10"),
    ("glomerular_onco_pressure",        "pg",                     "clamp",          "computed: 5*(cm-2)"),
    ("mean_plasma_protein",             "cm",                     "clamp",          "computed: (total_protein/ff)*log(1/(1-ff))"),
    ("aff_arteriolar_resistance",       "ra",                     "clamp",          "computed: ((map - glom_pressure) / rbf_s) * 1328"),
    ("eff_arteriolar_resistance",       "re",                     "clamp",          "computed: (gfr_s / (kfg * (rbf_s - gfr_s))) * 1328"),

    # --- MRI (procedure = bold_mri) ---
    ("asl_right",                       "pcasl3d_right",          "bold_mri",       "renamed"),
    ("asl_left",                        "pcasl3d_left",           "bold_mri",       "renamed"),
    ("bold_r_bl_cortex",                "bold_r_bl_cortex",       "bold_mri",       "as-is"),
    ("bold_r_bl_medulla",               "bold_r_bl_medulla",      "bold_mri",       "as-is"),
    ("bold_r_bl_kidney",                "bold_r_bl_kidney",       "bold_mri",       "as-is"),
    ("bold_r_pf_cortex",                "bold_r_pf_cortex",       "bold_mri",       "as-is"),
    ("bold_r_pf_medulla",               "bold_r_pf_medulla",      "bold_mri",       "as-is"),
    ("bold_r_pf_kidney",                "bold_r_pf_kidney",       "bold_mri",       "as-is"),
    ("bold_l_bl_cortex",                "bold_l_bl_cortex",       "bold_mri",       "as-is"),
    ("bold_l_bl_medulla",               "bold_l_bl_medulla",      "bold_mri",       "as-is"),
    ("bold_l_bl_kidney",                "bold_l_bl_kidney",       "bold_mri",       "as-is"),
    ("bold_l_pf_cortex",                "bold_l_pf_cortex",       "bold_mri",       "as-is"),
    ("bold_l_pf_medulla",               "bold_l_pf_medulla",      "bold_mri",       "as-is"),
    ("bold_l_pf_kidney",                "bold_l_pf_kidney",       "bold_mri",       "as-is"),
    ("adc_right",                       "adc_right",              "bold_mri",       "as-is"),
    ("adc_left",                        "adc_left",               "bold_mri",       "as-is"),
    ("volume_right",                    "volume_right",           "bold_mri",       "as-is"),
    ("volume_left",                     "volume_left",            "bold_mri",       "as-is"),
    ("volume_right_manual",             "volume_right_manual",    "bold_mri",       "as-is"),
    ("volume_left_manual",              "volume_left_manual",     "bold_mri",       "as-is"),

    # --- Biopsy ---
    ("gloms",                           "gloms",                  "kidney_biopsy",  "as-is"),
    ("gloms_gs",                        "gloms_gs",               "kidney_biopsy",  "as-is"),
    ("ifta",                            "ifta",                   "kidney_biopsy",  "as-is"),
    ("glom_enlarge",                    "glom_enlarge",           "kidney_biopsy",  "as-is"),
    ("fia",                             "fia",                    "kidney_biopsy",  "as-is"),
    ("glom_tuft_area",                  "glom_tuft_area",         "kidney_biopsy",  "as-is"),
    ("glom_volume_weibel",              "glom_volume_weibel",     "kidney_biopsy",  "as-is"),
    ("glom_volume_wiggins",             "glom_volume_wiggins",    "kidney_biopsy",  "as-is"),
    ("glom_volume_con",                 "glom_volume_con",        "kidney_biopsy",  "as-is"),
    ("mes_matrix_area",                 "mes_matrix_area",        "kidney_biopsy",  "as-is"),
    ("mes_index",                       "mes_index",              "kidney_biopsy",  "as-is"),
    ("mes_volume_weibel",               "mes_volume_weibel",      "kidney_biopsy",  "as-is"),
    ("mes_volume_wiggins",              "mes_volume_wiggins",     "kidney_biopsy",  "as-is"),
    ("mes_volume_con",                  "mes_volume_con",         "kidney_biopsy",  "as-is"),
    ("gbm_thick_artmean",               "gbm_thick_artmean",      "kidney_biopsy",  "as-is"),
    ("gbm_thick_harmmean",              "gbm_thick_harmmean",     "kidney_biopsy",  "as-is"),
    ("pod_nuc_density",                 "pod_nuc_density",        "kidney_biopsy",  "as-is"),
    ("pod_cell_volume",                 "pod_cell_volume",        "kidney_biopsy",  "as-is"),
    ("art_intima",                      "art_intima",             "kidney_biopsy",  "as-is"),
    ("art_media",                       "art_media",              "kidney_biopsy",  "as-is"),

    # --- MMTT ---
    ("mmtt_bun_base",                   "bun_base",               "mmtt",           "renamed"),
    ("mmtt_creat_base",                 "creat_base",             "mmtt",           "renamed"),

    # --- AZ urine metabolites (procedure = az_u_metab) ---
    ("az_creatine_p",                   "az_creatine_p",          "az_u_metab",     "as-is"),
    ("az_taurine_p",                    "az_taurine_p",           "az_u_metab",     "as-is"),
    ("az_kynurenic_acid_n",             "az_kynurenic_acid_n",    "az_u_metab",     "as-is"),
    ("az_kynurenine_p",                 "az_kynurenine_p",        "az_u_metab",     "as-is"),
    ("az_trimethylamine_n_oxide_tmao_p","az_trimethylamine_n_oxide_tmao_p","az_u_metab","as-is"),
    ("az_hippuric_acid_p",              "az_hippuric_acid_p",     "az_u_metab",     "as-is"),
    ("az_oxalic_acid_n",                "az_oxalic_acid_n",       "az_u_metab",     "as-is"),

    # --- Covariates (medications procedure) ---
    ("uric_acid_med",                   "uric_acid_med",          "medications",    "as-is"),
    ("mra_med",                         "mra_med",                "medications",    "as-is"),
    ("htn_med_type",                    "htn_med_type",           "medications",    "as-is"),
]

# ---------------------------------------------------------------------------
# Check which harmonized columns actually exist
# ---------------------------------------------------------------------------
all_cols = set(improve.columns)
not_found = []
for r_name, h_col, proc, note in FIELD_MAP:
    if note.startswith("dropped"):
        continue
    if h_col not in all_cols:
        not_found.append((r_name, h_col, proc, note))

if not_found:
    print("\n=== Columns NOT found in harmonized dataset ===")
    for r, h, p, n in not_found:
        print(f"  R field: {r:40s}  harmonized: {h:35s}  procedure: {p}  ({n})")

# ---------------------------------------------------------------------------
# Missingness by visit × procedure
# ---------------------------------------------------------------------------
visits = ["baseline", "3_months_post_surgery", "12_months_post_surgery"]

rows = []
for r_name, h_col, proc, note in FIELD_MAP:
    if note.startswith("dropped"):
        row = {
            "r_field": r_name,
            "harmonized_col": "(dropped — was checkbox/flag)",
            "procedure": proc,
            "note": note,
        }
        for v in visits:
            row[f"n_{v}"] = "—"
            row[f"pct_missing_{v}"] = "—"
        rows.append(row)
        continue

    row = {
        "r_field": r_name,
        "harmonized_col": h_col,
        "procedure": proc,
        "note": note,
    }
    if h_col not in all_cols:
        for v in visits:
            row[f"n_{v}"] = "col not found"
            row[f"pct_missing_{v}"] = "col not found"
        rows.append(row)
        continue

    for v in visits:
        sub = improve[(improve["visit"] == v) & (improve["procedure"] == proc)]
        n_total = len(sub)
        if n_total == 0:
            row[f"n_{v}"] = 0
            row[f"pct_missing_{v}"] = "no rows"
        else:
            n_missing = sub[h_col].isna().sum()
            pct = n_missing / n_total * 100
            row[f"n_{v}"] = n_total
            row[f"pct_missing_{v}"] = round(pct, 1)
    rows.append(row)

results = pd.DataFrame(rows)

# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------
print("\n" + "=" * 100)
print("IMPROVE Kidney Variable Missingness (% missing by visit)")
print("=" * 100)

grouped = results.groupby("procedure")
for proc, grp in grouped:
    print(f"\n{'─'*100}")
    print(f"  Procedure: {proc}")
    print(f"{'─'*100}")
    header = f"  {'R field':<42} {'Harmonized col':<35} {'Baseline':>12} {'3mo post':>12} {'12mo post':>12}  Note"
    print(header)
    print(f"  {'-'*40} {'-'*33} {'-'*12} {'-'*12} {'-'*12}")
    for _, r in grp.iterrows():
        b  = str(r.get("pct_missing_baseline", "—"))
        m3 = str(r.get("pct_missing_3_months_post_surgery", "—"))
        m12= str(r.get("pct_missing_12_months_post_surgery", "—"))
        print(f"  {r['r_field']:<42} {r['harmonized_col']:<35} {b:>12} {m3:>12} {m12:>12}  {r['note']}")

# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------
today = date.today().strftime("%Y%m%d")
out_html = os.path.join(REPORT_DIR, f"IMPROVE_kidney_missingness_{today}.html")

def color_cell(val):
    if not isinstance(val, (int, float)):
        return f'<td style="background:#f0f0f0;text-align:center">{val}</td>'
    if val == 0:
        return f'<td style="background:#d4edda;text-align:center">0%</td>'
    elif val <= 10:
        return f'<td style="background:#d4edda;text-align:center">{val}%</td>'
    elif val <= 30:
        return f'<td style="background:#fff3cd;text-align:center">{val}%</td>'
    else:
        return f'<td style="background:#f8d7da;text-align:center">{val}%</td>'

html_parts = ["""
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>IMPROVE Kidney Variable Missingness</title>
<style>
  body { font-family: Arial, sans-serif; font-size: 13px; margin: 24px; }
  h1 { font-size: 18px; }
  h2 { font-size: 15px; margin-top: 28px; border-bottom: 2px solid #444; padding-bottom: 4px; }
  table { border-collapse: collapse; width: 100%; margin-bottom: 18px; }
  th { background: #343a40; color: #fff; padding: 6px 10px; text-align: left; }
  td { border: 1px solid #dee2e6; padding: 5px 10px; }
  tr:nth-child(even) td { background: #f8f9fa; }
  .note { font-size: 11px; color: #666; }
  .legend { display: flex; gap: 18px; margin-bottom: 14px; font-size: 12px; }
  .leg { padding: 3px 10px; border-radius: 3px; }
</style>
</head>
<body>
""",
f"<h1>IMPROVE Kidney Variable Missingness Report — {date.today().strftime('%Y-%m-%d')}</h1>",
f"<p>Participants: <b>{improve['record_id'].nunique()}</b> &nbsp;|&nbsp; Total rows: <b>{len(improve)}</b></p>",
"""<div class="legend">
  <span class="leg" style="background:#d4edda">0–10% missing</span>
  <span class="leg" style="background:#fff3cd">11–30% missing</span>
  <span class="leg" style="background:#f8d7da">&gt;30% missing</span>
  <span class="leg" style="background:#f0f0f0">N/A or dropped</span>
</div>
"""]

proc_order = ["screening", "clamp", "bold_mri", "kidney_biopsy", "mmtt", "az_u_metab", "medications"]
proc_labels = {
    "screening": "Screening Labs",
    "clamp": "Clamp (raw + outcomes + computed)",
    "bold_mri": "Bold/ASL/ADC/Volume MRI",
    "kidney_biopsy": "Kidney Biopsy",
    "mmtt": "MMTT",
    "az_u_metab": "AZ Urine Metabolites",
    "medications": "Covariates (Medications)",
}

for proc in proc_order:
    grp = results[results["procedure"] == proc]
    if grp.empty:
        continue
    label = proc_labels.get(proc, proc)
    html_parts.append(f"<h2>{label}</h2>")
    html_parts.append("""<table>
<tr>
  <th>R field name</th>
  <th>Harmonized column</th>
  <th>Baseline<br><small>(n rows)</small></th>
  <th>3mo post-surgery<br><small>(n rows)</small></th>
  <th>12mo post-surgery<br><small>(n rows)</small></th>
  <th>Note</th>
</tr>""")
    for _, r in grp.iterrows():
        b_n  = r.get("n_baseline", "—")
        m3_n = r.get("n_3_months_post_surgery", "—")
        m12_n= r.get("n_12_months_post_surgery", "—")
        b    = r.get("pct_missing_baseline", "—")
        m3   = r.get("pct_missing_3_months_post_surgery", "—")
        m12  = r.get("pct_missing_12_months_post_surgery", "—")
        n_label_b   = f" (n={b_n})"  if isinstance(b_n, int)  else ""
        n_label_m3  = f" (n={m3_n})" if isinstance(m3_n, int) else ""
        n_label_m12 = f" (n={m12_n})"if isinstance(m12_n, int)else ""
        html_parts.append(
            f"<tr>"
            f"<td>{r['r_field']}</td>"
            f"<td><b>{r['harmonized_col']}</b></td>"
            f"{color_cell(b)}{'' if not isinstance(b_n, int) else ''}"
        )
        # rebuild row properly
        html_parts.pop()  # remove partial row
        b_cell   = color_cell(b).replace("</td>", f"<br><small>{n_label_b}</small></td>") if isinstance(b_n, int) else color_cell(b)
        m3_cell  = color_cell(m3).replace("</td>", f"<br><small>{n_label_m3}</small></td>") if isinstance(m3_n, int) else color_cell(m3)
        m12_cell = color_cell(m12).replace("</td>", f"<br><small>{n_label_m12}</small></td>") if isinstance(m12_n, int) else color_cell(m12)
        html_parts.append(
            f"<tr>"
            f"<td>{r['r_field']}</td>"
            f"<td><b>{r['harmonized_col']}</b></td>"
            f"{b_cell}{m3_cell}{m12_cell}"
            f'<td class="note">{r["note"]}</td>'
            f"</tr>"
        )
    html_parts.append("</table>")

html_parts.append("</body></html>")

with open(out_html, "w") as f:
    f.write("\n".join(html_parts))

print(f"\nHTML report written: {out_html}")

# Also save CSV
out_csv = os.path.join(REPORT_DIR, f"IMPROVE_kidney_missingness_{today}.csv")
results.to_csv(out_csv, index=False)
print(f"CSV written:         {out_csv}")
