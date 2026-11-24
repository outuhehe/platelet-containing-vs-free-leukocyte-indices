# Optional exact-style renderer: Python + matplotlib 3.11.2 + Pillow.
# All plotting values are read from public aggregates, never hard-coded.
"""Figure 5: robustness forest plot of paired DeltaAUC across all analyses.

All point estimates and 95% CIs are taken verbatim from the reported aggregate
files (available_case_primary, common_sample_sensitivity, cluster_bootstrap_ci,
first_stay, emergency_stratified, reviewer_all_absolute_subgroup_ci, IPW
weighted/bootstrap outputs, and eICU paired delong). Direction is harmonised to
platelet-free minus platelet-containing, so positive values favour the
platelet-free version.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["axes.unicode_minus"] = False

# (label, NLR-vs-SII delta, lo, hi, SIRI-vs-PIV delta, lo, hi)
import csv, os
from pathlib import Path
RA = Path(os.getenv("RA_ROOT", Path(__file__).resolve().parent.parent))
TAB = RA / "outputs" / "tables"
OUTDIR = Path(os.getenv("FRONTMED_FIGURE_DIR", RA / "outputs" / "figures"))
OUTDIR.mkdir(parents=True, exist_ok=True)
def table(name):
    with (TAB/name).open(encoding="utf-8-sig",newline="") as f: return list(csv.DictReader(f))
def inverse(row): return (-float(row['delta_AUC']),-float(row['CI95_upper']),-float(row['CI95_lower']))
a=table('available_case_primary.csv'); b=table('common_sample_sensitivity.csv')
cl=table('cluster_bootstrap_ci.csv'); fs=table('first_stay_sensitivity.csv')
em=table('emergency_stratified.csv'); aa=table('reviewer_all_absolute_subgroup_ci.csv')
ip=table('ipw_weighted_primary.csv'); ib=table('ipw_bootstrap_ci.csv')
iy=table('ipw_outcome_conditioned_primary.csv'); iyb=table('ipw_outcome_conditioned_bootstrap_ci.csv')
ei=table('eicu_sicu_paired_delong_bootstrap.csv')
labels=['Primary available-case','Common five-index sample','Subject-clustered bootstrap',
        'First stay per subject','Emergency/urgent stratum','All-absolute subgroup',
        'IPW, baseline covariates','IPW, outcome-conditioned','eICU directional consistency']
columns=[]
for j in range(2):
    columns.append([inverse(a[j]),inverse(b[j]),
        (-float(b[j]['delta_AUC']),-float(cl[j]['cluster_CI95_upper']),-float(cl[j]['cluster_CI95_lower'])),
        inverse(fs[j]),inverse(em[j]),inverse(aa[j]),
        (float(ip[j]['weighted_delta_AUC']),float(ib[j]['CI95_lower']),float(ib[j]['CI95_upper'])),
        (float(iy[j]['weighted_delta_AUC']),float(iyb[j]['CI95_lower']),float(iyb[j]['CI95_upper'])),inverse(ei[j])])
DATA=[(label,*columns[0][i],*columns[1][i]) for i,label in enumerate(labels)]


labels = [d[0] for d in DATA]
y = list(range(len(DATA)))[::-1]  # first entry on top

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 6.4), sharey=True)
fig.subplots_adjust(left=0.20, right=0.985, top=0.88, bottom=0.13, wspace=0.12)

for ax, idx, title in (
    (ax1, (1, 2, 3), "NLR vs SII"),
    (ax2, (4, 5, 6), "SIRI vs PIV"),
):
    for yi, d in zip(y, DATA):
        dlt, lo, hi = d[idx[0]], d[idx[1]], d[idx[2]]
        ax.plot([lo, hi], [yi, yi], color="#17365D", linewidth=1.4, zorder=2)
        ax.plot([lo, lo], [yi - 0.12, yi + 0.12], color="#17365D", linewidth=1.2)
        ax.plot([hi, hi], [yi - 0.12, yi + 0.12], color="#17365D", linewidth=1.2)
        ax.scatter(dlt, yi, s=52, color="#17365D", zorder=3)
        ax.text(hi + 0.006, yi, f"{dlt:.3f} ({lo:.3f}, {hi:.3f})",
                va="center", ha="left", fontsize=8.3)
    ax.axvline(0, color="#C83E4D", linewidth=1.1, linestyle="--", zorder=1)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("\u0394AUC (95% CI)", fontsize=11)
    ax.set_xlim(-0.02, 0.245)
    ax.grid(axis="x", color="#E2E5E9", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

ax1.set_yticks(y)
ax1.set_yticklabels(labels, fontsize=10)
fig.suptitle("Robustness of the paired \u0394AUC across analyses", fontsize=13, y=0.96)
fig.text(0.20, 0.03,
         "Positive values favour the platelet-free index. Every interval excludes zero; the eICU estimates use a patient-level general surgical ICU cohort.",
         ha="left", va="bottom", fontsize=8.8)

OUT = str(OUTDIR / "Figure3_robustness_forest.png")
fig.savefig(OUT, dpi=300, facecolor="white")
plt.close(fig)
print("Wrote", OUT)
from PIL import Image
with Image.open(OUT) as im:
    im.convert("RGB").save(OUTDIR/"Figure3_robustness_forest.tif",compression="tiff_lzw",dpi=(300,300))
with (OUTDIR/"Figure3_plotting_values.csv").open("w",newline="",encoding="utf-8") as f:
    w=csv.writer(f);w.writerow(["analysis","NLR_delta","NLR_lower","NLR_upper","SIRI_delta","SIRI_lower","SIRI_upper"]);w.writerows(DATA)
