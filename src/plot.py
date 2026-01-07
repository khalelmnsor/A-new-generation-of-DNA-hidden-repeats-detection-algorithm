import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import MultipleLocator

# ============================================================
# Setup
# ============================================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(BASE_DIR)

# ============================================================
# Helper functions
# ============================================================
def cdf(values):
    x = np.sort(values)
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y


def beautify_axes(
    ax,
    x_major=None, x_minor=None,
    y_major=None, y_minor=None
):
    if x_major is not None:
        ax.xaxis.set_major_locator(MultipleLocator(x_major))
    if x_minor is not None:
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor))

    if y_major is not None:
        ax.yaxis.set_major_locator(MultipleLocator(y_major))
    if y_minor is not None:
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor))

    ax.grid(True, which="major", alpha=0.4)
    ax.grid(True, which="minor", alpha=0.15)


def logP(df):
    return -np.log10(df[df["pvalue"] > 0]["pvalue"])


# ============================================================
# Load data
# ============================================================
avg_real  = pd.read_csv("segments_real_avg.csv")
avg_null  = pd.read_csv("segments_null_avg.csv")
pv_real   = pd.read_csv("segments_real_pvalue.csv")
pv_null   = pd.read_csv("segments_null_pvalue.csv")

# keep valid p-values
avg_real = avg_real[avg_real["pvalue"] > 0]
avg_null = avg_null[avg_null["pvalue"] > 0]
pv_real  = pv_real[pv_real["pvalue"] > 0]
pv_null  = pv_null[pv_null["pvalue"] > 0]

# ============================================================
# 1. CDF — AVG pipeline (Real vs Null)
# ============================================================
x_r, y_r = cdf(avg_real["pvalue"])
x_n, y_n = cdf(avg_null["pvalue"])

plt.figure(figsize=(8,6))
plt.plot(x_r, y_r, linewidth=2, label="Real segments (AVG)")
plt.plot(x_n, y_n, "--", linewidth=2, label="Null segments (AVG)", color="black")

plt.xlabel("P-value")
plt.ylabel("CDF")
plt.title("CDF of P-values — AVG-based segmentation")

ax = plt.gca()
beautify_axes(ax, x_major=0.05, x_minor=0.01, y_major=0.1, y_minor=0.05)

plt.legend()
plt.tight_layout()
plt.show()

# ============================================================
# 2. CDF — P-value pipeline (Real vs Null)
# ============================================================
x_r, y_r = cdf(pv_real["pvalue"])
x_n, y_n = cdf(pv_null["pvalue"])

plt.figure(figsize=(8,6))
plt.plot(x_r, y_r, linewidth=2, label="Real segments (P-value)")
plt.plot(x_n, y_n, "--", linewidth=2, label="Null segments (P-value)", color="black")

plt.xlabel("P-value")
plt.ylabel("CDF")
plt.title("CDF of P-values — P-value–based segmentation")

ax = plt.gca()
beautify_axes(ax, x_major=0.05, x_minor=0.01, y_major=0.1, y_minor=0.05)

plt.legend()
plt.tight_layout()
plt.show()

# ============================================================
# 3. CDF — Comparison of segmentation strategies (Real only)
# ============================================================
x_avg, y_avg = cdf(avg_real["pvalue"])
x_pv,  y_pv  = cdf(pv_real["pvalue"])

plt.figure(figsize=(8,6))
plt.plot(x_avg, y_avg, linewidth=2, label="AVG-based segmentation")
plt.plot(x_pv,  y_pv,  linewidth=2, label="P-value–based segmentation")

plt.xlabel("P-value")
plt.ylabel("CDF")
plt.title("Comparison of segmentation strategies (Real segments)")

ax = plt.gca()
beautify_axes(ax, x_major=0.05, x_minor=0.01, y_major=0.1, y_minor=0.05)

plt.legend()
plt.tight_layout()
plt.show()

# ============================================================
# 4. Histogram of significance (-log10 P)
# ============================================================
avg_real_lp = logP(avg_real)
avg_null_lp = logP(avg_null)
pv_real_lp  = logP(pv_real)
pv_null_lp  = logP(pv_null)

bins = np.linspace(
    0,
    max(avg_real_lp.max(), pv_real_lp.max()),
    40
)

plt.figure(figsize=(12,5))

# AVG pipeline
plt.subplot(1,2,1)
plt.hist(avg_real_lp, bins=bins, alpha=0.8, label="Real", edgecolor="black")
plt.hist(avg_null_lp, bins=bins, alpha=0.6, label="Null", edgecolor="black")

plt.xlabel("-log10(P-value)")
plt.ylabel("Number of segments")
plt.title("AVG-based segmentation")

ax = plt.gca()
beautify_axes(ax, x_major=2, x_minor=1, y_major=5, y_minor=1)
plt.legend()

# P-value pipeline
plt.subplot(1,2,2)
plt.hist(pv_real_lp, bins=bins, alpha=0.8, label="Real", edgecolor="black")
plt.hist(pv_null_lp, bins=bins, alpha=0.6, label="Null", edgecolor="black")

plt.xlabel("-log10(P-value)")
plt.ylabel("Number of segments")
plt.title("P-value–based segmentation")

ax = plt.gca()
beautify_axes(ax, x_major=2, x_minor=1, y_major=5, y_minor=1)
plt.legend()

plt.tight_layout()
plt.show()

# ============================================================
# 5. Scatter — AVG vs Significance
# ============================================================
avg_real["logP"] = -np.log10(avg_real["pvalue"])
pv_real["logP"]  = -np.log10(pv_real["pvalue"])

plt.figure(figsize=(8,6))

plt.scatter(
    avg_real["AVG"], avg_real["logP"],
    s=35, alpha=0.7, label="AVG-based pipeline"
)

plt.scatter(
    pv_real["AVG"], pv_real["logP"],
    s=35, alpha=0.7, label="P-value–based pipeline"
)

plt.xlabel("AVG match (%)")
plt.ylabel("-log10(P-value)")
plt.title("Relationship between AVG match and statistical significance")

ax = plt.gca()
beautify_axes(ax, x_major=2, x_minor=1, y_major=2, y_minor=1)

plt.legend()
plt.tight_layout()
plt.show()

print("\nAll figures generated successfully.")
# ============================================================
# 6. Segment length vs statistical significance (STRONG FIGURE)
# ============================================================

plt.figure(figsize=(8,6))

# compute logP if not exists
avg_real["logP"] = -np.log10(avg_real["pvalue"])
pv_real["logP"]  = -np.log10(pv_real["pvalue"])

plt.scatter(
    avg_real["length"], avg_real["logP"],
    s=35, alpha=0.7,
    label="AVG-based pipeline"
)

plt.scatter(
    pv_real["length"], pv_real["logP"],
    s=35, alpha=0.7,
    label="P-value–based pipeline"
)

plt.xlabel("Segment length (bp)")
plt.ylabel("-log10(P-value)")
plt.title("Segment length vs statistical significance")

ax = plt.gca()

# Axes control (קטן, ברור, אקדמי)
beautify_axes(
    ax,
    x_major=500, x_minor=250,   # תוכל לשנות לפי L
    y_major=2,   y_minor=1
)

plt.legend()
plt.tight_layout()
plt.show()
# ============================================================
# FIXED CDF DIFFERENCE — meaningful range only
# ============================================================

p_grid = np.logspace(-6, -1, 400)  # 1e-6 → 0.1

def cdf_interp(values, grid):
    x, y = cdf(values)
    return np.interp(grid, x, y, left=0, right=1)

diff_avg = (
    cdf_interp(avg_real["pvalue"], p_grid)
    - cdf_interp(avg_null["pvalue"], p_grid)
)

diff_pv = (
    cdf_interp(pv_real["pvalue"], p_grid)
    - cdf_interp(pv_null["pvalue"], p_grid)
)

plt.figure(figsize=(8,6))

plt.plot(p_grid, diff_avg, label="AVG-based pipeline")
plt.plot(p_grid, diff_pv, label="P-value–based pipeline")

plt.axhline(0, color="black", linestyle="--", linewidth=1)

plt.xscale("log")
plt.xlabel("P-value (log scale)")
plt.ylabel("CDF(real) − CDF(null)")
plt.title("Separation between real and null segments (significant range)")

ax = plt.gca()
beautify_axes(ax, y_major=0.05, y_minor=0.01)

plt.legend()
plt.tight_layout()
plt.show()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# ===============================
# Setup
# ===============================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(BASE_DIR)

# ===============================
# Load data
# ===============================
avg_real  = pd.read_csv("segments_real_avg.csv")
avg_null  = pd.read_csv("segments_null_avg.csv")
pv_real   = pd.read_csv("segments_real_pvalue.csv")
pv_null   = pd.read_csv("segments_null_pvalue.csv")

# keep valid p-values
def logP(df):
    return -np.log10(df[df["pvalue"] > 0]["pvalue"].values)

avg_real_lp = logP(avg_real)
avg_null_lp = logP(avg_null)
pv_real_lp  = logP(pv_real)
pv_null_lp  = logP(pv_null)

# ===============================
# Prepare Q–Q quantiles
# ===============================
def qq(real, null):
    n = min(len(real), len(null))
    return np.sort(null)[:n], np.sort(real)[:n]

x_avg, y_avg = qq(avg_real_lp, avg_null_lp)
x_pv,  y_pv  = qq(pv_real_lp,  pv_null_lp)

# ===============================
# Plot
# ===============================
plt.figure(figsize=(9,7))

plt.plot(
    x_avg, y_avg,
    marker="o", markersize=4,
    linewidth=2,
    alpha=0.8,
    label="AVG-based pipeline"
)

plt.plot(
    x_pv, y_pv,
    marker="s", markersize=4,
    linewidth=2,
    alpha=0.8,
    label="P-value–based pipeline"
)

# reference line
m = max(y_avg.max(), y_pv.max())
plt.plot([0, m], [0, m], "k--", label="y = x (null expectation)")

plt.xlabel("-log10(Null P-values)")
plt.ylabel("-log10(Real P-values)")
plt.title("Q–Q comparison of segmentation strategies\n(all segments length = 840)")

# axes detail
from matplotlib.ticker import MultipleLocator
ax = plt.gca()
ax.xaxis.set_major_locator(MultipleLocator(2))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(2))
ax.yaxis.set_minor_locator(MultipleLocator(1))

ax.grid(True, which="major", alpha=0.4)
ax.grid(True, which="minor", alpha=0.15)

plt.legend()
plt.tight_layout()
plt.show()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# ===============================
# Setup
# ===============================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(BASE_DIR)

# ===============================
# Load data
# ===============================
avg_real  = pd.read_csv("segments_real_avg.csv")
avg_null  = pd.read_csv("segments_null_avg.csv")
pv_real   = pd.read_csv("segments_real_pvalue.csv")
pv_null   = pd.read_csv("segments_null_pvalue.csv")

def get_logp(df):
    return -np.log10(df[df["pvalue"] > 0]["pvalue"].values)

avg_real_lp = get_logp(avg_real)
avg_null_lp = get_logp(avg_null)
pv_real_lp  = get_logp(pv_real)
pv_null_lp  = get_logp(pv_null)

# ===============================
# Enrichment computation
# ===============================
def enrichment_curve(real, null, thresholds):
    y = []
    for t in thresholds:
        r = np.sum(real >= t)
        n = np.sum(null >= t)
        y.append(r / max(n, 1))
    return np.array(y)

thresholds = np.linspace(0, max(
    avg_real_lp.max(),
    pv_real_lp.max()
), 60)

E_avg = enrichment_curve(avg_real_lp, avg_null_lp, thresholds)
E_pv  = enrichment_curve(pv_real_lp,  pv_null_lp,  thresholds)

# ===============================
# Plot
# ===============================
plt.figure(figsize=(9,6))

plt.plot(
    thresholds, E_avg,
    linewidth=2,
    label="AVG-based pipeline"
)

plt.plot(
    thresholds, E_pv,
    linewidth=2,
    label="P-value–based pipeline"
)

plt.axhline(1, color="black", linestyle="--", alpha=0.6)

plt.xlabel("-log10(P-value) threshold")
plt.ylabel("Enrichment (Real / Null)")
plt.title("Enrichment of statistically significant segments\n(segment length = 840)")

from matplotlib.ticker import MultipleLocator
ax = plt.gca()
ax.xaxis.set_major_locator(MultipleLocator(2))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.5))

ax.grid(True, which="major", alpha=0.4)
ax.grid(True, which="minor", alpha=0.15)

plt.legend()
plt.tight_layout()
plt.show()
import pandas as pd
import numpy as np
from scipy.stats import ks_2samp

# ===============================
# Load data
# ===============================
avg_real  = pd.read_csv("segments_real_avg.csv")
avg_null  = pd.read_csv("segments_null_avg.csv")
pv_real   = pd.read_csv("segments_real_pvalue.csv")
pv_null   = pd.read_csv("segments_null_pvalue.csv")

# keep valid p-values
def get_p(df):
    return df[df["pvalue"] > 0]["pvalue"].values

avg_real_p = get_p(avg_real)
avg_null_p = get_p(avg_null)
pv_real_p  = get_p(pv_real)
pv_null_p  = get_p(pv_null)

# ===============================
# KS tests
# ===============================
ks_avg = ks_2samp(avg_real_p, avg_null_p)
ks_pv  = ks_2samp(pv_real_p,  pv_null_p)

print("\n===============================")
print("KOLMOGOROV–SMIRNOV RESULTS")
print("===============================\n")

print("AVG-based pipeline:")
print(f"KS statistic = {ks_avg.statistic:.4f}")
print(f"KS p-value   = {ks_avg.pvalue:.2e}\n")

print("P-value–based pipeline:")
print(f"KS statistic = {ks_pv.statistic:.4f}")
print(f"KS p-value   = {ks_pv.pvalue:.2e}\n")

# ===============================
# Simple bar plot (academic & clean)
# ===============================
import matplotlib.pyplot as plt

plt.figure(figsize=(6,4))
plt.bar(
    ["AVG-based", "P-value–based"],
    [ks_avg.statistic, ks_pv.statistic],
    edgecolor="black"
)

plt.ylabel("KS distance")
plt.title("Separation between Real and Null distributions")

plt.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.show()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import MultipleLocator

# ============================================================
# Setup
# ============================================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(BASE_DIR)

# ============================================================
# Load REAL segments (only REAL has class)
# ============================================================
avg = pd.read_csv("segments_real_avg.csv")
pv  = pd.read_csv("segments_real_pvalue.csv")

# keep valid
avg = avg[avg["pvalue"] > 0].copy()
pv  = pv[pv["pvalue"] > 0].copy()

avg["logP"] = -np.log10(avg["pvalue"])
pv["logP"]  = -np.log10(pv["pvalue"])

classes = ["strong", "weak", "noise"]
colors  = {"strong":"red", "weak":"orange", "noise":"blue"}

# ============================================================
# Helper for clean axes
# ============================================================
def beautify(ax, xM=None, xm=None, yM=None, ym=None):
    if xM: ax.xaxis.set_major_locator(MultipleLocator(xM))
    if xm: ax.xaxis.set_minor_locator(MultipleLocator(xm))
    if yM: ax.yaxis.set_major_locator(MultipleLocator(yM))
    if ym: ax.yaxis.set_minor_locator(MultipleLocator(ym))
    ax.grid(True, which="major", alpha=0.4)
    ax.grid(True, which="minor", alpha=0.15)

# ============================================================
# 1. NUMBER OF SEGMENTS PER CLASS (comparison)
# ============================================================
cnt_avg = avg["class"].value_counts().reindex(classes, fill_value=0)
cnt_pv  = pv["class"].value_counts().reindex(classes, fill_value=0)

x = np.arange(len(classes))
w = 0.35

plt.figure(figsize=(7,5))
plt.bar(x - w/2, cnt_avg.values, w, label="AVG-based", edgecolor="black")
plt.bar(x + w/2, cnt_pv.values,  w, label="P-value-based", edgecolor="black")

plt.xticks(x, [c.upper() for c in classes])
plt.xlabel("Class")
plt.ylabel("Number of segments")
plt.title("Class distribution (REAL segments)")
plt.legend()
plt.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.show()

# ============================================================
# 2. AVG DISTRIBUTION PER CLASS (comparison)
# ============================================================
plt.figure(figsize=(9,5))

data_avg = [avg[avg["class"] == c]["AVG"] for c in classes]
data_pv  = [pv[pv["class"]  == c]["AVG"] for c in classes]

plt.boxplot(
    data_avg,
    positions=np.arange(len(classes)) - 0.2,
    widths=0.3,
    showfliers=False,
    labels=[c.upper() for c in classes]
)

plt.boxplot(
    data_pv,
    positions=np.arange(len(classes)) + 0.2,
    widths=0.3,
    showfliers=False
)

plt.xlabel("Class")
plt.ylabel("AVG match (%)")
plt.title("AVG distribution per class (REAL segments)")
plt.legend(
    [plt.Line2D([0],[0],color="black")]*2,
    ["AVG-based", "P-value-based"]
)
plt.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.show()

# ============================================================
# 3. SIGNIFICANCE DISTRIBUTION PER CLASS (comparison)
# ============================================================
plt.figure(figsize=(9,5))

data_avg = [avg[avg["class"] == c]["logP"] for c in classes]
data_pv  = [pv[pv["class"]  == c]["logP"] for c in classes]

plt.boxplot(
    data_avg,
    positions=np.arange(len(classes)) - 0.2,
    widths=0.3,
    showfliers=False
)

plt.boxplot(
    data_pv,
    positions=np.arange(len(classes)) + 0.2,
    widths=0.3,
    showfliers=False
)

plt.xticks(np.arange(len(classes)), [c.upper() for c in classes])
plt.xlabel("Class")
plt.ylabel("-log10(P-value)")
plt.title("Statistical significance per class (REAL segments)")
plt.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.show()

# ============================================================
# 4. AVG vs SIGNIFICANCE (comparison, colored by class)
# ============================================================
plt.figure(figsize=(8,6))

for c in classes:
    plt.scatter(
        avg[avg["class"] == c]["AVG"],
        avg[avg["class"] == c]["logP"],
        s=40, alpha=0.6,
        color=colors[c],
        marker="o",
        label=f"{c.upper()} (AVG)"
    )
    plt.scatter(
        pv[pv["class"] == c]["AVG"],
        pv[pv["class"] == c]["logP"],
        s=40, alpha=0.6,
        color=colors[c],
        marker="x",
        label=f"{c.upper()} (P-value)"
    )

plt.xlabel("AVG match (%)")
plt.ylabel("-log10(P-value)")
plt.title("AVG vs statistical significance (REAL segments)")
ax = plt.gca()
beautify(ax, xM=2, xm=1, yM=2, ym=1)
plt.legend(ncol=2, fontsize=8)
plt.tight_layout()
plt.show()

# ============================================================
# 5. CLASS PERCENTAGES (simple & strong)
# ============================================================
pct_avg = cnt_avg / cnt_avg.sum() * 100
pct_pv  = cnt_pv  / cnt_pv.sum()  * 100

plt.figure(figsize=(7,5))
plt.bar(x - w/2, pct_avg.values, w, label="AVG-based", edgecolor="black")
plt.bar(x + w/2, pct_pv.values,  w, label="P-value-based", edgecolor="black")

plt.xticks(x, [c.upper() for c in classes])
plt.ylabel("Percentage of segments (%)")
plt.title("Class proportions (REAL segments)")
plt.legend()
plt.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.show()

print("All comparative class-based figures generated successfully.")
