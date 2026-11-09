"""Plot helpers shared by the hub-exclusion figure scripts."""

import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D

import hxcommon  # noqa: F401  (sets sys.path, CH2_VARIANT)
from common import (COND_COLORS, COND_SHORT, CONDS, GROUP_COLORS, GROUPS, PT,
                    axes_in, color_ticklabels)
from hxcommon import PCT, VCOLOR, VKEYS, VLABEL, studies

_p = Path(__file__).resolve().parents[1] / "04_jaccard_heatmaps.py"
_s = importlib.util.spec_from_file_location("jacc", _p)
jacc = importlib.util.module_from_spec(_s)
_s.loader.exec_module(jacc)

XPOS = {v: i for i, v in enumerate(VKEYS)}      # thresholds at even spacing
GREY = "#dcdbd6"


def threshold_axis(ax, label=True):
    ax.set_xticks(range(len(VKEYS)))
    ax.set_xticklabels([str(PCT[v]) for v in VKEYS])
    ax.set_xlim(-0.35, len(VKEYS) - 0.65)
    if label:
        ax.set_xlabel("% of regulators excluded")


def cond_handles():
    return [Line2D([], [], color=COND_COLORS[c], marker="o", ms=3, label=COND_SHORT[c])
            for c in CONDS]


def variant_handles(marker="o"):
    return [Line2D([], [], ls="", marker=marker, color=VCOLOR[v], ms=4,
                   label=VLABEL[v]) for v in VKEYS]


# ---- overlap heatmaps (fig_h1 A and its supplement) ----------------------

def deg0_order(pairs, annot):
    """Trait order of 04_jaccard_heatmaps.py (--order shared) on deg0:
    average-linkage clustering of the mean observed Jaccard over conditions."""
    st = studies()
    ms = [matrix(pairs, "deg0", c, annot, "observed", st) for c in CONDS]
    return [st[i] for i in jacc.cluster_order(np.nan_to_num(np.mean(ms, axis=0)))]


def matrix(pairs, variant, cond, annot, col, order):
    idx = {s: i for i, s in enumerate(order)}
    m = np.full((len(order), len(order)), np.nan)
    sub = pairs[(pairs.variant == variant) & (pairs.condition == cond) & (pairs.annot == annot)]
    for a, b, v in zip(sub.study_a, sub.study_b, sub[col]):
        m[idx[a], idx[b]] = m[idx[b], idx[a]] = v
    np.fill_diagonal(m, 1.0)
    return m


def _heatmap(ax, m, order, diseases, vmin, vmax, cmap, ylabels, outline):
    n = len(order)
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad(GREY)
    shown = m.copy()
    np.fill_diagonal(shown, np.nan)
    im = ax.imshow(np.ma.masked_invalid(shown), cmap=cm, vmin=vmin, vmax=vmax,
                   interpolation="nearest", aspect="auto")
    # diagonal: neutral, not "missing"
    ax.imshow(np.where(np.eye(n, dtype=bool), 1.0, np.nan),
              cmap=ListedColormap(["#f1f1ef"]), vmin=0, vmax=1, aspect="auto")
    ax.set_xticks([])
    ax.set_yticks(range(n))
    if ylabels:
        ax.set_yticklabels([diseases.loc[s, "abbrev"] for s in order], fontsize=5.5)
        color_ticklabels(ax.get_yticklabels(), order, diseases)
    else:
        ax.set_yticklabels([])
    ax.tick_params(length=0, pad=1)
    for sp in ax.spines.values():
        sp.set_visible(False)
    codes = np.array([[GROUPS.index(diseases.loc[s, "group"]) for s in order]])
    strip = ListedColormap([GROUP_COLORS[g] for g in GROUPS])
    ax.imshow(codes, cmap=strip, vmin=0, vmax=1, aspect="auto",
              extent=(-0.5, n - 0.5, -0.85, -0.6))
    ax.imshow(codes.T, cmap=strip, vmin=0, vmax=1, aspect="auto",
              extent=(-0.85, -0.6, n - 0.5, -0.5))
    ax.set_xlim(-0.9, n - 0.5)
    ax.set_ylim(n - 0.5, -0.9)
    from matplotlib.patches import Rectangle
    for i, j in outline:
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="#0b0b0b",
                               lw=0.5, zorder=5))
    return im


def heatmap_block(fig, x0, y0, side, gap, pairs, sizes, cond, annot, diseases):
    """Two rows x four variants: observed (top) and excess Jaccard (bottom),
    deg0 trait order, one shared colour scale per row. Empty sets are grey;
    target sets that lost >= 1 seeding anchor get an open circle on the
    diagonal. Returns the bottom of the block (inches)."""
    order = deg0_order(pairs, annot)
    pos = {s: i for i, s in enumerate(order)}
    sub = pairs[(pairs.condition == cond) & (pairs.annot == annot)]
    sz = sizes[(sizes.condition == cond) & (sizes.annot == annot)]
    obs_max = sub.observed.max()
    span = np.nanmax(np.abs(sub.excess))
    outline = jacc.same_gwas_pairs(order, diseases)
    rows = [("observed", 0, obs_max, "Blues", "Observed Jaccard"),
            ("excess", -span, span, "RdBu_r", "Observed − expected")]
    for r, (col, vmin, vmax, cmap, label) in enumerate(rows):
        y = y0 - r * (side + 0.42)
        for k, v in enumerate(VKEYS):
            ax = axes_in(fig, x0 + k * (side + gap), y, side, side)
            m = matrix(pairs, v, cond, annot, col, order)
            empty = sz[(sz.variant == v) & (sz.n_genes == 0)].study
            for s in empty:                     # empty / lost set -> grey
                m[pos[s], :] = m[:, pos[s]] = np.nan
                m[pos[s], pos[s]] = 1.0
            im = _heatmap(ax, m, order, diseases, vmin, vmax, cmap, k == 0, outline)
            lost = sz[(sz.variant == v) & sz.lost_anchor].study
            ax.scatter([pos[s] for s in lost], [pos[s] for s in lost], s=7,
                       facecolor="none", edgecolor="#0b0b0b", lw=0.5, zorder=6)
            if r == 0:
                ax.set_title(f"{VLABEL[v]} excluded", fontsize=PT["tick"], pad=5)
        cax = axes_in(fig, x0 + 4 * (side + gap) + 0.0, y + 0.15 * side, 0.06, 0.7 * side)
        cb = fig.colorbar(im, cax=cax)
        cb.outline.set_linewidth(0.5)
        cb.ax.tick_params(labelsize=PT["tick"] - 1, length=2, width=0.5, pad=1)
        cb.set_label(label, fontsize=PT["tick"] - 1, labelpad=1)
        cb.ax.yaxis.set_major_locator(plt.MaxNLocator(4))
    return y0 - side - 0.42


# ---- forest columns (fig_h2 and its supplement) -------------------------

DODGE = {v: (i - 1.5) * 0.18 for i, v in enumerate(VKEYS)}


def forest(ax, h, order, value, se, sig, xlabel, diseases, ylabels=True,
           sig_label=True, zero=0.0):
    """Trait rows (top to bottom = order), one dodged point per variant with a
    95 % CI. Filled = ``sig``; diamond = target set that lost >= 1 seeding
    anchor. Non-fitted cells: grey square at the left edge of the slot."""
    from common import shade_groups
    from hcommon import flagged_name  # noqa: F401
    n = len(order)
    ypos = {s: n - 1 - i for i, s in enumerate(order)}
    shade_groups(ax, order[::-1], diseases, label=False)
    ax.axvline(zero, color="#9a9994", lw=0.5, zorder=1)
    for v in VKEYS:
        d = h[h.variant == v]
        for _, r in d.iterrows():
            if r.study not in ypos:
                continue
            y = ypos[r.study] - DODGE[v]
            if not np.isfinite(r[value]):
                ax.scatter(0.015, y, s=6, marker="s", color=GREY, lw=0, zorder=2,
                           transform=ax.get_yaxis_transform(), clip_on=False)
                continue
            c = VCOLOR[v]
            if np.isfinite(r[se]):
                ax.plot([r[value] - 1.96 * r[se], r[value] + 1.96 * r[se]], [y, y],
                        color=c, lw=0.6, zorder=2, solid_capstyle="butt")
            filled = bool(r[sig]) if sig else False
            ax.scatter(r[value], y, s=9, marker="D" if r.get("lost_anchor", False) else "o",
                       facecolor=c if filled else "white", edgecolor=c, lw=0.6, zorder=3)
    ax.set_ylim(-0.5, n - 0.5)
    ax.set_yticks(range(n))
    if ylabels:
        ax.set_yticklabels([diseases.loc[s, "display_name"] for s in order[::-1]])
        color_ticklabels(ax.get_yticklabels(), order[::-1], diseases)
    else:
        ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel(xlabel)
