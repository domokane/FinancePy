# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 18:42:20 2026

@author: Dominic
"""

import matplotlib.pyplot as plt
from cycler import cycler

cc = [
    "#4477AA",
    "#EE6677",
    "#228833",
    "#CCBB44",
    "#66CCEE",
    "#AA3377",
    "#BBBBBB",
]

plt.rcParams.update(
    {
        "font.size": 16,
        "figure.figsize": (12, 6),
        "figure.dpi": 150,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.7,
        "axes.labelsize": 18,
        "axes.titlesize": 16,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "grid.linewidth": 0.6,
        "grid.alpha": 0.3,
        "legend.frameon": False,
        "legend.fontsize": 16,
        "lines.linewidth": 3,
        "lines.markersize": 10,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "figure.constrained_layout.use": True,
        "axes.prop_cycle": cycler(color=cc),
    }
)
