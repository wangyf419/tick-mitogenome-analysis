#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# 03_mcmc_diagnostic_figures.py
#
# Four-chain convergence diagnostics for the nine MCMCTREE datasets
# (Figure S1 and the local diagnostic series Figure_S1-S7). The script reads
# stored statistics and trace caches; no chain is sampled or summarised here.
#
# Run from the module directory:
#   python scripts/03_mcmc_diagnostic_figures.py
# ------------------------------------------------------------------------------

import json
import logging
import string
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, NullFormatter

import fitz

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/mcmc'
FIG = ROOT / 'outputs/mcmc'
MEDIA = FIG / 'previews'
FIG.mkdir(parents=True, exist_ok=True)
MEDIA.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 8,
    'font.weight': 'normal',
    'axes.labelsize': 8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'axes.linewidth': .75,
    'lines.linewidth': .75,
    'pdf.fonttype': 42,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# Subsetting messages from the font tools are irrelevant for these figures.
logging.getLogger('fontTools.subset').setLevel(logging.ERROR)

# Letter page and the drawing frame used by the original layout, in points.
W, H = 612., 792.
L, R, T, B = 51.0234, 561.259785821385, 51.0234, 738.510346558398
GAP = 18

COLORS = ['#222222', '#0072B2', '#E69F00', '#CC79A7']
STYLES = ['-', '--', '-.', ':']
MARKERS = ['o', 's', '^', 'D']
LABELS = ['Original', 'Chain 1', 'Chain 2', 'Chain 3']

# Stored statistics and trace caches.
TABLE = OUT / 'statistics_csv'
summary = pd.read_csv(TABLE / 'dataset_convergence_summary.csv')
diag = pd.read_csv(TABLE / 'all_parameter_diagnostics.csv')
chains = pd.read_csv(TABLE / 'per_chain_posterior_summary.csv')
comparison = pd.read_csv(TABLE / 'original_vs_added_chains.csv')
manifest = json.loads((OUT / 'manifest.json').read_text(encoding='utf8'))

# Figure registry written next to the PDFs.
registry = []

# Two columns inside the text frame.
COLS = [(L, (L + R - GAP) / 2), ((L + R + GAP) / 2, R)]


def newfig():
    """Return a new Letter-sized figure."""
    return plt.figure(figsize=(8.5, 11))


def legend(fig, markers=False):
    """Add the shared four-chain legend at the top of the page."""
    handles = [
        Line2D(
            [0], [0],
            color=COLORS[i], ls=STYLES[i], lw=.8,
            marker=MARKERS[i] if markers else None,
            markersize=3, label=LABELS[i],
        )
        for i in range(4)
    ]
    fig.legend(
        handles=handles, loc='upper center',
        bbox_to_anchor=((L + R) / 2 / W, 1 - (T + 1) / H),
        ncol=4, frameon=False, borderaxespad=0,
    )


def fit(fig, ax, cell, index, label=None, italic=False):
    """Place an axis inside its cell and shrink it until nothing overlaps."""
    x0, y0, x1, y1 = cell
    bounds = [x0 + 1, H - y1 + 1, x1 - 1, H - y0 - 18]
    ax.set_position([
        (bounds[0] + 40) / W,
        (bounds[1] + 25) / H,
        (bounds[2] - bounds[0] - 45) / W,
        (bounds[3] - bounds[1] - 29) / H,
    ])
    for _ in range(5):
        renderer = fig.canvas.get_renderer()
        scale = 72 / fig.dpi
        tight = ax.get_tightbbox(renderer)
        pos = ax.get_position()
        over = [
            max(0, pos.x0 * W - tight.x0 * scale),
            max(0, pos.y0 * H - tight.y0 * scale),
            max(0, tight.x1 * scale - pos.x1 * W),
            max(0, tight.y1 * scale - pos.y1 * H),
        ]
        nl = bounds[0] + over[0] + 1
        nb = bounds[1] + over[1] + 1
        nr = bounds[2] - over[2] - 1
        nt = bounds[3] - over[3] - 1
        assert nr - nl > 60 and nt - nb > 24, (index, label, nr - nl, nt - nb)
        ax.set_position([nl / W, nb / H, (nr - nl) / W, (nt - nb) / H])
    ax.tick_params(width=.75, length=3, pad=3)
    for spine in ax.spines.values():
        spine.set_linewidth(.75)
    fig.text(x0 / W, 1 - y0 / H, f'({string.ascii_uppercase[index]})',
             fontsize=10, va='top')
    if label:
        fig.text(
            (x0 + 26) / W, 1 - (y0 + 1) / H, label,
            fontsize=8, fontstyle='italic' if italic else 'normal', va='top',
        )


def save(fig, title, caption, kind, **extra):
    """Write one PDF, crop a PNG preview and append a registry record."""
    number = len(registry) + 1
    path = FIG / f'Figure_S{number}.pdf'
    fig.savefig(path, format='pdf')
    plt.close(fig)
    with fitz.open(path) as doc:
        page = doc[0]
        boxes = [
            fitz.Rect(box) for box_type, box in page.get_bboxlog()
            if box_type == 'fill-text'
        ]
        bounds = boxes[0]
        for other in boxes[1:]:
            bounds |= other
        for drawing in page.get_drawings():
            rect = drawing['rect']
            if rect.width < 610 or rect.height < 790:
                bounds |= rect
        clip = fitz.Rect(L, max(T, bounds.y0 - 3), R, min(B, bounds.y1 + 3))
        assert clip.height < 610, (number, clip)
        png = MEDIA / f'Figure_S{number}.png'
        page.get_pixmap(
            matrix=fitz.Matrix(300 / 72, 300 / 72),
            clip=clip, alpha=False,
        ).save(png)
        assert tuple(page.rect) == (0, 0, 612, 792) and path.stat().st_size < 10_000_000
    registry.append({
        'number': number,
        'figure': f'Figure S{number}',
        'file': path.name,
        'kind': kind,
        'title': title,
        'caption': caption,
        'word_png': str(png),
        'word_width_in': clip.width / 72,
        'word_height_in': clip.height / 72,
        'layout_frame_pt': [L, T, R, B],
        **extra,
    })
    (FIG / 'figure_registry.json').write_text(
        json.dumps(registry, ensure_ascii=False, indent=2), encoding='utf8'
    )
    print('FIGURE_DONE', number, title, flush=True)


# ------------------------------------------------------------------------------
# Figure S1: every dataset, four-chain numerical diagnostics
# ------------------------------------------------------------------------------

fig = newfig()
axes = [fig.add_axes([.1, .1, .8, .8]) for _ in range(3)]
rows = np.arange(len(summary))

axes[0].scatter(summary.rhat_max, rows, s=16, color=COLORS[1])
axes[0].axvline(1.01, color='black', ls='--', lw=.75)
axes[0].set(
    xlabel='Maximum rank-normalized split R-hat',
    yticks=rows, yticklabels=summary.display_name, xlim=(.999, 1.012),
)

axes[1].scatter(summary.ess_bulk_min, rows - .12, s=14,
               color=COLORS[1], label='Bulk ESS')
axes[1].scatter(summary.ess_tail_min, rows + .12, s=14,
               color=COLORS[2], marker='s', label='Tail ESS')
axes[1].axvline(400, color='black', ls='--', lw=.75)
axes[1].set(
    xscale='log', xlabel='Minimum multi-chain ESS',
    yticks=rows, yticklabels=summary.display_name,
)
axes[1].legend(frameon=False, loc='upper right')

axes[2].scatter(summary.single_chain_ess_mean_min, rows, s=16, color=COLORS[3])
axes[2].axvline(200, color='black', ls='--', lw=.75)
axes[2].set(
    xscale='log', xlabel='Minimum single-chain mean ESS',
    yticks=rows, yticklabels=summary.display_name,
)

for ax, values in [(axes[1], [400, 1000, 4000, 10000]), (axes[2], [200, 400, 800])]:
    ax.set_xticks(values)
    ax.set_xticklabels([str(v) for v in values])
    ax.xaxis.set_minor_formatter(NullFormatter())

cells = [
    (L, T, R, T + 165),
    (COLS[0][0], T + 190, COLS[0][1], T + 425),
    (COLS[1][0], T + 190, COLS[1][1], T + 425),
]
for index, (ax, cell) in enumerate(zip(axes, cells)):
    ax.invert_yaxis()
    ax.grid(axis='x', lw=.5, color='#E5E5E5')
    for tick, name in zip(ax.get_yticklabels(), summary.display_name):
        tick.set_fontstyle('normal' if name == 'Major tick' else 'italic')
    fit(fig, ax, cell, index)

save(
    fig,
    'Four-chain convergence diagnostics for the nine datasets',
    'Each group comprises the single chain of the original figure plus the three added '
    'chains 1-3, 36 chains in total. (A) Maximum rank-normalized split R-hat across all '
    'node ages and mu, sigma2 and lnL of each group; dashed line at 1.01. (B) Minimum joint '
    'bulk and tail ESS; dashed line at 400. (C) Minimum single-chain mean ESS; dashed line '
    'at 200. ESS, effective sample size. Complete values are kept as reproducible data.',
    'overview',
)


# ------------------------------------------------------------------------------
# Figures S2-S4: all 229 nodes, four chains, on a comfortable 8.5-pt row pitch
# ------------------------------------------------------------------------------

groups = [
    ['01-tick', 'Alectorobius', 'Amblyomma'],
    ['Dermacentor', 'Haemaphysalis1', 'Hyalomma'],
    ['Ixodes', 'Ornithodoros', 'Rhipicephalus'],
]

for group in groups:
    fig = newfig()
    legend(fig, markers=True)

    # Balance the panels across the two columns without splitting a dataset
    # over two pages.
    heights = {
        ds: len(chains[
            (chains.dataset == ds) & (chains.chain_index == 0)
            & chains.parameter.str.startswith('t_n')
        ]) * 8.5 + 50
        for ds in group
    }
    colload = [0., 0.]
    allocation = {}
    for ds in sorted(group, key=lambda k: heights[k], reverse=True):
        column = int(np.argmin(colload))
        allocation[ds] = (column, colload[column])
        colload[column] += heights[ds] + 15
    assert max(colload) < 590, colload
    if allocation[group[0]][0] == 1:
        allocation = {k: (1 - c, y) for k, (c, y) in allocation.items()}
    ordered = sorted(group, key=lambda k: (allocation[k][1], allocation[k][0]))

    for index, ds in enumerate(ordered):
        df = chains[
            (chains.dataset == ds) & chains.parameter.str.startswith('t_n')
        ]
        params = df[df.chain_index == 0].parameter.tolist()
        display = df.display_name.iloc[0]
        ax = fig.add_axes([.1, .1, .8, .8])
        for chain in range(4):
            q = df[df.chain_index == chain].set_index('parameter').loc[params]
            ys = np.arange(len(params)) + (chain - 1.5) * .19
            ax.hlines(ys, q.HPD95_low_Ma, q.HPD95_high_Ma,
                      color=COLORS[chain], lw=.65)
            ax.scatter(q.age_mean_Ma, ys, s=5.0,
                       color=COLORS[chain], marker=MARKERS[chain], zorder=3)
        ax.set(yticks=np.arange(len(params)), yticklabels=params,
               xlabel='Node age (Ma)')
        ax.invert_yaxis()
        ax.grid(axis='x', color='#E5E5E5', lw=.5)
        column, offset = allocation[ds]
        cell = (COLS[column][0], T + 27 + offset,
                COLS[column][1], T + 27 + offset + heights[ds])
        fit(fig, ax, cell, index, display, display != 'Major tick')

    labels = [summary.set_index('dataset').loc[k, 'display_name'] for k in ordered]
    save(
        fig,
        ', '.join(labels) + ': four-chain node ages and 95% HPD',
        '; '.join(
            f'({string.ascii_uppercase[i]}) {label}'
            for i, label in enumerate(labels)
        )
        + '. All internal nodes of each group are shown; black is the chain of the original '
        'figure and blue, orange and purple are the added chains 1-3. Dots are per-chain '
        'posterior means and the horizontal bars are 95% highest posterior density intervals; '
        'node indices and the corresponding descendant species sets are kept in the '
        'per-node reproducible data. Ma, million years; HPD, highest posterior density. '
        'The four chains are used for diagnostic comparison; final ages were not re-estimated '
        'from the four chains combined.',
        'all_node_hpd',
        datasets=ordered,
    )


# ------------------------------------------------------------------------------
# Figure S5: original-figure result versus the three already added runs
# ------------------------------------------------------------------------------

fig = newfig()
top = T
row_height = 104.
vgap = 15

for index, (ds, q) in enumerate(comparison.groupby('dataset', sort=False)):
    ax = fig.add_axes([.1, .1, .8, .8])
    limit = max(q.original_mean_Ma.max(), q.new_three_chain_mean_Ma.max()) * 1.06
    ax.plot([0, limit], [0, limit], color='black', ls='--', lw=.75)
    ax.scatter(q.original_mean_Ma, q.new_three_chain_mean_Ma,
               s=9, color=COLORS[1])
    ax.set(
        xlim=(0, limit), ylim=(0, limit),
        xlabel='Original-chain mean age (Ma)',
        ylabel='Added three-chain mean\nage (Ma)',
    )
    ax.xaxis.set_major_locator(MaxNLocator(3))
    ax.yaxis.set_major_locator(MaxNLocator(3))
    column = index % 2
    row = index // 2
    cell = (
        COLS[column][0], top + row * (row_height + vgap),
        COLS[column][1], top + row * (row_height + vgap) + row_height,
    )
    fit(fig, ax, cell, index, q.display_name.iloc[0], ds != '01-tick')

save(
    fig,
    'Node-age agreement between the original chain and the three added chains',
    '(A-I) Major tick, Alectorobius, Amblyomma, Dermacentor, Haemaphysalis, Hyalomma, '
    'Ixodes, Ornithodoros and Rhipicephalus, in this order. Each point is an internal '
    'node matched by its set of descendant species; the dashed line marks equal means. '
    'The x axis is the mean of the original-figure chain and the y axis the combined '
    'posterior mean of the three added chains, both in Ma. This panel is a repeatability '
    'check and should be read together with Figure S1 and the four-chain diagnostics. '
    'Per-node ages and 95% HPD are kept as reproducible data.',
    'original_vs_added',
)


# ------------------------------------------------------------------------------
# Figures S6-S7: one conservative, predetermined parameter per dataset, namely
# the parameter with the maximum four-chain R-hat
# ------------------------------------------------------------------------------

choices = []
for ds in manifest['datasets']:
    q = diag[diag.dataset == ds['dataset']]
    worst = q.loc[q.rhat_rank_split.idxmax()]
    choices.append((ds, worst.parameter, float(worst.rhat_rank_split)))

for start, count in [(0, 5), (5, 4)]:
    fig = newfig()
    legend(fig)
    details = []
    for row, (ds, param, _rhat) in enumerate(choices[start:start + count]):
        name = ds['dataset']
        cache = np.load(OUT / f'{name}_plot.npz')
        scale = 100 if param.startswith('t_n') else 1
        axes = [fig.add_axes([.1, .1, .8, .8]) for _ in range(2)]
        for chain in range(4):
            axes[0].plot(
                cache['iteration'][chain], cache['trace'][chain],
                color=COLORS[chain], ls=STYLES[chain], lw=.5, alpha=.7,
            )
            axes[1].stairs(
                cache['density'][chain], cache['edges'],
                color=COLORS[chain], ls=STYLES[chain], lw=.75,
            )
        axes[0].set(
            xlabel='Post-burn-in iteration (millions)',
            ylabel=param + (' (Ma)' if scale == 100 else ''),
        )
        axes[0].ticklabel_format(axis='y', style='plain', useOffset=False)
        axes[1].set(
            xlabel=param + (' (Ma)' if scale == 100 else ''),
            ylabel='Posterior density',
        )
        for column, ax in enumerate(axes):
            height = 105 if count == 5 else 120
            y_start = T + 27 + row * (height + 14)
            cell = (COLS[column][0], y_start, COLS[column][1], y_start + height)
            fit(
                fig, ax, cell, row * 2 + column,
                ds['display_name'] if column == 0 else None,
                ds['display_name'] != 'Major tick',
            )
        details.append(
            f'({string.ascii_uppercase[row * 2]}, '
            f'{string.ascii_uppercase[row * 2 + 1]}) {ds["display_name"]}: {param}'
        )
        cache.close()
    save(
        fig,
        'Four-chain traces and posterior densities for representative parameters',
        '; '.join(details)
        + '. The parameter with the highest four-chain R-hat is used for each group; '
        'left, trace; right, posterior density. Black, blue, orange and purple are the '
        'chain of the original figure and the added chains 1-3. At most 2400 evenly spaced '
        'points per chain are drawn, while the densities and diagnostics use all stored '
        'samples. Ages are given in Ma and the other parameters keep their original units. '
        'All parameters enter the numerical diagnostics.',
        'representative_traces',
        datasets=[d['dataset'] for d, p, r in choices[start:start + count]],
    )

assert len(registry) == 7
print('COMPACT_FIGURES_COMPLETE', len(registry), flush=True)
