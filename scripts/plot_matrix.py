from pathlib import Path
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter
STAGE = Path(os.environ["REPRO_WORK"]) / "BMC geriatrics/V1/figures/panels"
PANEL_B = "Figure3B_matched_case_estimate_matrix"
HRS = "#D55E00"; CHARLS = "#0072B2"; INK = "#242628"; NEUTRAL = "#8B9095"
def draw_panel(rows):
    STAGE.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(180 / 25.4, 90 / 25.4))
    fig.text(.025, .965, "B", fontweight="bold", fontsize=11, va="top", color=INK)
    fig.text(.055, .965, "Blood-marker adjustment within matched complete cases",
             fontsize=9, va="top", color=INK)
    canvas = fig.add_axes([0, 0, 1, 1], zorder=0)
    canvas.set_xlim(0, 1)
    canvas.set_ylim(0, 1)
    canvas.axis("off")
    left, middle, right = .28, .64, .985
    centres = [.465, .825]
    for x, label in zip(centres, ["Combined eGFR\n+ clinical covariates", "+ blood markers\n(same complete cases)"]):
        fig.text(x, .84, label, ha="center", va="top", fontsize=9, color=INK)
    for y in [.75, .405, .075]:
        canvas.plot([.025, right], [y, y], color="#C6C9CB", lw=.7)
    for x in [left, middle]:
        canvas.plot([x, x], [.075, .75], color="#E2E4E5", lw=.6)
    for cohort_index, (name, colour, marker) in enumerate([("HRS", HRS, "o"), ("CHARLS", CHARLS, "s")]):
        y_top = .70 - .34 * cohort_index
        group = rows[cohort_index * 2:cohort_index * 2 + 2]
        fig.text(.035, y_top, name, ha="left", va="top", fontsize=10, fontweight="bold", color=colour)
        fig.text(.035, y_top - .10,
                 f"{int(group[0].n_people):,} participants\n{int(group[0].n_intervals):,} intervals\n"
                 f"{int(group[0].events):,} recoveries", ha="left", va="top", fontsize=8.5,
                 linespacing=1.45, color=INK)
        for model_index, row in enumerate(group):
            x_start = .32 + .36 * model_index
            ax = fig.add_axes([x_start, y_top - .16, .29, .08], zorder=1)
            ax.set_xscale("log")
            ax.set_xlim(.95, 1.30)
            ax.set_ylim(0, 1)
            ax.errorbar(row.odds_ratio, .5,
                        xerr=[[row.odds_ratio - row.ci_lower], [row.ci_upper - row.odds_ratio]],
                        fmt=marker, color=colour, markerfacecolor="white" if model_index == 0 else colour,
                        markeredgewidth=1.3, markersize=5.3, capsize=2.5, linewidth=1.4, zorder=4)
            ax.axvline(1, color=NEUTRAL, ls="--", linewidth=.8, zorder=2)
            ax.xaxis.set_major_locator(FixedLocator([1, 1.1, 1.2, 1.3]))
            ax.xaxis.set_major_formatter(FixedFormatter(["1.0", "1.1", "1.2", "1.3"]))
            ax.xaxis.set_minor_locator(FixedLocator([]))
            ax.tick_params(axis="x", length=2.5, color=NEUTRAL, pad=2,
                           labelbottom=cohort_index == 1, bottom=cohort_index == 1, labelsize=7.5)
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            fig.text(centres[model_index], y_top, f"OR {row.odds_ratio:.3f}",
                     ha="center", va="top", fontsize=10.5, fontweight="bold", color=colour)
            fig.text(centres[model_index], y_top - .24,
                     f"95% CI {row.ci_lower:.3f}\u2013{row.ci_upper:.3f}",
                     ha="center", va="top", fontsize=8.7, color=INK)
    fig.text(.645, .025, "Recovery OR (logarithmic scale)", ha="center", va="bottom", fontsize=8, color=INK)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bounds = fig.bbox
    for text in fig.findobj(lambda artist: isinstance(artist, plt.Text)):
        if not text.get_visible() or not text.get_text():
            continue
        box = text.get_window_extent(renderer)
        assert box.x0 >= bounds.x0 - 1 and box.y0 >= bounds.y0 - 1
        assert box.x1 <= bounds.x1 + 1 and box.y1 <= bounds.y1 + 1, text.get_text()
    for extension in ["pdf", "svg", "png"]:
        fig.savefig(STAGE / f"{PANEL_B}.{extension}", dpi=450, facecolor="white")
    plt.close(fig)


