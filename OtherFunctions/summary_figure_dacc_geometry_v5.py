import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import ListedColormap

# Optional neuroimaging imports
try:
    import nibabel as nib
    from nilearn import plotting, image
    HAVE_NILEARN = True
except Exception:
    HAVE_NILEARN = False


# =========================
# User-adjustable inputs
# =========================
MASK_FILE = "/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/Manuscript/Revision5/SA/SuppFigures/Pipeline/m01_mu_p_FWEc227_mask.nii"
OUT_PNG = "/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/Manuscript/Revision5/SA/SuppFigures/Pipeline/dacc_geometry_summary_science_style_final_v8.png"
OUT_PDF = "/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/Manuscript/Revision5/SA/SuppFigures/Pipeline/dacc_geometry_summary_science_style_final_v8.pdf"
OUT_SVG = "/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/Manuscript/Revision5/SA/SuppFigures/Pipeline/dacc_geometry_summary_science_style_final_v8.svg"
DPI = 300

THETA_DEG = 84
SCALE = 0.9
R_BEHAV = 0.494
P_BEHAV = 0.006
N_SUBJ = 29

# Optional: replace with your real data
ORTHO_X = None
ACCURACY_Y = None


# =========================
# Style
# =========================
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 14,
    "axes.linewidth": 1.0,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none"
})

COLORS = {
    "error": "#3E7CB1",
    "difficulty": "#D9911A",
    "readout": "#6B59B3",
    "switch": "#2E9C57",
    "dacc_magenta": "#D81B60",
    "gray1": "#222222",
    "gray2": "#666666",
    "gray3": "#B8B8B8",
    "gray4": "#E8E8E8",
    "bg": "#FAFAF8",
    "white": "#FFFFFF"
}


# =========================
# Helpers
# =========================
def rotation_matrix(theta_deg: float) -> np.ndarray:
    theta = np.deg2rad(theta_deg)
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s], [s, c]])


def style_panel(ax):
    ax.set_facecolor(COLORS["bg"])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.patch.set_alpha(0.35)

    for s in ax.spines.values():
        s.set_visible(False)


def add_panel_header(ax, label, title, subtitle=None):
    ax.text(0.00, 1.03, label,
            transform=ax.transAxes,
            ha='left', va='bottom',
            fontsize=8, fontweight='bold',
            color=COLORS["gray1"])

    ax.text(0.12, 1.01, title,
            transform=ax.transAxes,
            ha='left', va='bottom',
            fontsize=8, fontweight='bold',
            color=COLORS["gray1"])

    if subtitle is not None:
        ax.text(0.05, 0.972, subtitle,
                transform=ax.transAxes,
                ha='left', va='top',
                fontsize=6,
                color=COLORS["gray2"])


def add_between_panel_arrow(fig, x0, x1, y):
    arr = FancyArrowPatch(
        (x0, y), (x1, y),
        transform=fig.transFigure,
        arrowstyle='-|>',
        mutation_scale=18,
        lw=1.8,
        color=COLORS["gray1"]
    )
    fig.add_artist(arr)


def prettify_float_p(p):
    if p < 0.001:
        return f"{p:.1e}"
    return f"{p:.3f}"

from matplotlib.patches import Circle, Rectangle, FancyBboxPatch

def add_task_icon_biorender(ax, x=0.16, y=0.84, s=0.18,
                            line="#7F8C8D",
                            body="#DDE6F2",
                            accent="#AFC8E8",
                            screen="#F3F6FA",
                            alpha=0.95):
    """
    BioRender-like task icon drawn in axes-fraction coordinates.

    Parameters
    ----------
    ax : matplotlib axis
    x, y : anchor position in ax.transAxes coordinates
    s : overall size
    """
    tr = ax.transAxes
    lw = 1.0

    # ----- chair -----
    ax.add_patch(Rectangle(
        (x - 0.10*s, y - 0.28*s), 0.03*s, 0.16*s,
        transform=tr, facecolor=accent, edgecolor=line, lw=lw, alpha=alpha
    ))
    ax.add_patch(Rectangle(
        (x - 0.10*s, y - 0.12*s), 0.10*s, 0.025*s,
        transform=tr, facecolor=accent, edgecolor=line, lw=lw, alpha=alpha
    ))

    # ----- desk -----
    ax.add_patch(Rectangle(
        (x + 0.10*s, y - 0.14*s), 0.22*s, 0.025*s,
        transform=tr, facecolor=accent, edgecolor=line, lw=lw, alpha=alpha
    ))
    ax.add_patch(Rectangle(
        (x + 0.28*s, y - 0.14*s), 0.02*s, 0.18*s,
        transform=tr, facecolor=accent, edgecolor=line, lw=lw, alpha=alpha
    ))

    # ----- monitor -----
    ax.add_patch(FancyBboxPatch(
        (x + 0.17*s, y - 0.07*s), 0.12*s, 0.08*s,
        boxstyle="round,pad=0.002,rounding_size=0.01",
        transform=tr, facecolor=screen, edgecolor=line, lw=lw, alpha=alpha
    ))
    ax.add_patch(Rectangle(
        (x + 0.222*s, y - 0.10*s), 0.018*s, 0.03*s,
        transform=tr, facecolor=accent, edgecolor=line, lw=lw, alpha=alpha
    ))

    # ----- head -----
    ax.add_patch(Circle(
        (x, y), 0.045*s,
        transform=tr, facecolor=body, edgecolor=line, lw=lw, alpha=alpha
    ))

    # ----- torso -----
    ax.add_patch(FancyBboxPatch(
        (x - 0.045*s, y - 0.15*s), 0.09*s, 0.09*s,
        boxstyle="round,pad=0.002,rounding_size=0.02",
        transform=tr, facecolor=body, edgecolor=line, lw=lw, alpha=alpha
    ))

    # ----- arm to desk -----
    ax.plot(
        [x + 0.03*s, x + 0.11*s],
        [y - 0.11*s, y - 0.13*s],
        color=line, lw=1.2, alpha=alpha, transform=tr, clip_on=False,
        solid_capstyle='round'
    )

    # ----- thigh -----
    ax.plot(
        [x + 0.02*s, x + 0.10*s],
        [y - 0.15*s, y - 0.22*s],
        color=line, lw=1.3, alpha=alpha, transform=tr, clip_on=False,
        solid_capstyle='round'
    )

    # ----- lower leg -----
    ax.plot(
        [x + 0.10*s, x + 0.05*s],
        [y - 0.22*s, y - 0.30*s],
        color=line, lw=1.3, alpha=alpha, transform=tr, clip_on=False,
        solid_capstyle='round'
    )

    # ----- back hint -----
    ax.plot(
        [x - 0.01*s, x - 0.06*s],
        [y - 0.08*s, y - 0.14*s],
        color=line, lw=1.1, alpha=alpha, transform=tr, clip_on=False,
        solid_capstyle='round'
    )

# =========================
# Panel A
# =========================
def add_brain_panel(fig, panel_pos, mask_file):
    axA = fig.add_axes(panel_pos)
    style_panel(axA)
    add_panel_header(
        axA,
        "A",
        "dACC locus",
        "Confidence-related dACC region"
    )

    x, y, w, h = panel_pos

    # Main brain axis
    brain_ax = fig.add_axes([x + 0.06 * w, y + 0.17 * h, 0.88 * w, 0.64 * h])
    brain_ax.set_facecolor("white")
    brain_ax.set_xticks([])
    brain_ax.set_yticks([])
    for s in brain_ax.spines.values():
        s.set_visible(False)

    if HAVE_NILEARN and os.path.exists(mask_file):
        try:
            mask_img = nib.load(mask_file)

            # ensure mask is binary 0/1
            mask_data = mask_img.get_fdata()
            mask_bin = (mask_data > 0).astype(np.uint8)
            mask_bin_img = nib.Nifti1Image(mask_bin, mask_img.affine, mask_img.header)

            # color map: whole mask rendered as crimson
            dacc_cmap = ListedColormap([COLORS["dacc_magenta"]])

            # glass brain, white/transparent background feel
            display = plotting.plot_glass_brain(
                stat_map_img=None,
                display_mode='r',
                axes=brain_ax,
                black_bg=False,
                colorbar=False,
                annotate=False,
                plot_abs=False,
                alpha=0.35
            )
            # from nilearn import datasets
            # template = datasets.load_mni152_template()
            # display = plotting.plot_roi(
            #     roi_img=mask_bin_img,
            #     bg_img=template,
            #     display_mode='x',
            #     cut_coords=[8],          # 可调
            #     axes=brain_ax,
            #     cmap=dacc_cmap,
            #     black_bg=False,
            #     alpha=0.95,
            #     annotate=False,
            #     draw_cross=False,
            #     dim=0.4
            # )
            # overlay the binary mask as full crimson
            display.add_overlay(
                mask_bin_img,
                cmap=dacc_cmap,
                threshold=0.5,
                alpha=1.0
            )

            # crisp boundary
            try:
                display.add_contours(
                    mask_bin_img,
                    levels=[0.5],
                    colors=[COLORS["dacc_magenta"]],
                    linewidths=1.4
                )
            except Exception:
                pass

        except Exception as e:
            brain_ax.text(0.5, 0.58, "Brain rendering failed",
                          ha='center', va='center',
                          fontsize=8, fontweight='bold', color=COLORS["gray1"])
            brain_ax.text(0.5, 0.40, str(e),
                          ha='center', va='center',
                          fontsize=8, color=COLORS["gray2"], wrap=True)
    else:
        brain_ax.add_patch(plt.Rectangle((0, 0), 1, 1,
                                         facecolor="white",
                                         edgecolor=COLORS["gray4"],
                                         linewidth=0.8))
        brain_ax.text(0.5, 0.58, "Brain rendering unavailable",
                      ha='center', va='center',
                      fontsize=8, fontweight='bold', color=COLORS["gray1"])
        brain_ax.text(0.5, 0.40, "Require nibabel + nilearn + valid NIfTI mask",
                      ha='center', va='center',
                      fontsize=8, color=COLORS["gray2"])

    axA.text(
        0.50, 0.08,
        "Confidence-related dACC locus",
        transform=axA.transAxes,
        ha='center', va='center',
        fontsize=8, color=COLORS["dacc_magenta"], fontweight='bold'
    )

    return axA


# =========================
# Panel B
# =========================
def add_geometry_panel(fig, panel_pos):
    axB = fig.add_axes(panel_pos)
    style_panel(axB)

    col_err_axis = "#D8A031"
    col_diff_axis = "#6FAF63"
    col_switch_axis = COLORS["dacc_magenta"]

    add_panel_header(
        axB,
        "B",
        "Neural geometry",
        "Switch evidence is read out along a third axis \nspanning error and difficulty"
    )

    x, y, w, h = panel_pos

    geom_ax = fig.add_axes([x + 0.08 * w, y + 0.04 * h, 0.82 * w, 0.74 * h])
    style_panel(geom_ax)
    geom_ax.set_xlim(0, 1)
    geom_ax.set_ylim(0, 1)

    # common origin
    x0, y0 = 0.16, 0.20

    # ---------------------------------
    # Base axes: error and difficulty
    # ---------------------------------
    geom_ax.annotate(
        "", xy=(0.92, y0), xytext=(x0, y0),
        arrowprops=dict(arrowstyle="-|>", lw=2.7, color=col_err_axis)
    )
    geom_ax.annotate(
        "", xy=(x0, 0.90), xytext=(x0, y0),
        arrowprops=dict(arrowstyle="-|>", lw=2.7, color=col_diff_axis)
    )

    geom_ax.text(
        0.5, 0.05, "Number of errors",
        ha='center', va='center',
        fontsize=6, color=col_err_axis, fontweight='bold'
    )
    geom_ax.text(
        0.07, 0.58, "Difficulty",
        ha='center', va='center',
        fontsize=6, color=col_diff_axis, fontweight='bold',
        rotation=90
    )

    # ---------------------------------
    # Grid in error-difficulty plane
    # ---------------------------------
    y_levels = [0.34, 0.56, 0.78]
    x_levels = [0.32, 0.56, 0.80]

    for yy, lab in zip(y_levels, ["1", "2", "3"]):
        geom_ax.plot([x0, 0.92], [yy, yy],
                     ls=(0, (4, 3)), lw=1.7,
                     color=col_diff_axis, alpha=0.45, zorder=1)
        geom_ax.text(0.135, yy,
                     lab,
                     ha='right', va='center',
                     fontsize=6.0, color=col_diff_axis)

    for xx, lab in zip(x_levels, ["0", "1", "2"]):
        geom_ax.plot([xx, xx], [y0, 0.90],
                     ls=(0, (4, 3)), lw=1.7,
                     color=col_err_axis, alpha=0.45, zorder=1)
        geom_ax.text(xx, 0.145,
                     lab,
                     ha='center', va='top',
                     fontsize=6.0, color=col_err_axis)

    # ---------------------------------
    # Projection of switch-evidence axis
    # onto the error-difficulty plane
    # ---------------------------------
    proj_x = np.array([0.25, 0.38, 0.51, 0.64, 0.77, 0.84])
    proj_y = np.array([0.25, 0.38, 0.51, 0.64, 0.76, 0.82])

    geom_ax.plot(
        proj_x, proj_y,
        color=col_switch_axis, alpha=0.45,
        lw=2.3,
        ls=(0, (5, 3)),
        solid_capstyle='round',
        zorder=3,
    )

    geom_ax.scatter(proj_x[0], proj_y[0], s=28, color=col_switch_axis, alpha=0.45, zorder=4)
    geom_ax.scatter(proj_x[-1], proj_y[-1], s=28, color=col_switch_axis, alpha=0.45, zorder=4)

    geom_ax.text(
        0.56, 0.45,
        "Projection onto\nerror-difficulty plane",
        fontsize=5.6,
        color=col_switch_axis, alpha=0.65,
        fontweight='bold',
        ha='left',
        va='bottom',
        linespacing=1.35
    )

    # ---------------------------------
    # True switch-evidence axis
    # starts from SAME origin
    # ---------------------------------
    sw_end_x, sw_end_y = 0.73, 0.92

    geom_ax.annotate(
        "",
        xy=(sw_end_x, sw_end_y), xytext=(x0, y0),
        arrowprops=dict(arrowstyle="-|>", lw=3.0, color=col_switch_axis)
    )

    #geom_ax.scatter(sw_end_x, sw_end_y, s=34, color=col_switch_axis, zorder=5)

    geom_ax.text(
        0.55, 0.95,
        "Switch evidence axis",
        fontsize=6,
        color=col_switch_axis,
        fontweight='bold',
        ha='left',
        va='bottom'
    )

    # ---------------------------------
    # Projection helper: from true-axis endpoint
    # to projected endpoint
    # ---------------------------------
    geom_ax.plot(
        [sw_end_x-0.015, proj_x[-1]], [sw_end_y-0.015, proj_y[-1]],
        color=COLORS["gray3"], lw=2, ls=(0, (2, 2)), zorder=3
    )

    # ---------------------------------
    # Endpoint annotations
    # ---------------------------------
    geom_ax.text(
        0.865, 0.845,
        "2 errors\nDifficulty 3",
        ha='left', va='center',
        fontsize=6.0,
        color=COLORS["gray1"],
        linespacing=1.35
    )

    geom_ax.text(
        0.170, 0.200,
        "0 errors\nDifficulty 1",
        ha='right', va='top',
        fontsize=6,
        color=COLORS["gray2"],
        linespacing=1.35
    )

    geom_ax.text(
        0.12, 1.05,
        "Near-orthogonal geometry",
        ha='left', va='bottom',
        fontsize=8, color=COLORS["gray1"], fontweight='bold'
    )

    return axB


# =========================
# Panel C
# =========================
def add_readout_panel(fig, panel_pos):
    axC = fig.add_axes(panel_pos)
    style_panel(axC)
    add_panel_header(
        axC,
        "C",
        "Readout",
        "Confidence drives switch probability \nvia threshold crossing"
    )

    # -------- left: latent switch confidence axis --------
    x0, x1, y0 = 0.14, 0.55, 0.50
    axC.plot([x0, x1], [y0, y0], lw=2.4, color=COLORS["readout"])
    axC.annotate(
        "",
        xy=(x1, y0), xytext=(x0, y0),
        arrowprops=dict(arrowstyle="-|>", lw=2.4, color=COLORS["readout"])
    )

    axC.text(
        (x0 + x1) / 2, y0 + 0.12,
        "Latent switch confidence",
        ha='center', va='bottom',
        fontsize=5.2, color=COLORS["readout"], fontweight='bold'
    )

    axC.text(
        x0, y0 - 0.08, "low",
        ha='center', va='top',
        fontsize=8, color=COLORS["gray2"]
    )
    axC.text(
        x1, y0 - 0.08, "high",
        ha='center', va='top',
        fontsize=8, color=COLORS["gray2"]
    )

    dots = [
        ("hard +\nfew errors", 0.10, COLORS["gray3"]),
        ("mid", 0.42, COLORS["gray2"]),
        ("easy +\nrepeated errors", 0.86, COLORS["readout"])
    ]

    for txt, frac, col in dots:
        xp = x0 + frac * (x1 - x0)
        axC.scatter(
            xp, y0, s=62,
            color=col, edgecolor='white',
            linewidth=0.8, zorder=3
        )
        axC.text(
            xp, y0 + 0.06, txt,
            ha='center', va='bottom',
            fontsize=4, color=col, linespacing=1.35
        )

    axC.annotate(
        "", xy=(0.68, y0), xytext=(0.59, y0),
        arrowprops=dict(arrowstyle="-|>", lw=1.6, color=COLORS["gray2"])
    )

    # -------- right: threshold-crossing schematic --------
    dist_ax_x = np.linspace(0.74, 0.94, 500)
    mu = 0.84
    sigma = 0.035
    baseline = 0.28
    amp = 0.30

    # Gaussian-like confidence distribution
    yy = baseline + amp * np.exp(-0.5 * ((dist_ax_x - mu) / sigma) ** 2)

    # threshold
    theta_x = 0.875

    axC.plot(dist_ax_x, yy, color=COLORS["dacc_magenta"], alpha=0.6, lw=2.4)

    # Fill tail above threshold = P(switch)
    mask = dist_ax_x >= theta_x
    axC.fill_between(
        dist_ax_x[mask], baseline, yy[mask],
        color=COLORS["dacc_magenta"], alpha=0.28, zorder=1
    )

    axC.plot([0.74, 0.94], [baseline, baseline], lw=1.2, color=COLORS["gray3"])
    axC.plot([theta_x, theta_x], [baseline, baseline + amp * 1.06],
             lw=1.4, ls=(0, (3, 2)), color=COLORS["gray2"])

    axC.text(
        theta_x, baseline + amp * 1.11,
        r"$\theta$",
        ha='center', va='bottom',
        fontsize=7.5, color=COLORS["gray1"], fontweight='bold'
    )

    axC.text(
        0.875, 0.67, "P(switch)",
        ha='center', va='bottom',
        fontsize=6.0, color=COLORS["dacc_magenta"], fontweight='bold'
    )
    #axC.text(
    #    0.835, 0.245, "confidence",
    #    ha='center', va='top',
    #    fontsize=6, color=COLORS["gray2"]
    #)
    axC.text(
        0.905, 0.245, "threshold",
        ha='center', va='top',
        fontsize=5.4, color=COLORS["gray2"]
    )

    # subtle arrow from distribution body to tail
    axC.annotate(
        "", xy=(0.865, 0.6), xytext=(0.95, 0.6),
        arrowprops=dict(arrowstyle="<|-", lw=1.2, color=COLORS["dacc_magenta"])
    )

    axC.text(
        0.02, 0.02,
        "Readout aligned with the error-\ndifficulty span yields graded confidence",
        transform=axC.transAxes,
        ha='left', va='bottom',
        fontsize=6, color=COLORS["gray2"], linespacing=1.35
    )

    axC.set_xlim(0, 1)
    axC.set_ylim(0, 1)
    return axC


# =========================
# Panel D
# =========================
def generate_behavior_data():
    if ORTHO_X is not None and ACCURACY_Y is not None:
        x = np.asarray(ORTHO_X).reshape(-1)
        y = np.asarray(ACCURACY_Y).reshape(-1)
        return x, y

    rng = np.random.default_rng(7)
    x = np.linspace(-0.016, 0.018, N_SUBJ)
    y = 0.492 + 2.75 * x + rng.normal(0, 0.010, size=x.size)
    return x, y


def add_behavior_panel(fig, panel_pos):
    axD = fig.add_axes(panel_pos)
    axD.set_facecolor(COLORS["bg"])
    axD.patch.set_alpha(0.35)

    add_panel_header(
        axD,
        "D",
        "Behavior",
        "Orthogonal geometry predicts \nbetter inference"
    )

    x, y = generate_behavior_data()
    coef = np.polyfit(x, y, 1)
    yfit = np.polyval(coef, x)

    px, py, pw, ph = panel_pos
    plot_ax = fig.add_axes([px + 0.18 * pw, py + 0.25 * ph, 0.74 * pw, 0.50 * ph])
    plot_ax.set_facecolor(COLORS["bg"])

    plot_ax.scatter(
        x, y, s=42,
        color=COLORS["error"],
        alpha=0.8,
        edgecolor='white',
        linewidth=0.6,
        zorder=3
    )
    plot_ax.plot(x, yfit, lw=2.4, color=COLORS["error"], alpha=0.8, zorder=2)

    plot_ax.spines['top'].set_visible(False)
    plot_ax.spines['right'].set_visible(False)
    plot_ax.spines['left'].set_color(COLORS["gray2"])
    plot_ax.spines['bottom'].set_color(COLORS["gray2"])

    plot_ax.set_xticks([])
    plot_ax.set_yticks([])

    plot_ax.set_xlabel("Orthogonality index", fontsize=6, color=COLORS["gray1"], labelpad=6)
    plot_ax.set_ylabel("Inference accuracy", fontsize=6, color=COLORS["gray1"], labelpad=6)

    xr = np.max(x) - np.min(x)
    yr = np.max(y) - np.min(y)
    plot_ax.set_xlim(np.min(x) - 0.08 * xr, np.max(x) + 0.08 * xr)
    plot_ax.set_ylim(np.min(y) - 0.08 * yr, np.max(y) + 0.12 * yr)

    add_task_icon_biorender(
    plot_ax,
    x=0.18, y=0.9, s=1.1,
    line="#8C96A3",
    body="#E8EEF7",
    accent="#C9D9EE",
    screen="#F7F9FC",
    alpha=0.95
    )

    axD.set_xticks([])
    axD.set_yticks([])

    for s in axD.spines.values():
        s.set_visible(True) #False
        s.set_color(COLORS["gray4"])
        s.set_linewidth(1.0)

    return axD


# =========================
# Main
# =========================
def main():
    fig = plt.figure(figsize=(9.4, 3.0), facecolor="white")

    panelA = [0.015, 0.05, 0.24, 0.76]
    panelB = [0.295, 0.05, 0.24, 0.76]
    panelC = [0.575, 0.05, 0.18, 0.76]
    panelD = [0.795, 0.05, 0.185, 0.76]

    add_brain_panel(fig, panelA, MASK_FILE)
    add_geometry_panel(fig, panelB)
    add_readout_panel(fig, panelC)
    add_behavior_panel(fig, panelD)

    y_arrow = 0.50
    add_between_panel_arrow(fig, panelA[0] + panelA[2] + 0.006, panelB[0] - 0.008, y_arrow)
    add_between_panel_arrow(fig, panelB[0] + panelB[2] + 0.006, panelC[0] - 0.008, y_arrow)
    add_between_panel_arrow(fig, panelC[0] + panelC[2] + 0.006, panelD[0] - 0.008, y_arrow)

    fig.text(
        0.5, 0.92,
        "dACC geometry links error and difficulty coding to confidence readout and behavior",
        ha='center', va='center',
        fontsize=12, fontweight='bold', color=COLORS["gray1"]
    )

    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    #os.makedirs(os.path.dirname(OUT_PDF), exist_ok=True)
    #os.makedirs(os.path.dirname(OUT_SVG), exist_ok=True)

    plt.savefig(OUT_PNG, dpi=DPI, facecolor='white')
    #plt.savefig(OUT_PDF, dpi=DPI, facecolor='white')
    #plt.savefig(OUT_SVG, facecolor='white', format='svg')

    print(f"Saved: {OUT_PNG}")
    #print(f"Saved: {OUT_PDF}")
    #print(f"Saved: {OUT_SVG}")


if __name__ == "__main__":
    main()