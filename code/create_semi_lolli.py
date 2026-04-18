#########################################################################
#                                      ###     ###    #######   ###     #
#                                      ###     ###   ###        ###     #
#                                      ###     ###   ###        ###     #
#                                      ###     ###   ###        ###     #
#                                       #########     #######   #########
#                                                                       #
#                                                                       #
#                      HypoSubstrateTesting                             #
# Testing hypothesized neural substrates with lesions                   #
#                                                                       #
# This script reads a single CSV file produced by simulations.m and     #
# creates a semi-lollipop plot of significance ratios for each noise    #
# bin across the four base images per mode (e.g. tracts:true substrate) #
#                                                                       #
#                                                                       #
# last update: 2026/03/07.                                              #
#                                                                       #
#                                                                       #
#########################################################################




########################################
#                                      #
#             LOAD LIBRARIES           #
#                                      #
########################################


# ---- load libraries ---- #
import os, argparse, numpy as np, matplotlib.pyplot as plt, pandas as pd
from matplotlib.markers import MarkerStyle
from matplotlib.transforms import Affine2D
from matplotlib.font_manager import FontProperties


#####################################################
#                                                   #
#                DEFINE FUNCTIONS                   #
#                                                   #
#####################################################



def read_csv(csv_path, add_spacers=False):
    '''
    Read a SignificanceRatios.csv produced by sim_analysis.m.

    When *add_spacers* is True an empty spoke (label '', values NaN) is
    inserted between every pair of consecutive base-image groups.  This
    creates visual breathing room in circular charts so that the last
    bin of one base image does not sit right next to the first bin of
    the next.

    Returns:
        labels  - list of spoke labels (baseImage + bin range)
        slope   - array of pctSurviving_slopeComparison values
        covar   - array of pctSurviving_withVolumeCovariate values
    '''
    df = pd.read_csv(csv_path)

    # Clean base image names
    df['cleanBase'] = df['baseImage']

    seen_order = []
    _seen = set()
    for b in df['cleanBase']:
        if b not in _seen:
            _seen.add(b)
            seen_order.append(b)
    group_rank = {name: idx for idx, name in enumerate(seen_order)}
    df['_group_rank'] = df['cleanBase'].map(group_rank)
    df = df.sort_values(['_group_rank', 'noiseBinStart']).reset_index(drop=True)
    df = df.drop(columns=['_group_rank'])

    # Build spoke labels: "MDN [0.20-0.25]"
    labels = [
        f"{row['cleanBase']} [{row['noiseBinStart']:.2f}-{row['noiseBinEnd']:.2f}]"
        for _, row in df.iterrows()
    ]

    slope = df['pctSurviving_slopeComparison'].values
    covar = df['pctSurviving_withVolumeCovariate'].values

    if add_spacers:
        # Detect group boundaries (where base image name changes)
        bases = df['cleanBase'].values
        new_labels = []
        new_slope  = []
        new_covar  = []
        for k in range(len(labels)):
            new_labels.append(labels[k])
            new_slope.append(slope[k])
            new_covar.append(covar[k])
            # create block design
            if k < len(labels) - 1 and bases[k] != bases[k + 1]:
                new_labels.append('')
                new_slope.append(np.nan)
                new_covar.append(np.nan)

        if len(labels) > 0:
            new_labels.append('')
            new_slope.append(np.nan)
            new_covar.append(np.nan)

        labels = new_labels
        slope  = np.array(new_slope)
        covar  = np.array(new_covar)

    return labels, slope, covar


def plot_semi_lollipop_chart(data, metrics, title, colors=['darkorange', 'teal'],
                             ylim=None, semicircle_radius=10, font_prop=None,
                             mode='true', spoke_color='split'):
    '''
    Circular semi-lollipop chart.

    Each spoke has a single light-grey radial line:
      - solid from the origin to max(model0, model1)
      - dashed from max(model0, model1) to the outer rim (and text label)

    At each model value a filled half-circle is drawn on the
    corresponding side of the spoke:
      - model 0 (colors[0], e.g. darkorange): clockwise side
        (negative angular offset)
      - model 1 (colors[1], e.g. teal): counter-clockwise side
        (positive angular offset)

    The outer text label shows the full bin range, e.g. "0.30-0.35".

    :param data:              2-D numpy array (n_spokes x 2)
    :param metrics:           list of spoke labels
    :param title:             chart title
    :param colors:            list of two colors for model 0 and model 1
    :param ylim:              tuple (min, max) for radial axis limits
    :param semicircle_radius: marker size for the half-circle markers
                              (in points, passed to markersize)
    :param font_prop:         matplotlib FontProperties instance for spoke
                              labels, or None for the default font
    :param mode:              label for axis annotation (true/false positive)
    :param spoke_color:       color for the solid spoke segments. Use
                              'split' to color the segment from origin to
                              min(v0,v1) with the lower model's color and
                              the segment from min to max with the higher
                              model's color.  Any other string is treated
                              as a matplotlib color name and draws the
                              entire solid portion in that color.
    '''
    N = len(metrics)
    theta = np.linspace(0, 2 * np.pi, N, endpoint=False)

    fig, ax = plt.subplots(figsize=(12, 12),
                           subplot_kw={'projection': 'polar'})
    ax.set_title(title, y=1.15, fontsize=20)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rlabel_position(90)
    ax.spines['polar'].set_visible(False)
    ax.xaxis.grid(False)
    ax.yaxis.grid(False)

    if ylim is not None:
        ax.set_ylim(ylim)
        r_max = ylim[1]
    else:
        r_max = np.nanmax(data)


    def _nice_interval(x):
        '''Round *x* down to the nearest value in {1,2,5}*10^n.'''
        if x <= 0:
            return 1
        exp = np.floor(np.log10(x))
        base = 10 ** exp
        frac = x / base
        if frac >= 5:
            return int(5 * base)
        elif frac >= 2:
            return int(2 * base)
        else:
            return int(base)

    interval = _nice_interval(r_max / 5)
    ring_angles = np.linspace(0, 2 * np.pi, 200)

    # Find the first spacer spoke to use for axis label placement
    spacer_idx = None
    for j in range(N):
        if np.any(np.isnan(data[j])):
            spacer_idx = j
            break
    label_angle = theta[spacer_idx] if spacer_idx is not None else 0.0

    font_kw_ring = {'fontproperties': font_prop} if font_prop else {}

    # Rotation so text runs along the spacer spoke
    _p0 = ax.transData.transform((label_angle, 0))
    _p1 = ax.transData.transform((label_angle, r_max))
    _dx, _dy = _p1[0] - _p0[0], _p1[1] - _p0[1]
    spoke_rot = np.degrees(np.arctan2(_dy, _dx))

    if spoke_rot < -90 or spoke_rot > 90:
        spoke_rot += 180
        va_below = 'top'
        va_above = 'bottom'
    else:
        va_below = 'top'
        va_above = 'bottom'

    for k in range(1, 5):
        r = interval * k
        if r >= r_max:
            break
        ax.plot(ring_angles, np.full_like(ring_angles, r),
                color='lightgrey', linewidth=0.75, linestyle='dashed',
                zorder=0)
        ax.text(label_angle, r, str(int(r)),
                ha='center', va=va_below,
                fontsize=14, color='black',
                rotation=spoke_rot, rotation_mode='anchor',
                zorder=2, **font_kw_ring)

    # ---- half lollipop ---- #

    dashed_color = 'darkgrey'


    fig.canvas.draw()
    _r_sample = r_max * 0.5
    _p_lo = ax.transData.transform((0, _r_sample))
    _p_hi = ax.transData.transform((0, _r_sample + 1))
    _px_per_unit = np.hypot(*(_p_hi - _p_lo))
    _marker_radius_pts = semicircle_radius / 2.0 
    _marker_radius_px  = _marker_radius_pts * fig.dpi / 72.0
    marker_radius_data = _marker_radius_px / _px_per_unit * 0.80

    for j in range(N):
        if np.any(np.isnan(data[j])):
            continue  # skip spacer spokes

        v0 = data[j, 0]  # slope comparison   (model 0 / colors[0])
        v1 = data[j, 1]  # volume adjusted   (model 1 / colors[1])
        v_min = min(v0, v1)
        v_max = max(v0, v1)
        ang = theta[j]

        v_spoke_end = v_max + marker_radius_data

        if spoke_color == 'split':
            # Determine which model has the lower / higher value
            if v0 <= v1:
                color_lo, color_hi = colors[0], colors[1]
            else:
                color_lo, color_hi = colors[1], colors[0]

            # --- Solid spoke: lower-value color from origin to v_min ---
            ax.plot([ang, ang], [0, v_min],
                    color=color_lo, linewidth=1.5, linestyle='solid',
                    zorder=1)
            # --- Solid spoke: higher-value color from v_min to v_spoke_end ---
            ax.plot([ang, ang], [v_min, v_spoke_end],
                    color=color_hi, linewidth=1.5, linestyle='solid',
                    zorder=1)
        else:
            # --- Solid spoke: single color from origin to v_spoke_end ---
            ax.plot([ang, ang], [0, v_spoke_end],
                    color=spoke_color, linewidth=1.5, linestyle='solid',
                    zorder=1)

        # --- Dashed spoke from spoke end to outer rim ---
        ax.plot([ang, ang], [v_spoke_end, r_max],
                color=dashed_color, linewidth=0.75, linestyle='dashed',
                zorder=1)

        # --- Half-circle markers ---
        # Use MarkerStyle with fillstyle='left'/'right' and rotate
        # the marker so its split line aligns with the spoke.
        marker_rotation = -np.degrees(ang)

        for model_idx, value in enumerate([v0, v1]):
            color = colors[model_idx]
            # Model 0 -> clockwise side  -> fillstyle='right'
            # Model 1 -> counter-clockwise side -> fillstyle='left'
            fs = 'right' if model_idx == 0 else 'left'
            m = MarkerStyle('o', fillstyle=fs)
            m = m.transformed(Affine2D().rotate_deg(marker_rotation))

            ax.plot(ang, value, marker=m, markersize=semicircle_radius,
                    color=color, markerfacecoloralt='none',
                    markeredgecolor='none', markeredgewidth=0,
                    zorder=5)

        lbl = metrics[j]
        bin_range = ''
        if '[' in lbl and '-' in lbl:
            bin_range = lbl.split('[')[1].split(']')[0]

        if bin_range:
            display_angle = 90.0 - np.degrees(theta[j])
            rot = display_angle - 90.0          # perpendicular to spoke

            rot_norm = (rot + 180) % 360 - 180
            if rot_norm < -90 or rot_norm > 90:
                rot_norm += 180

            font_kw = {'fontproperties': font_prop} if font_prop else {}
            ax.text(theta[j], r_max * 1.02, bin_range,
                    ha='center', va='center',
                    fontsize=16, color='black',
                    rotation=rot_norm, rotation_mode='anchor',
                    **font_kw)


    font_kw_axis = {'fontproperties': font_prop} if font_prop else {}
    ax.text(label_angle, r_max * 0.55, f'{mode} positive %',
            ha='center', va=va_above,
            fontsize=12, color='black', fontstyle='italic',
            rotation=spoke_rot, rotation_mode='anchor',
            zorder=2, **font_kw_axis)


    centre_angle = None
    for j in range(N):
        if np.any(np.isnan(data[j])):
            continue
        lbl = metrics[j]
        if '[' in lbl and '-' in lbl:
            bs = lbl.split('[')[1].split('-')[0]
            if bs.strip() == '0.45':
                centre_angle = theta[j]
                break
    if centre_angle is not None:
        curved_text = 'noise levels'
        arc_r = r_max * 1.12
        char_span = 0.022
        total_span = (len(curved_text) - 1) * char_span
        start_angle = centre_angle - total_span / 2

        for ci, ch in enumerate(curved_text):
            a = start_angle + ci * char_span
            disp = 90.0 - np.degrees(a)
            rot_ch = disp - 90.0
            rot_ch = (rot_ch + 180) % 360 - 180
            if rot_ch < -90 or rot_ch > 90:
                rot_ch += 180
            ax.text(a, arc_r, ch,
                    ha='center', va='center',
                    fontsize=16, color='black', fontstyle='italic',
                    rotation=rot_ch, rotation_mode='anchor',
                    zorder=2, **font_kw_axis)

    ax.set_xticks(theta)
    ax.set_xticklabels([''] * N)
    ax.set_yticklabels([])

    return fig


#####################################################
#                                                   #
#                     MAIN                          #
#                                                   #
#####################################################


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Create a semi-lollipop chart of significance ratios '
                    'from a single CSV file produced by simulations.m. '
                    'Each spoke represents a base image at a specific '
                    'noise bin, with half-circles indicating the two '
                    'model values.')
    parser.add_argument('csv', type=str,
                        help='Path to the SignificanceRatios.csv file '
                             'produced by simulations.m.')
    parser.add_argument('--title', type=str, default='Significance Ratios',
                        help='Title for the chart '
                             '(default: "Significance Ratios").')
    parser.add_argument('--colors', type=str, nargs=2,
                        default=['#0072be', '#77ad30'],
                        help='Two colors for slope comparison and '
                             'volume-covariate half-circles '
                             '(default: #0072be #77ad30).')
    parser.add_argument('--output', type=str, default=None,
                        help='Path to save the figure. If not set, the '
                             'figure is shown interactively.')
    parser.add_argument('--format', type=str, default='png',
                        choices=['png', 'pdf', 'svg'],
                        help='Output figure format (default: png).')
    parser.add_argument('--limit', type=int, default=100,
                        help='Y axis limit to use (default: 100).')
    parser.add_argument('--font', type=str, default=None,
                        help='Path to a .ttf or .otf font file to use '
                             'for spoke labels (e.g. DINPro-Regular.otf).')
    parser.add_argument('--mode', type=str, default='true',
                        help='Define modus of simulation '
                             '(either true or false positive).')
    parser.add_argument('--semicircle-radius', type=float, default=10.0,
                        help='Marker size for the half-circle markers '
                             'in points (default: 10).')
    parser.add_argument('--spoke-color', type=str, default='split',
                        help='Color for the solid spoke segments. Use '
                             '"split" to color each segment by the '
                             'corresponding model color (default). Any '
                             'other value is treated as a matplotlib '
                             'color name (e.g. "darkgrey", "#888888").')
    args = parser.parse_args()

    if not os.path.isfile(args.csv):
        raise FileNotFoundError(f"CSV file not found: {args.csv}")

    # ---- load font ---- #
    font_prop = None
    if args.font is not None:
        if not os.path.isfile(args.font):
            raise FileNotFoundError(f"Font file not found: {args.font}")
        font_prop = FontProperties(fname=args.font)

    labels, slope, covar = read_csv(args.csv, add_spacers=True)

    # ---- base names ---- #
    seen = set()
    basenames_order = []
    for lbl in labels:
        if not lbl:
            continue
        base = lbl.split(' [')[0]
        if base not in seen:
            seen.add(base)
            basenames_order.append(base)

    # ---- build title ---- #
    title = ""

    # ---- stack as (n_spokes x 2): col 0 = slope comparison, col 1 = volume covariate ---- #
    data = np.column_stack([slope, covar])

    fig = plot_semi_lollipop_chart(
        data, labels, title,
        colors=args.colors, ylim=(0, args.limit),
        semicircle_radius=args.semicircle_radius,
        font_prop=font_prop,
        mode=args.mode,
        spoke_color=args.spoke_color)
    plt.tight_layout()

    # ---- save ---- #
    if args.output is not None:
        out_dir = os.path.dirname(args.output)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        fig.savefig(args.output, dpi=300, bbox_inches='tight')
        print(f"Saved: {args.output}")
        plt.close('all')
    else:
        plt.show()
