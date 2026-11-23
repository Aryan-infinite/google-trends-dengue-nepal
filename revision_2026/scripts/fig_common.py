"""Shared settings for manuscript figures (Epidemiology & Infection specification).

Journal rules used here (Cambridge 'Preparing your materials', read 23 Nov 2026):
TIF at 100% of final size; 1000-1200 dpi for line drawings; CMYK; 9 pt Arial.
The journal does not publish a column width, so figures are 170 mm wide
(Cambridge double-column width); change WIDTH_MM if production asks for another size.
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from PIL import Image, ImageCms

DPI = 1200
WIDTH_MM = 170
MM = 1 / 25.4
FONT_PT = 9
# Colour profiles: set CMYK_ICC to a FOGRA39 (ISO Coated v2) profile if yours is elsewhere.
# The submitted figures used the FOGRA39L_coated.icc profile shipped with TeX Live 2023.
CMYK = os.environ.get('CMYK_ICC', '/usr/share/texlive/texmf-dist/tex/generic/colorprofiles/FOGRA39L_coated.icc')

INK = '#0b0b0b'        # observed cases
INK2 = '#52514e'       # axes, secondary text
GRID = '#d9d8d3'
BLUE = '#2a78d6'       # "dengue" search series
ORANGE = '#eb6834'     # control / second model


def pick_font():
    """Arial if installed, otherwise Liberation Sans (metrically identical to Arial)."""
    names = {f.name for f in font_manager.fontManager.ttflist}
    return 'Arial' if 'Arial' in names else 'Liberation Sans'


def setup():
    font = pick_font()
    plt.rcParams.update({
        'font.family': font, 'font.size': FONT_PT, 'axes.labelsize': FONT_PT,
        'xtick.labelsize': FONT_PT, 'ytick.labelsize': FONT_PT, 'legend.fontsize': FONT_PT,
        'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
        'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
        'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4,
        'axes.spines.top': False, 'axes.spines.right': False,
        'svg.fonttype': 'none', 'pdf.fonttype': 42,
    })
    return font


def save_tiff_cmyk(fig, path):
    """Render at DPI, convert sRGB -> CMYK (FOGRA39, perceptual), write LZW-compressed TIFF."""
    tmp = path.replace('.tif', '_rgb.png')
    fig.savefig(tmp, dpi=DPI, facecolor='white')
    rgb = Image.open(tmp).convert('RGB')
    if os.path.exists(CMYK):
        xf = ImageCms.buildTransform(ImageCms.createProfile('sRGB'), CMYK, 'RGB', 'CMYK',
                                     renderingIntent=ImageCms.Intent.PERCEPTUAL)
        cmyk = ImageCms.applyTransform(rgb, xf)
    else:
        print('WARNING: CMYK profile not found; using an uncalibrated conversion. Set CMYK_ICC.')
        cmyk = rgb.convert('CMYK')
    cmyk.save(path, compression='tiff_lzw', dpi=(DPI, DPI))
    preview = path.replace('.tif', '_preview.png')
    rgb.resize((rgb.width // 6, rgb.height // 6), Image.LANCZOS).save(preview)
    os.remove(tmp)
    return preview
