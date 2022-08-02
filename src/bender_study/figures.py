"""Plot manuscript figures."""

from pkg_resources import resource_filename
import matplotlib.pyplot as plt


def set_style(style_path=None):
    """Set default plot style."""
    if style_path is None:
        style_path = resource_filename('bender_study', 'style/figures.mplstyle')
    plt.style.use(style_path)
