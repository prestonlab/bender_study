"""Plot manuscript figures."""

import os
import re
import glob
from pkg_resources import resource_filename
import matplotlib.pyplot as plt


def set_style(style_path=None):
    """Set default plot style."""
    if style_path is None:
        style_path = resource_filename('bender_study', 'style/figures.mplstyle')
    plt.style.use(style_path)


def read_images(pool, im_dir):
    """Read image files for a set of items."""
    images = {}
    for i, item in pool.iterrows():
        # file names should have underscores instead of spaces
        file_base = item.stim.replace(' ', '_')
        im_base = os.path.join(im_dir, item.subcategory, file_base)

        # search for a matching image file for this item
        res = [f for f in glob.glob(im_base + '.*') if re.search(r'\w+\.(png|jpg)', f)]
        if not res:
            raise IOError(f'No file found matching: {im_base}')
        elif len(res) > 1:
            raise IOError(f'Multiple matches for: {im_base}')
        im_file = res[0]

        # read the image
        image = plt.imread(im_file)
        images[item.stim] = image
    return images
