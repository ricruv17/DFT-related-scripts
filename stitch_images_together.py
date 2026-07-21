"""
stitch_images_together.py

Program to stitch multiple images into a grid while preserving the exact
pixel resolution of each individual image.

Requirements:
    Python 3
    matplotlib
    numpy
    pillow (PIL)

This will create a 2×3 grid (as defined in the script) with layout:

    panel1  panel2  panel3
    panel4  panel5  panel6

and save it as:
    stitched.png

Notes:
    - The grid size is currently hardcoded as 2 rows × 3 columns:
          stitch_images_together(2, 3, files)
      You can change this at the bottom of the script.

    - All images are placed without resizing.

    - Output resolution will be exactly:
          width  = columns × max(image widths)
          height = rows    × max(image heights)

Usage:
    python stitch_images_together.py image1.png ... imageN.png output.png

Author: Ricardo Ruvalcaba Briones
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import PIL
import sys

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

files = sys.argv[1:-1]
out_file = sys.argv[-1]


def stitch_images_together(rows, columns, filenames, dpi=100):

    # Load all images first
    images = [PIL.Image.open(f) for f in filenames]

    # Assume all images same size (or use max if different)
    widths  = [img.size[0] for img in images]
    heights = [img.size[1] for img in images]

    max_w = max(widths)
    max_h = max(heights)

    # Total output resolution in pixels
    total_width  = columns * max_w
    total_height = rows    * max_h

    # Convert pixels → inches for matplotlib
    fig_w = total_width / dpi
    fig_h = total_height / dpi

    fig, axes = plt.subplots(rows, columns,
                             figsize=(fig_w, fig_h),
                             dpi=dpi)

    axes = np.atleast_1d(axes).flatten()

    for ax, img in zip(axes, images):

        ax.imshow(img,
                  interpolation='nearest',  # no smoothing
                  resample=False)           # preserve pixels

        ax.set_xlim(0, img.size[0])
        ax.set_ylim(img.size[1], 0)

        ax.axis('off')

    # Remove ALL padding
    plt.subplots_adjust(
        left=0,
        right=1,
        bottom=0,
        top=1,
        wspace=0,
        hspace=0
    )

    fig.savefig(out_file,
                dpi=dpi,
                bbox_inches=None,
                pad_inches=0)

    plt.close(fig)


stitch_images_together(3, 5, files)