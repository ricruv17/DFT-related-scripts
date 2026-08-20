"""
Rotate and/or translate atomic structures using ASE.

Supports sequential rotations:
    --rotate x 10 y 20 z 30

or arbitrary axes:
    --rotate 0 0 1 45  1 1 0 30
"""

import numpy as np
from ase.io import read, write
import os
import argparse
import ast

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# ------------------ ARGPARSE SETUP ------------------
parser = argparse.ArgumentParser(
    description="Rotate and/or translate an atomic structure."
)

parser.add_argument(
    "filename",
    help="Input structure file"
)

parser.add_argument(
    "--rotate",
    nargs='+',
    help=(
        "Sequential rotations.\n"
        "Examples:\n"
        "  --rotate x 30 y 15 z 90\n"
        "  --rotate 0 0 1 45  1 1 0 30"
    )
)

parser.add_argument(
    "--translate",
    nargs='+',
    help='Translation vector: x y z OR "[x,y,z]"'
)

args = parser.parse_args()

prefix, suffix = os.path.splitext(args.filename)

# ------------------ LOAD STRUCTURE ------------------
atoms = read(args.filename)

output_filename = f"{prefix}_modified{suffix}"

# ------------------ ROTATIONS ------------------
if args.rotate:

    tokens = args.rotate
    i = 0

    while i < len(tokens):

        # ---------------- AXIS STRING ----------------
        if tokens[i] in ["x", "y", "z"]:

            if i + 1 >= len(tokens):
                raise ValueError(
                    f"Missing angle after axis '{tokens[i]}'"
                )

            axis = tokens[i]
            angle = float(tokens[i + 1])

            i += 2

        # ---------------- VECTOR AXIS ----------------
        else:

            if i + 3 >= len(tokens):
                raise ValueError(
                    "Vector rotations require: ax ay az angle"
                )

            axis = list(map(float, tokens[i:i + 3]))
            angle = float(tokens[i + 3])

            i += 4

        atoms.rotate(
            angle,
            axis,
            center=atoms.get_center_of_mass()
        )

        print(f"Rotated {angle}° around axis {axis}")

# ------------------ TRANSLATION ------------------
if args.translate:

    if len(args.translate) == 1:
        translation = ast.literal_eval(args.translate[0])

    elif len(args.translate) == 3:
        translation = list(map(float, args.translate))

    else:
        raise ValueError("Invalid --translate format")

    atoms.translate(translation)

    print(f"Translated by {translation}")

# ------------------ SAVE ------------------
write(output_filename, atoms)

print(f"Saved to {output_filename}")
