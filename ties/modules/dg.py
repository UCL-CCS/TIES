"""
Calculate the dG from ddG using cinnabar
"""

import argparse
import warnings
from pathlib import Path

from openff.units import unit
from cinnabar.measurements import Measurement, ReferenceState
from cinnabar import FEMap

"""
Generate data. 

Idea 1: Globs: 
Take a random number of points. We know the centres
and can sample a distributions from each centre mu 
with different gaussian sigma. 

Then, feeding it directly into here, we should be able to recover 
the original mu centres. 
"""

kcal_mol = unit.kilocalorie_per_mole


def apply_correction(df, exp_value=None, exp_label=None):
    if exp_value is None:
        exp_rows = df[not df["computational"]]
        if len(exp_rows) == 0:
            if len(df[df["label"] == "0"]) == 1:
                warnings.warn(
                    "No experimental point found. Found label '0'. "
                    "Using the shift to the label '0' to have the value 0"
                )
                exp_value = 0
                exp_label = "0"
            else:
                warnings.warn("No experimental point found. Returning only MLE")
                return df
        else:
            assert len(exp_rows) == 1

            # grab the original value
            exp_value = exp_rows["DG (kcal/mol)"].values[0]
            exp_label = exp_rows["label"].values[0]

    # cinnabar applied a shift (to 0 apparently)
    # see #111

    # use the MLE value of your experimental point to calculate the shift
    exp_val_mle = df[(df["label"] == exp_label) & (df["source"] == "MLE")][
        "DG (kcal/mol)"
    ].values[0]

    # calculate the shift
    shift = exp_value - exp_val_mle

    # apply shift to all MLE points
    df.loc[df["source"] == "MLE", "DG (kcal/mol)"] += shift

    return df


def new_interface():
    exp_value = -6  # kcal_mol
    exp_label = "0"

    fe = FEMap()
    # add the experimental point
    fe.add_measurement(
        Measurement(
            labelA=ReferenceState(),
            labelB="0",
            DG=exp_value * kcal_mol,
            uncertainty=0.10 * kcal_mol,
            computational=False,
        )
    )
    fe.add_measurement(
        Measurement(
            labelA=exp_label,
            labelB="1",
            DG=-1 * kcal_mol,
            uncertainty=0.10 * kcal_mol,
            computational=True,
        )
    )
    fe.add_measurement(
        Measurement(
            labelA="1",
            labelB="2",
            DG=-1 * kcal_mol,
            uncertainty=0.10 * kcal_mol,
            computational=True,
        )
    )

    fe.generate_absolute_values()
    df = fe.get_absolute_dataframe()
    df = apply_correction(df, exp_value, exp_label)
    print(df)


def from_csv_ref(csv):
    fe = FEMap.from_csv(csv)
    fe.generate_absolute_values()
    df = fe.get_absolute_dataframe()

    df = apply_correction(df)
    print(df.to_string())
    return df


parser = argparse.ArgumentParser()
parser.add_argument(
    "-csv",
    metavar="str",
    dest="csv",
    type=Path,
    required=True,
    help="A CSV file with energies",
)
parser.add_argument(
    "-out",
    metavar="str",
    dest="corrected_csv_name",
    type=Path,
    required=False,
    default="dG.csv",
    help="A CSV file with energies",
)

if __name__ == "__main__":
    args = parser.parse_args()

    df = from_csv_ref(args.csv)
    df.to_csv(args.corrected_csv_name)
