"""
Calculate the dG from ddG using cinnabar
"""

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

    print(df)


def from_csv_ref():
    fe = FEMap.from_csv("ref.csv")
    fe.generate_absolute_values()
    assert fe.check_weakly_connected()
    df = fe.get_absolute_dataframe()
    print(df.to_string())


new_interface()
