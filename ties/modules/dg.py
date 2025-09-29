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


def basic():
    kcal_mol = unit.kilocalorie_per_mole
    g = ReferenceState(label="0")
    experimental_result1 = Measurement(
        labelA=g,
        labelB="CAT-13a",
        DG=-8.83 * kcal_mol,
        uncertainty=0.10 * kcal_mol,
        computational=False,
    )
    experimental_result2 = Measurement(
        labelA=g,
        labelB="CAT-17g",
        DG=-9.73 * kcal_mol,
        uncertainty=0.10 * kcal_mol,
        computational=False,
    )
    # Load/create calculated results
    calculated_result = Measurement(
        labelA="CAT-13a",
        labelB="CAT-17g",
        DG=0.36 * kcal_mol,
        uncertainty=0.11 * kcal_mol,
        computational=True,
    )
    # Incrementally created FEMap
    fe = FEMap()
    fe.add_measurement(experimental_result1)
    fe.add_measurement(experimental_result2)
    fe.add_measurement(calculated_result)

    df = fe.get_absolute_dataframe()
    print(df)


def better():
    g = ReferenceState()
    experimental_result1 = Measurement(
        labelA=g,
        labelB="0",
        DG=-6 * kcal_mol,
        uncertainty=0.10 * kcal_mol,
        computational=False,
    )

    # Incrementally created FEMap
    fe = FEMap()
    fe.add_measurement(experimental_result1)
    fe.add_relative_calculation(
        labelA="0", labelB="1", value=-1 * kcal_mol, uncertainty=0.10 * kcal_mol
    )

    fe.generate_absolute_values()
    df = fe.get_absolute_dataframe()
    print(df)


def fine_but_old_interface():
    exp_value = -6
    exp_label = "0"

    fe = FEMap()
    fe.add_experimental_measurement(
        exp_label, value=exp_value * kcal_mol, uncertainty=0.10 * kcal_mol
    )
    fe.add_relative_calculation(
        labelA=exp_label, labelB="1", value=-1 * kcal_mol, uncertainty=0.10 * kcal_mol
    )
    fe.add_relative_calculation(
        labelA="1", labelB="2", value=-1 * kcal_mol, uncertainty=0.10 * kcal_mol
    )

    fe.generate_absolute_values()
    df = fe.get_absolute_dataframe()

    # look up the experimental value to find the shift
    exp_val_mle = df[(df["label"] == exp_label) & (df["source"] == "MLE")][
        "DG (kcal/mol)"
    ].values[0]

    # calculate the needed shift
    shift = exp_value - exp_val_mle
    df["DG (kcal/mol)"] += shift

    print(df)


def new_interface():
    exp_value = -6
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
print("hi")
