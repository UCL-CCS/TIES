"""
A case where the O1-O3 in the main joint area
are pointing in opposite direction.

This creates a side effect when we align the
Z molecule to A molecule, leading
to the alchemical regions pointing in other direction.

Project:
Ideally, we'd rotate the alchemical region
with the understanding that it should remain
in the original position.

For this we should have a flag
which conveys that both molecules
were docked and we do not want to disturb
where the alchemical region is.

Furthermore, because this issues exists,
we should always check if we are rotating the alchemical
region and issue a warning that we are moving it out
of the pocket.

Temporary solution used was FEgrow to recapture the
original binding.
"""

from ties import Pair

pair = Pair("0.mol2", "9806.mol2")

hybrid = pair.superimpose(superimposition_starting_pairs="O1-O3")

assert len(hybrid) == 26
