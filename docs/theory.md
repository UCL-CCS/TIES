# TIES Protocol

## Superimposition and defining the alchemical region

Any two pairs are superimposed using a recursive joint
traversal of two molecules starting from any two pairs.

A heuristics (on by default) reduces the search space
by selecting the rarer atoms that are present across the two
molecules as the starting points for the traversal,
decreasing substantially the computational cost.


## Charge treatment

TIES 20 supports the transformation between ligands
that have **the same net charge**.

We employ a dual topology approach which divides the atoms
in each transformation into three groups:

1. **Joint region**. This is the region of the molecule where
   the atoms are the same
   meaning that they are shared across the two ligands in the
   transformation.
2. **Disappearing region**. Atoms present only in the starting ligand
   of the transformation which are fully represented
   at lambda=0 and which will be scaled accordingly during the
   lambda progression.
3. **Appearing region**. Atoms present only in the ending ligand
   of the transformation and therefore not present at
   lambda=0. These atoms start appearing during the
   lambda progression and are fully represented at
   lambda=1.

When the two ligands in a transformation
are superimposed together, the treatment of charges
depends on which group they belong to.


## Joint region: matched atoms and their charges


In the joint region of the transformation,
first **--q-pair-tolerance** is used to determine
whether the two original atoms are truly the same atoms.
If their charges differ by more than this value (default 0.1e),
then the two atoms will be added to the alchemical regions
(Disappearing and appearing).

It is possible that a lot of matched atoms
in the joint region, with each pair being within 0.1e of each other,
cumulatively have rather different charges between
the starting and the ending ligand. For this reason, TIES 20 sums the
differences between the starting and the ending atoms in the joint region,
and if the total is larger than **-netqtol** (default 0.1e)
then we further expand the alchemical region until
the "appearing" and "disappearing" regions in the joint region
are of a sufficiently similar net charge.

Abiding by **-netqtol** rule has the further effect that,
inversely, the alchemical regions (disappearing and appearing regions),
will have very *similar* net charges - which is a necessary
condition for the calculation of the partial derivative of the potential energy
with respect to the lambda.

If **-netqtol** rule is violated, different schemes
for the removal of the matched atoms in the joint region
are tried to satisfy the net charge limit. The
scheme that removes fewest matched pairs,
is used. In other words, TIES 20 is trying to
use the smallest alchemical region possible while
satisfying the rule.

Note that we are not summing together the
absolute differences in charges in the joint region.
This means that if one atom pair has 0.02e charge difference,
and another pair has -0.02e charge difference, then their total is zero.
In other words, we are **not worried about the distribution
of the differences in charges** in the joint region.

The hydrogen charges are considered by absorbing them
into the heavy atoms.

The charges in the joint region for each pair are averaged.

The last step is `redistribution`, where the final goal
is that the net charge is the same in the Appearing and
in the Disappearing alchemical region. After
averaging the charges in the joint region, its overall
charge summed with the charge of each alchemical region
should be equal to the whole molecule net charge:
:math:`q_{joint} + q_{appearing} == q_{joint} + q_{disappearing} == q_{molecule}`.
Therefore, after averaging the charges, :math:`q_{molecule} - q_{joint} - q_{appearing}`
is distributed equally in the region :math:`q_{appearing}`.
The same rule is applied in :math:`q_{disappearing}`.

