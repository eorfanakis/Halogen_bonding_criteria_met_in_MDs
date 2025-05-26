# Halogen_bonding_criteria_met_in_MDs
Quantificate cases that geometrical criteria are met for XBs in MDs performed with GROMACS

A crucial type for various biochemical interactions, such in protein-ligand binding and in membrane permeability or solvation properties that influences ADME an Pharmacokinetic properties is gaining more and more importance. Conventional pairwise-additive Force Fields being used in most Molecular Dynamics simulations, underestimate or totally ignore the contribution of these less common interactions. Parametrization of these Force Fields or the us of polarizable Force Fields that treat accurately these, possibly, strong bonds. Furthermore, Virtual Sites or pseudoatoms are being used to examine with a more proper way the sigma hole.
With these being said, its becoming essential to spot cases where geometrical criteria are being met, if pairwise additive Force Fields are being used.

A case of an Alprazolam molecule, which is biased moving through a lipid membrane (DOPC molecules) at Umbrella Sampling simulations with GROMACS is being examined.

Lewis bases of phosphate groups and lipid tails oxygen atoms were selected, 8 atoms at total.
Firstly, the halogen atom (CL) of the alprazolam molecule (Chlorine) is being accounted as number 20 at .gro files, of the ALP molecule (Alprazolam).
Also oxygen atoms of DOPC lipid molecules are the atoms: O11, O12, O13, O14, O21, O22, O31 and O32 from:
$ gmx make_ndx -f umbrella_window.gro and make the files XB_O11.ndx, XB_012.ndx, XB_O13.ndx, XB_O14.ndx, XB_O21.ndx, XB_O22.ndx, XB_O31.ndx and XB_O32.ndx uploaded. Note the [ distance_pairs ] section and the [ angles ] section composed by the Chlorine atom, the Carbon atom with its being connected (19 in .gro files) and the Oxygen atom.

