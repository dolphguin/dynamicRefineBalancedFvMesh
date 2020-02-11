# dynamicRefineBalancedFvMesh

# Description
Allows dynamic load balancing when working with dynamic meshes (not very stable
yet...)

A multi-field version has been merged, which extends the standard dynamic mesh
refinement capabilities to allow combining multiple refinement conditions based
on different fields.

# Requirements
Ubuntu-14.04LTS or later
(not tested for other systems)

OpenFOAM-6
(the scripts are prepared to work sourcing the default ubuntu-repository
version; othewise the line "WM_PROJECT_DIR=/opt/openfoam7" must be accordingly
edited)

# Installation
Save folder "dynamicRefineBalancedFvMesh" (the one with the "Allwmake/Allwclean"
scripts inside) into $HOME/OpenFOAM/\<username\>-7/applications/utilities/ (suggestion) and use
"Allwmake" script as usual e.g. "Allwmake -j".

PS: the flag "j" enables multi-thread compilation.

# Library linking
The utility can be linked in runtime by addding the corresponding entry in the
file "system/controlDict":

[...]
libs
(
    "libdynamicRefineBalancedFvMesh-of7.so"
    "libdynamicMultiFieldRefineFvMesh-of7.so"
    "libdynamicMultiFieldRefineBalancedFvMesh-of7.so"
);
[...]

# Use
The decomposition method for the domain redistribution (load balancing) is read
from the dictionary "system/balanceParDict". Ensure this file is present. It can
be initially copied from the template in the library.
Similarly, the dynamic mesh refinement method as well as the load balancing
switch are read from the "constant/dynamicMeshDict" file. It is suggested to use
the template provided within the library.

# Notes
When selecting the multi-field version "dynamicMultiFieldRefineBalancedFvMesh"
the fields shall be given in descending order of "maxRefinement" in order to get
an a proper refinement/unrefinement behaviour.

# TO-DO list
Additional refinement controls based on field gradients, field curls and mesh
regions that were present in the original implementation from Tyler V. (see
https://github.com/tgvoskuilen/meshBalancing) were kept by the moment, because
they may become useful in some cases. But not enough attention has been paid to
make them compatible witht the multi-field option...
EDIT: actually, only the regions have been kept, though they are by the moment
incompatible with the balancing option.

The current implementation lacks support for lagrangian particles.

Crashes have been reported when using high amounts of processes as well as high
number of cells (several millions). Some of them seem to be directly related to
Scotch library (methods "scotch" and "ptscotch", respectively), so we must
generate a simple test-case to report at their bug-reporting site. In the
meanwhile, using different decomposition methods e.g. "hierarchical",
"multiLevel", etc. can overcome those problems. Nevertheless, there seem to be
other crashes that are non related to Scotch library too. We shall get a simple
case to reproduce this latter kind of crashes and debug them.

Check for potential GNU-GPLv3 issues (are file headers OK?)


# Authors
Rodrigo Gómez Vázquez
Current development.
(https://tiss.tuwien.ac.at/person/64973.html)

Francisco Javier Codina Álvarez
Initial attempt to harmonize region-based refinement with the previously implemented
multi-field refinement.
(-)

Robert Feichtenschlager
Multi-field functionality.
(https://tiss.tuwien.ac.at/person/61423.html)

Tyler V.
Original load balancing implementation.
(https://github.com/tgvoskuilen)

OpenFOAM creators/developers
Dynamic mesh refinement, whole software infrastructure.
(https://openfoam.org)
