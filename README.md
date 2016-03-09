# qwrap
A fast PBC-wrapping function for VMD

To install under Linux and the like:
```
$ tar xf qwrap.tar.gz (in own directory) 
$ make
```

To use in VMD:
```
> load path/to/qwrap.so 
> qwrap [first n] [last n] [compound none|res|beta [refatoms occ|none]] [center <seltext>]
```
* `first`, `last`: frames to be wrapped
* `compound`: wrap groups of atoms instead of individual atoms. Groups can be VMD's `res`idues, or custom groups defined by a common value of the `beta` parameter. 
* `refatoms`: if `compound` is set, each group may be wrapped according to its geometric center, or to the center of a set of reference atoms within the group. Then reference atoms are defined by nonzero `occ`upancy. An example use is to wrap lipid molecules as residues, and give all phosphorus atoms nonzero occupancy so that each lipid has its phosphorus atom in the center unit cell.
* `center`: the geometric center of the given selection text will be translated to (0, 0, 0). 

Compared with PBCTools it's a coarse and rigid tool, mostly built to answer my own needs (e.g. orthorhombic cells only, the center of the cell is (0,0,0)...). But in my hands, it is 10 to 30 times faster. 

The compound-based wrapping doesn't use a specific reference atom, but the center of the wrapping group. Which by the way can be defined either as residues, or by setting custom flags in the beta field, so you can keep a macromolecule in one piece. The flags should be integers, and could be alternating 0s and 1s or anything else, they only need to change when the next atom belongs to a different wrapping group. 

A more mature version of this could eventually provide a fast back-end for the pure-Tcl pbc wrap. 
It's work in progress, with the usual caveats as to bugs etc. Please let me know if it is useful, and if you improve it, share it back! 
