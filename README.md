# qwrap — Fast PBC wrapping and unwrapping for VMD

## Installation
To install under Linux and the like:
```
$ tar xf qwrap.tar.gz (in own directory) 
$ make
# update the PLUGINDIR variable in Makefile
$ make install
# use 'sudo make install' if VMD is not installed within your user directory
```

## Usage
To use in VMD:
```
> package require qwrap
> qwrap [first n] [last n] [compound none|res|beta|fragment [refatoms occ|none]] [center <seltext>] [sel <seltext>]
> qunwrap [first n] [last n] [compound none|res|beta|fragment [refatoms occ|none]] [sel <seltext>]
```
* `sel`: selection text indicating atoms to be wrapped
* `first`, `last`: frames to be wrapped (defaults: 0 to -1, which means the last frame).
* `compound`: wrap groups of atoms instead of individual atoms. Groups can be VMD's `res`idues (default), or custom groups defined by a common *integer* value of the `beta` parameter. If VMD has connectivity information (i.e. bonds are correctly set), the option `fragment` can be useful.
* `refatoms`: if `compound` is set, each group may be wrapped according to its geometric center, or to the center of a set of reference atoms within the group. Then reference atoms are defined by nonzero `occ`upancy (when converted to integer, that is, *greater than 1*). An example use is to wrap lipid molecules as residues, and give all phosphorus atoms nonzero occupancy so that each lipid has its phosphorus atom in the center unit cell.
* `center`: the geometric center of the given selection text will be translated to (0, 0, 0). Not supported by qunwrap.

## Examples
A classic use to wrap the solvent around a protein would be:
```
qwrap sel "not protein" center protein
```

To repair a protein complex that has been split across PBC, load a trajectory where **the first frame has the oligomer in one piece**, then run:
```
qunwrap sel protein compound fragment
```
You can then run the command above to wrap the solvent around the complex.

## Advanced definition of wrapping groups

Wrapping groups can be defined either as residues, by fragment, or by setting custom flags in the beta field.
The flags should be integers, and could be alternating 0s and 1s or anything else, they only need to change when the next atom belongs to a different wrapping group.
Compound-based wrapping can use a custom set of reference atoms; the reference position of each wrapping group (compound) is the center of the reference atoms contained within this block (it may be a single atom).
If refatoms is not defined, the reference position is the center of the whole group.

## Limitations

Compared with PBCTools it's a less flexible tool, mostly built to answer my own needs (e.g. orthorhombic cells only, the center of the cell is (0,0,0)...). But in my hands, it is 10 to 30 times faster. 

Please let me know if it is useful, and if you improve it, share it back! 
