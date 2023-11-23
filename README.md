# qwrap — Fast PBC wrapping and unwrapping for VMD
[![DOI](https://zenodo.org/badge/31314121.svg)](https://zenodo.org/badge/latestdoi/31314121)

## Development news

### Version 1.6, 2023-04-30

Remove dynamically linked Tcl library and link tclstub library statically instead. This should make the shared object more portable.

### Version 1.5, 2022-01-03

This version reverts the change of unwrapping algorithm of the previous version, which was found to lead to distortion of molecules when unwrapping Gromacs trajectories atom-wise (`compound none`). Shifts are again computed as integer factors of the current cell basis vectors.
It turns out that the von Bülow unwrapping algorithm is not well-suited to qwrap's special treatment of groups of atoms (compounds), and I expect the bias in long-time diffusion from the qunwrap algorithm to be negligible.

### Version 1.4, 2021-05-25

As of version 1.4, unwrapping is done according to the algorithm of von Bülow, Bullerjahn and Hummer (https://arxiv.org/pdf/2003.09205.pdf) to obtain correct diffusion coefficient at long times.
Note that the algorithm used up to version 1.3 was already similar to this one, and gave only small, bounded deviations from the true long-time diffusion.


## Getting the library

### Downloading the pre-built binary for Linux
On Linux, you can just download the file [Linux/qwrap.so](https://github.com/jhenin/qwrap/raw/master/Linux/qwrap.so), save it somewhere and remember its full path.

### Building on Unix-like systems
First make sure Tcl development packages are installed (tcl-dev or tcl-devel).
```
tar xf qwrap.tar.gz (in own directory) 
make
```

For Macs with ARM64 processors (M1, M2, M3 etc.), try adapting the paths in `Makefile.OSX_ARM64` and running:
```
make -f Makefile.OSX_ARM64
```

Optionally, install as a VMD plugin:
```
# update the PLUGINDIR variable in Makefile
make install
# use 'sudo make install' if VMD is installed system-wide rather than in your user directory
```

## Loading the library into VMD

In the VMD terminal or TkConsole, type:
```
load <path>/qwrap.so
```
You can use the absolute or relative path, including '.' for the local directory.
Or, if you have installed qwrap as a VMD plugin:
```
package require qwrap
```
You can now use the `qwrap` and `qunwrap` commands.


## Usage
To use in VMD:
```
qwrap [first n] [last n] [compound none|res|beta|fragment [refatoms occ|none]] [center <seltext>] [sel <seltext>]
qunwrap [first n] [last n] [compound none|res|beta|fragment [refatoms occ|none]] [sel <seltext>]
```
* `sel`: selection text indicating atoms to be wrapped
* `first`, `last`: frames to be wrapped (defaults: 0 to -1, which means the last frame).
* `compound`: wrap groups of atoms instead of individual atoms. Groups can be VMD's `res`idues or bonded `frag`ments (default), or custom groups defined by a common *integer* value of the `beta` parameter. The default option (fragment) corresponds to the wrapping behavior of NAMD; however it only works correctly if VMD has accurate connectivity information (eg. read from a PSF file).
To undo GROMACS-style, atom-wise wrapping, use `compound none`.
* `refatoms`: if `compound` is set, each group may be wrapped according to its geometric center, or to the center of a set of reference atoms within the group. Then reference atoms are defined by nonzero `occ`upancy (when converted to integer, that is, *greater than 1*). An example use is to wrap lipid molecules as residues, and give all phosphorus atoms nonzero occupancy so that each lipid has its phosphorus atom in the center unit cell.
* `center`: the geometric center of the given selection text will be translated to (0, 0, 0). Not supported by qunwrap.

## Examples
A classic use to wrap the solvent around a protein would be:
```
qwrap sel "not protein" center protein
```

To repair a protein complex that has been split across PBC by NAMD (wrapAll), load a trajectory where **the first frame has the oligomer in one piece**, then run:
```
qunwrap sel protein
```
You can then run the command above to wrap the solvent around the complex.

To unwrap a GROMACS trajectory, make sure that the first frame has all molecules in one piece and run:
```
qunwrap compound none
```

## Advanced definition of wrapping groups

Wrapping groups can be defined either as residues, by fragment, or by setting custom flags in the beta field.
The flags should be integers, and could be alternating 0s and 1s or anything else, they only need to change when the next atom belongs to a different wrapping group.
Compound-based wrapping can use a custom set of reference atoms; the reference position of each wrapping group (compound) is the center of the reference atoms contained within this block (it may be a single atom).
If refatoms is not defined, the reference position is the center of the whole group.

## Limitations

Compared with PBCTools, qwrap is a less flexible tool (orthorhombic cells only, center of the cell fixed at (0,0,0)...), but it is many times faster.

Please let me know if it is useful, and if you improve it, share it back! 
