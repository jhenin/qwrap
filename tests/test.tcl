# Test relying on VMD API

load qwrap.so

# Check nearest-image wrapping in a triclinic cell.

set triclinic_mol [mol new 2_equ.pdb]
molinfo $triclinic_mol set {a b c alpha beta gamma} {1 1 1 60 60 60}
set triclinic_sel [atomselect $triclinic_mol "index 0"]
$triclinic_sel set {x y z} {{1.0 0.577350269 0.40824829}}
qwrap nearest yes compound none sel "index 0"
set wrapped [lindex [$triclinic_sel get {x y z}] 0]
set distance2 [veclength2 $wrapped]
if {abs($distance2 - 0.5) > 1e-4} {
	error "triclinic nearest-image test failed: squared distance is $distance2"
} else {
    puts "triclinic nearest-image test passed"
}
$triclinic_sel delete
mol delete $triclinic_mol

#exit
