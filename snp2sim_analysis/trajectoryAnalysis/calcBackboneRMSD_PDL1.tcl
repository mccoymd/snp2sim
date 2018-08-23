set trajFiles [glob *.pdb]
set rmsdFile [open "trajectoryRMSD.txt" "w"]
puts $rmsdFile "Protein\tVariant\tSimID\tIndex\tRMSD_BindingRes\tRMSD_Backbone"
foreach indTraj $trajFiles {
    mol new $indTraj waitfor all
    set stats [split $indTraj .]
    set nf [molinfo top get numframes]
    set refStruct [atomselect top backbone frame 0]
    set refPDL1res [atomselect top "backbone and resid 19 20 54 56 66 68 115 116 117 121 122 123 124 125" frame 0]
    for {set i 0} {$i < $nf} {incr i} {
	set curStruct [atomselect top backbone frame $i]
	set curPDL1res [atomselect top "backbone and resid 19 20 54 56 66 68 115 116 117 121 122 123 124 125" frame $i]
	set M [measure fit $curStruct $refStruct]
	$curStruct move $M
	puts -nonewline $rmsdFile "[lindex $stats 0]\t[lindex $stats 1]\t[lindex $stats 2]\t$i\t"
	puts -nonewline $rmsdFile "[measure rmsd $curPDL1res $refPDL1res]\t"
	puts $rmsdFile [measure rmsd $curStruct $refStruct]
    }
}

close $rmsdFile
