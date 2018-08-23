set trajFiles [glob *.pdb]
set rmsdFile [open "trajectoryRMSD.txt" "w"]
puts $rmsdFile "Protein\tVariant\tSimID\tIndex\tRMSD"
foreach indTraj $trajFiles {
    mol new $indTraj waitfor all
    set stats [split $indTraj .]
    set nf [molinfo top get numframes]
    set refStruct [atomselect top backbone frame 0]
    for {set i 0} {$i < $nf} {incr i} {
	set curStruct [atomselect top backbone frame $i]
	set M [measure fit $curStruct $refStruct]
	$curStruct move $M
	puts -nonewline $rmsdFile "[lindex $stats 0]\t[lindex $stats 1]\t[lindex $stats 2]\t$i\t"
	puts $rmsdFile [measure rmsd $curStruct $refStruct]
    }
}

close $rmsdFile
