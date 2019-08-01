import os
out = "searchSpace.txt"
scaffRes = "54 117 115 56 123 113 66 68 58 121 63 122 76 73"
boxTCL = "box.tcl"
templatePDB = "/Users/vikram/Documents/georgetown/summer_2019/run_workflow/snp2simresults/variantSimulations/PDL1/structures/PDL1.template.pdb"
clustTCL = open(boxTCL,"w+")
clustTCL.write("set output [open %s w]\n" % out)

clustTCL.write("mol new %s\n" % templatePDB)
clustTCL.write("set pocket [atomselect top \"resid "+scaffRes+"\"]\n")
clustTCL.write("set cen [measure center $pocket]\n")
clustTCL.write("set bound [measure minmax $pocket]\n")

clustTCL.write("puts -nonewline $output \"center_x = \"\n")
clustTCL.write("puts $output [lindex $cen 0]\n")
clustTCL.write("puts -nonewline $output \"center_y = \"\n")
clustTCL.write("puts $output [lindex $cen 1]\n")
clustTCL.write("puts -nonewline $output \"center_z = \"\n")
clustTCL.write("puts $output [lindex $cen 2]\n")

clustTCL.write("puts -nonewline $output \"size_x = \"\n")
clustTCL.write("puts $output [expr abs([lindex [lindex $bound 1] 0] - [lindex [lindex $bound 0] 0])]\n")
clustTCL.write("puts -nonewline $output \"size_y = \"\n")
clustTCL.write("puts $output [expr abs([lindex [lindex $bound 1] 1] - [lindex [lindex $bound 0] 1])]\n")
clustTCL.write("puts -nonewline $output \"size_z = \"\n")
clustTCL.write("puts $output [expr abs([lindex [lindex $bound 1] 2] - [lindex [lindex $bound 0] 2])]\n")

clustTCL.write("close $output\n")
clustTCL.write("quit\n")
clustTCL.close()

VMDpath = "/Applications/VMD\\ 1.9.3.app/Contents/Resources/VMD.app/Contents/MacOS/VMD"
vmdCommand = "%s -dispdev text -e box.tcl" % (VMDpath)
os.system(vmdCommand)  