import os
import sys

variants = []
writeDir = "/Users/vikram/Documents/georgetown/summer_2019/trajs"
readDir = sys.argv[1].replace(" ", "\\ ")
for file in os.listdir(sys.argv[1]):
	variants.append(file.split(".")[1])

variants = set(variants)

for file in os.listdir(writeDir):
	if file.endswith("dcd"):
		variants.remove(file.split(".")[1])

for v in set(variants):
	with open(writeDir + "/add_" + v, "w") as tcl:
		for file in os.listdir(sys.argv[1]):
			if v == file.split(".")[1]:
				tcl.write("mol addfile " + readDir + "/" + file + " waitfor all\n")
		tcl.write("animate write dcd %s/PDL.%s.dcd beg 0 end -1 sel [atomselect top all]\n" %(writeDir, v))
		tcl.write("quit")
	os.system("/Applications/VMD\\ 1.9.3.app/Contents/Resources/VMD.app/Contents/MacOS/VMD -e " + writeDir + "/add_" + v)

