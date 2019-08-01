import os
import sys

variants = []
vmd = "\"/Applications/VMD 1.9.3.app/Contents/Resources/VMD.app/Contents/MacOS/VMD\""
for file in os.listdir(sys.argv[1]):
	if file != "analysis" and not file.startswith("."):
		variants.append(file)

for v in variants:
	os.system("python /Users/vikram/Documents/georgetown/summer_2019/snp2sim/snp2sim.py \
	--runDIR /Users/vikram/Documents/georgetown/summer_2019/run_workflow/snp2simresults \
	--protein PDL1 --varAA %s --varResID %s --genStructures --VMDpath /Applications/VMD_1.9.3.app/Contents/Resources/VMD.app/Contents/MacOS/VMD \
	--mode varMDsim \
	--newStruct /Users/vikram/Documents/georgetown/summer_2019/snp2sim/example/PDL1.Vtype.pdb" %(v[-1], v[:-1]))

