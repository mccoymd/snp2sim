import os
import sys
import csv

runDIR = sys.argv[1]
protein = sys.argv[2]
legacy = sys.argv[3]
variants = []
results = "%s/variantSimulations/%s/results" % \
								   (runDIR, protein)
for file in os.listdir(results):
	if file != "analysis" and not file.startswith("."):
		variants.append(file)
for variant in variants:
	scaff = "%s/variantSimulations/%s/config" % \
									   (runDIR, protein)
	oldscaff = "%s/variantSimulations/%s/results/%s/old_scaffold/%s.%s.bindingRes.log" % \
									   (runDIR, protein, variant, protein, variant)
	for file in os.listdir(scaff):
		if file.endswith(".featureTable.csv") and file.startswith(protein + "." + variant + "."):
			table = scaff + "/" + file
			break
	if not table or not os.path.exists(oldscaff):
		print("no featureTable or log")
	else:
		if legacy in ["true", "True", "t"]:
			os.system("/Users/vikram/Documents/georgetown/summer_2019/snp2sim/snp2sim_analysis/scaffoldAnalysis/legacy_scaffold.sh %s" %oldscaff)

		clustLogfile = open(oldscaff, "r")
		clusterMembership = [x.split(",") for x in clustLogfile.readlines()]
		total_frames = sum(len(x) for x in clusterMembership)
		del clusterMembership[-1]
		clusterMembership = [x for x in clusterMembership if x and len(x) > .09 * total_frames]
		frame_clust = [["frame","cluster"]]
		for x in range(len(clusterMembership)):
			for frame in clusterMembership[x]:
				frame_clust.append([int(frame) + 1, x + 1])
		clust = "%s/variantSimulations/%s/results/%s/old_scaffold/cluster_by_frame.csv" % \
													(runDIR, protein, variant)
		with open(clust, "w") as file:
			csvWriter = csv.writer(file, delimiter=',')
			csvWriter.writerows(frame_clust)
		print(table)
		print(clust)
		os.system("Rscript ~/Documents/georgetown/summer_2019/snp2sim/snp2sim_analysis/scaffoldAnalysis/old_cluster_visualize.R \
					%s %s" %(table, clust))


