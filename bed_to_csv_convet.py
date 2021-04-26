file_bed = open("data/GM12878.CTCF.singletons_cluster_PET_2_3.txt", "r")
file_csv = open("data/GM12878.CTCF.singletons_cluster_PET_2_3.csv", "w")

for line in file_bed:
    line = line.split("\n")[0]
    node_a = line.split("\t")[0] + ":" + line.split("\t")[1] + "-" + line.split("\t")[2]
    node_b = line.split("\t")[3] + ":" + line.split("\t")[4] + "-" + line.split("\t")[5]
    link_wt = line.split("\t")[6]
    wt_line = node_a + "," + node_b + "," + link_wt + "\n"
    file_csv.write(wt_line)

file_bed.close()
file_csv.close()


