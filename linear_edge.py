import re


def sorting(l):
    """
    :type l: object
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


file_ccd = open("data/GM12878.CTCF.clusters_all_edges.csv", "r")
ccd_edge = []
linear_edge = []
for k, line2 in enumerate(file_ccd):
    if k == 0:
        continue
    line2 = line2.split("\n")[0]
    ccd_edge.append(line2.split(",")[0])
    ccd_edge.append(line2.split(",")[1])
    # print(line2)
    linear_edge.append(line2)

ccd_edge = list(dict.fromkeys(ccd_edge))
x = sorting(ccd_edge)
ccd_edge = x
print(len(ccd_edge))
ccd_edge = list(dict.fromkeys(ccd_edge))
print(len(ccd_edge))

for i in range(len(ccd_edge) - 1):
    # print(ccd_edge[i])
    # print(line.split("\t")[0] + "," + ccd_edge[i] + "," + ccd_edge[i+1] + ",1,Linear_edge,linear_edge")
    linear_edge.append(ccd_edge[i] + "," + ccd_edge[i + 1] + ",1")
file_ccd_wr = open("data/GM12878_RNAPII_PET2+_chr8_edges_liner_edge.csv", "w")
file_ccd_wr.write("anchor_id_A,anchor_id_B,link_score\n")
for line_write in linear_edge:
    print(line_write)
    file_ccd_wr.write(line_write + "\n")
file_ccd_wr.close()
