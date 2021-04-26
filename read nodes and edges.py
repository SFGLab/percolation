import re


def sorting(l):
    """

    :type l: object
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


nodes = []
file_whole_genome = open("GM12878.CTCF.clusters_PET4+_nodes.csv", "r")
edge = []
for k, link in enumerate(file_whole_genome):
    if k == 0:
        continue
    link = link.split("\n")[0]
    nodes.append(link)

file_whole_genome.close()

nodes = list(dict.fromkeys(nodes))
x = sorting(nodes)
nodes = x
nodes = list(dict.fromkeys(nodes))

edges = []
file_whole_genome = open("GM12878.CTCF.clusters_PET4+_edges.csv", "r")
for k, link in enumerate(file_whole_genome):
    if k == 0:
        continue
    link = link.split("\n")[0]
    edges.append(link)

file_whole_genome.close()

print(len(nodes))
print(len(edges))

file_edge_prob = open("edge_probability.csv", "w")
edge_prob = []
node_prev = "chr1"
node_new = ""
print(node_prev)
for node in nodes:
    node_new = node.split(":")[0]
    if node_prev != node_new:
        print(node.split(":")[0])
    node_prev = node_new
    wt = 0
    prob = []
    for edge in edges:
        sub = edge.split(",")
        if sub[1] == node or sub[2] == node:
            wt = wt + int(sub[3])
            prob.append(edge)
            # print("\t" + edge)
    # print(wt)edge_prob
    for edge in prob:
        sub = edge.split(",")
        prob_wt = float(int(sub[3]) / wt)
        edge = sub[1] + "," + sub[2] + "," + str(prob_wt)
        edge_prob.append(edge)
        file_edge_prob.write(edge + "\n")
        # print("\t", edge)

file_edge_prob.close()
