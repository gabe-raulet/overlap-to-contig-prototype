from igraph import *

# @func mm_to_igraph - read a read overlap matrix market
#                      file into an igraph object
def mm_to_igraph(filename):
    attributes = ("rc", "b0", "e0", "b1", "e1", "l0", "l1")
    with open(filename, "r") as f:
        next(f)
        n1,n2,m = (int(v) for v in next(f).rstrip().split())
        assert n1 == n2
        n = n1
        G = Graph(n, directed=True)
        G.vs['idx'] = [i for i in range(1,n+1)]
        arrays = [[] for i in range(len(attributes)+1)]
        for line in f.readlines():
            vals = [int(v) for v in line.rstrip().split()]
            arrays[0].append((vals[0]-1, vals[1]-1))
            for i in range(len(attributes)):
                arrays[i+1].append(vals[i+2])
    G.add_edges(arrays[0])
    for i,attribute in enumerate(attributes):
        G.es[attribute] = arrays[i+1]
    return G


