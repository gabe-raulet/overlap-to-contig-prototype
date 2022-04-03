from igraph import *
import numpy as np

# @func mm_to_igraph - read a read overlap matrix market
#                      file into an igraph object
def mm_to_igraph(filename):
    attributes = ("rc", "begV", "endV", "begH", "endH", "lenV", "lenH")
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

def tr_igraph_tofile(G, filename):
    with open(filename, "w") as f:
        f.write("row\tcol\tdir\tsuffix\tprefix\n")
        for e in G.es:
            u,v = e.vertex_tuple
            f.write("{}\t{}\t{}\t{}\t{}\n".format(
                    u['idx'],
                    v['idx'],
                    e['dir'],
                    e['suffix'],
                    e['prefix']
                ))

def add_tr_info(G, maxOverhang=25):
    G.es['dir'] = G.es['suffix'] = G.es['prefix'] = -1
    for e in G.es:
        if e['begV'] > 0 and e['begH'] < maxOverhang:
            e['prune'] = False
            e['suffix'] = e['lenH'] - e['endH']
            e['prefix'] = e['begV']
            e['dir'] = 1 if not e['rc'] else 3
        elif e['begH'] > 0 and e['begV'] < maxOverhang:
            e['prune'] = False
            e['suffix'] = e['begH']
            e['prefix'] = e['lenV'] - e['endV']
            e['dir'] = 2 if not e['rc'] else 0
        else:
            e['prune'] = True

def main():
    G = mm_to_igraph("toy2.mtx")
    add_tr_info(G)
    tr_igraph_tofile(G, "toy2.tr.txt")

###############################################
#                                             #
# (!rc && (begV > 0) && (begH < maxOverhang)) #
#                                             #
#        seqV                                 #
# ------------------->                        #
#     **|||||||||||**                         #
#     ---------------------->                 #
#              seqH                           #
#                                             #
# suffix(V, H) = lenH - endH = prefix(H, V)   #
# prefix(V, H) =        begV = suffix(H, V)   #
#                                             #
#    dir(V, H) = (V >--> H) (1)               #
#    dir(H, V) = (H <--< V) (2)               #
#                                             #
###############################################
#                                             #
# (!rc && (begH > 0) && (begV < maxOverhang)) #
#                                             #
#               seqV                          #
#        ------------------->                 #
#        **|||||||||**                        #
#  ------------------->                       #
#        seqH                                 #
#                                             #
# suffix(V, H) =        begH = prefix(H, V)   #
# prefix(V, H) = lenV - endV = suffix(H, V)   #
#                                             #
#    dir(V, H) = (V <--< H) (2)               #
#    dir(H, V) = (H >--> V) (1)               #
#                                             #
###############################################
#                                             #
# (rc && (begV > 0) && (begH < maxOverhang))  #
#                                             #
#        seqV                                 #
# ------------------->                        #
#      **||||||||||**                         #
#     <----------------------                 #
#              seqH                           #
#                                             #
# suffix(V, H) = lenH - endH = prefix(H, V)   #
# prefix(V, H) =        begV = suffix(H, V)   #
#                                             #
#    dir(V, H) = (V >--< H) (3)               #
#    dir(H, V) = (H >--< V) (3)               #
#                                             #
###############################################
#                                             #
# (rc && (begV > 0) && (begH < maxOverhang))  #
#                                             #
#              seqV                           #
#       ------------------->                  #
#       **|||||||||**                         #
#  <-----------------                         #
#        seqH                                 #
#                                             #
# suffix(V, H) =        begH = prefix(H, V)   #
# prefix(V, H) = lenV - endV = suffix(H, V)   #
#                                             #
#    dir(V, H) = (V <--> H) (0)               #
#    dir(H, V) = (H <--> V) (0)               #
#                                             #
###############################################


# @func mm_to_adj - read a read overlap matrix market
#                   file into an adjacency matrix
def mm_to_adj(filename):
    with open(filename, "r") as f:
        next(f)
        n1,n2,m = (int(v) for v in next(f).rstrip().split())
        assert n1 == n2
        n = n1
        A = np.zeros((n,n),dtype=int)
        for line in f.readlines():
            u,v = (int(v) for v in line.rstrip().split()[:2])
            A[v-1,u-1] = A[u-1,v-1] = 1
    return A

#  @func adj_to_igraph - create igraph object from adjacency matrix
def adj_to_igraph(A):
    G = Graph.Adjacency((A>0).tolist())
    return G



