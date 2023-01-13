import networkx as nx
from itertools import chain

def extend_tree(G, starting_edges):
    """ Create a spanning tree of G starting with certain edges.

    The `starting_edges` list of edges are automatically included in the tree.
    The tree is formed first using a depth first search. Then each edge in
    `starting_edges` is checked. If not in the tree, then the cycle formed by
    adding it to the tree is searched for an edge with an inductor. That inductor
    edge is removed from the tree and the starting_edge is added.

    The circuit should be a connected graph -- though it could be only
    connected via a ground node.

    As this is designed for circuits, G is assumed to be MultiDiGraph and
    the orientation of edges in G is maintained in the tree.

    """
    # make a neighbors function that works both edge directions for MultiDiGraphs
    def nbrs(G, node):
        return chain(iter(G.edges(node, keys=True)), iter(G.in_edges(node, keys=True)))

    # create undirected tree and add starting edges
    T = nx.MultiDiGraph()
    print("G:::",G.edges)
    print("G:::",G.nodes)

    # do a DFS to add edges to make a spanning tree
    visited = set()
    for start in G:
        if start in visited:
            continue
        print("Starting node:", start)
        visited.add(start)
        stack = [(start, nbrs(G, start))]
        print("nbrs of start:",list(nbrs(G, start)))
        print("Edges of start:", start, G.edges(start),G.in_edges(start))
        print("orig Stack now",stack)
        while stack:
            parent, next_edges = stack[-1]
            T.add_node(parent)
            try:
                print("Going for parent edges",parent)
                u, v, ekey = next(next_edges)
                print("Got the u,v,ekey",u,v,ekey)
                # get direction of edge
                child = v if (u == parent) else u
                assert(u == parent or v == parent)
                print("Looking at edge:", u,v,ekey,"child is:",child)
                if child not in visited:
                    stack.append((child, nbrs(G, child)))
                    visited.add(child)
                    T.add_edge(u, v, ekey)
                    print("Adding edge:", u,v,ekey,"child is:",child)
            except StopIteration:
                print("Stack popped", parent)
                stack.pop()
            print("Stack now",stack)

    # now check for starting_edges
    print("starting edges:",starting_edges)
    print("T.edges:",T.edges)
    print("T.nodes:",T.nodes)
    for u,v,k in starting_edges:
        if T.has_edge(u,v,k):
            continue
        TTT = nx.graphviews.generic_graph_view(T, create_using=nx.MultiGraph)
        print("TTT for u,v,k:",(u,v,k),TTT.nodes,TTT.edges)
        path = nx.shortest_path(TTT, u, v)
        for n,nbr in nx.utils.pairwise(path):
            # get direction of edge and key in T
            (parent, child) = (n,nbr) if nbr in T[n] else (nbr,n)
            # get edge key from the key dict T[parent][child]
            ekey = next(iter(T[parent][child]))
            assert(T.has_edge(parent,child,ekey))

            if 'L' in G[parent][child][ekey]:
                # found an inductor to trade for our edge
                T.remove_edge(parent, child, ekey)
                T.add_edge(u, v, k)
                break
        if not T.has_edge(u, v, k):
            print("Error:  There is a cycle in with no inductor.")
#            print("path:", path, "   checking edge:",(u,v,k))
#            print("starting_edges:",starting_edges)
#            print("tree edges:",T.edges)
#            print("circuit edges:",G.edges(keys=True, data=True))
            raise ValueError("Cycle with no inductor identified")

    return T

@nx.utils.not_implemented_for('graph')
@nx.utils.not_implemented_for('undirected')
def circuit_tree(G):
    """ Return a tree that splits the circuit G.

    Split the circuit so that capacitors and junctions
    are in the tree while resistors, current sources and
    voltage sources are not in the tree.  Inductors may
    or may not be in the tree.

    Assumes G is a MultiDiGraph.
    """
    resistor_edges = [(u,v,k) for (u,v,k,dd) in G.edges(data=True, keys=True) if 'R' in dd]
    current_source_edges = [(u,v,k) for (u,v,k,dd) in G.edges(data=True, keys=True) if 'I' in dd]
    voltage_source_edges = [(u,v,k) for (u,v,k,dd) in G.edges(data=True, keys=True) if 'V' in dd]
    capacitor_edges = [(u,v,k) for (u,v,k,dd) in G.edges(data=True, keys=True) if 'C' in dd]
    junction_edges = [(u,v,k) for (u,v,k,dd) in G.edges(data=True, keys=True) if 'J' in dd]
    inductor_edges = [(u,v,k) for (u,v,k,dd) in G.edges(data=True, keys=True) if 'L' in dd]
    none_edges = [(u,v,k) for (u,v,k,dd) in G.edges(data=True, keys=True) if not ({'R', 'I', 'V', 'C', 'J', 'L'} & dd.keys())]

    edges_toss = resistor_edges + current_source_edges  # ??? voltage sources?
    edges_keep = capacitor_edges + junction_edges + none_edges

    G_bare = nx.restricted_view(G, nodes=[], edges=edges_toss)
    print("edges_keep",edges_keep)
    print("G_bare nodes",G_bare.nodes)
    print("G_bare edges",G_bare.edges)
    tree = extend_tree(G_bare, edges_keep)
    # G_bare.show()
    print("THE TREE IS", tree)

    return tree




#####   Test the script  (the "if" makes it so this won't run if you import this file)
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    plt.ion()     # turn interactive mode on if it isnt already

    # First "simple" example
    LG=nx.MultiDiGraph()
    pos = {1:(1,1), 2:(2,1), 3:(3,1), 4:(2,0)}
    LG.add_edge(1,2,R=10, name='R')
    LG.add_edge(2,3,L=10, name='L')
    LG.add_edge(3,4,J=10, name='J1')
    LG.add_edge(4,1,J=10, name='J2')

    tree = circuit_tree(LG)
    print("Tree edges: {}".format(tree.edges))
    print("Should be: [(2, 3, 0), (3, 4, 0), (4, 1, 0)]")

    plt.figure(1)
    plt.clf()
    nx.draw_networkx(LG, pos, with_labels=True)
    enames= nx.get_edge_attributes(LG, 'name')
    elabels = {(u, v): name for ((u, v, ekey), name) in enames.items()}
    nx.draw_networkx_edge_labels(LG, pos, elabels)
    nx.draw_networkx(tree, pos, width=3)

    # Another test graph
    #                          Ib
    #                          |
    #         Iin ->1--Lc-2-Jc-3
    #               |          Lp
    #               Ls         4--Lsyn--5--Rsyn-6
    #               |          Jp      Csyn     |
    #               7----------8--------9-------10
    #                         grd
    LG=nx.MultiDiGraph()
    pos = {1:(1,3), 2:(2,3), 3:(3,3), 4:(3,2), 5:(4.5,2), 6:(6,2),
           7:(1,1), 8:(3,1), 9:(4.5,1), 10:(6,1), 11:(0,3), 12:(3,4)}
    LG.add_edge(1,2,L=10, name="Lc")
    LG.add_edge(2,3,J=10, name="Jc")
    LG.add_edge(3,4,L=10, name="Lp")
    LG.add_edge(4,8,J=10, name="Jp")
    LG.add_edge(1,7,L=10, name="Ls")
    LG.add_edge(4,5,L=10, name="Lsyn")
    LG.add_edge(5,6,R=10, name="Rsyn")
    LG.add_edge(5,9,C=10, name="Csyn")
    LG.add_edge(6,10)
    LG.add_edge(7,8)
    LG.add_edge(8,9)
    LG.add_edge(9,10)

    LG.add_edge(11,1,I=1, name="Iin")
    LG.add_edge(12,3,I=1, name="Ib")
    LG.add_edge(7, 11)
    LG.add_edge(7, 12)

    tree = circuit_tree(LG)
    print("Tree edges: {}".format(tree.edges))
    true_edges = [(1, 2, 0), (2, 3, 0), (3, 4, 0), (4, 8, 0),
                  (8, 9, 0), (9, 10, 0), (6, 10, 0), (5, 9, 0),
                  (7, 8, 0), (7, 11, 0), (7, 12, 0)]
    assert(set(tree.edges) == set(true_edges))
    print("Tree edges correctly found")

    plt.figure(2)
    plt.clf()
    nx.draw_networkx(LG, pos, with_labels=True)
    enames= nx.get_edge_attributes(LG, 'name')
    elabels = {(u, v): name for ((u, v, ekey), name) in enames.items()}
    nx.draw_networkx_edge_labels(LG, pos, elabels)
    nx.draw_networkx(tree, pos, width=3)

    # test cycle of junctions
    print()
    print("Testing if inductor-less cycles are found:")
    G = nx.MultiDiGraph()
    G.add_edge(1,2,J=1)
    G.add_edge(1,3,J=1)
    G.add_edge(2,3,J=1)
    try:
        tree = circuit_tree(G)
    except ValueError:
        print("Cycle correctly found!")
    else:
        print("Cycle not found even though there was one!")
        print("Error!")
    print()
