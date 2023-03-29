import math
import itertools
from collections import defaultdict
# from typing import Concatenate
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.tri as tri
import os
from os import path
import re

comp_type = {'B': 'J', 'L': 'L', 'R': 'R', 'C': 'C', 'I': 'I'}
units = dict(u=1.e-6, p=1.e-12, n=1e-9, m=1e-3, k=1e3)

def convert_value(value):
    # print(value)
    topvalue = value.split()[0]
    print("topvalue",topvalue)
    last = topvalue[-1]
    if last in "0123456789":
        # print(last)
        return float(topvalue)
    try:
        return float(topvalue[:-1]) * units[last]
    except ValueError:
        if value[:6] == "pulse(":
            value = make_pulse(value)
        return value

def make_pulse(input_string):
    assert input_string[:6] == "pulse("
    values = input_string[6:-1].split(',')
    float_values = (convert_value(v) for v in values)
    vinit, vpulse, delay, rise_time, fall_time, width, period = float_values

    def pulse(t):
        t = t % period if period > 0 else t
        if t < delay:
            return vinit
        if t < delay + rise_time:
            dt = t - delay
            return vinit + (vpulse - vinit)*dt/rise_time
        if t < delay + rise_time + width:
            return vpulse
        if t < delay + rise_time + width + fall_time:
            dt = t - delay - rise_time - width
            return vpulse + (vinit - vpulse)*dt/fall_time
        return vinit
    return pulse

def IsShuntResistor(lst,el): # Unnecessary
    elementclass = el[0][0]
    for pair in lst:
        if elementclass == "R" and [el(1),el(2)] in lst:
            return True
        else:
            return False

def JJ_models(line,maindict):
    lst = line.split()
    print("MODEL",lst)
    # adding model to dictionary
    for i in range(len(lst)):
        modelname = lst[1]
        print(modelname)
        maindict[modelname] ={}
        parameters = ''.join(lst[2:])
        mod_para = re.search('\(([^)]+)', parameters).group(1).split(",")
        for para in mod_para:
            para_sep = para.split("=")
            maindict[modelname][para_sep[0]] = convert_value(para_sep[1])


def cir_to_networkx(filename):
    G = nx.MultiDiGraph()
    G.mutual_inductance = defaultdict(dict)
    G.models = defaultdict(dict)
    G.model_name = defaultdict(dict)
    G.mod_of_junc = defaultdict(dict)
    edges_by_name = {}
    nodes_by_pos = {0: (0,0)}  # start with ground at origin
    pos_so_far = set([])
    rightmostpos = 0
    # models = {}
    with open(filename) as f:
        lines = f.readlines()
        title = lines.pop(0)
        mutualinductor = []
        coupling = []
        junction_nodes = []
        for line in lines:
            print(line, end='')
            # print("Nodes so far:",list(G))
            if line[0] in "*":
                continue
            if line[:5] == ".tran":
                _, dt, tmax, uic = line.split()
                continue
            if line[:5] == ".plot" or line[:5] == ".save":
                continue
            if line[:6] == ".model":
                JJ_models(line,G.model_name)
                continue
            # if line[:5] == ".plot" or line[:5] == ".save":
                # print("Found")
                # continue
            fields = line.split()
            # print(fields)
            tp =fields[0][0]
            if tp == "K":
                inductors = fields[1:-1]
                mutualinductor.append(fields[1:-1])
                # print("Mutual Ind",mutualinductor)
                print(fields)
                value = convert_value(fields[-1])
                coupling.append(convert_value(fields[-1]))

                continue
                # for L1, L2 in itertools.combinations(inductors, 2):
                #     print(L1,L2) # The issue is that K could come before L
                #     assert L1 in edges_by_name
                #     assert L2 in edges_by_name
                #     L1nodes = edges_by_name[L1]
                #     L2nodes = edges_by_name[L2]
                #     L1self_inductance = G.edges[L1nodes]['L']
                #     L2self_inductance = G.edges[L2nodes]['L']
                #     scale = math.sqrt(L1self_inductance * L2self_inductance)
                #     G.mutual_inductance[L1][L2] = value * scale
                # continue
            name = fields[0]

            if fields[1] == "GROUND":
                node = 1
            elif fields[2] == "GROUND":
                nbr = 1
            else:
                node = int(fields[1])
                nbr = int(fields[2])

            if tp == "B":
                # print(fields)
                phase = fields[3]
                jjmodel = fields[4].lower()
                value = " ".join(fields[5:])
                junction_nodes.append([node,nbr])
                G.mod_of_junc[fields[0]] = jjmodel
                G.models[fields[0]] = G.model_name[jjmodel]
                if value[:5] == "area=":
                    value = convert_value(value[5:])
            elif tp == "R":
                value = convert_value(" ".join(fields[3:]))
            elif tp in comp_type:
                value = convert_value(" ".join(fields[3:]))
                print("value = ",value)
            else:
                print("Line not clear: ",line)
                print("Ignored")
                continue

#           add nodes with positions
            # node_in_dict = node in nodes_by_pos
            # nbr_in_dict = nbr in nodes_by_pos
            # if node == "0":
            #     if not nbr_in_dict:
            #         nodes_by_pos[nbr] = (rightmostpos, 1)
            #         pos_so_far.add((rightmostpos, 1))
            #         rightmostpos += 1
            #     #if nbr_in_dict: then both nodes already positioned
            # elif nbr == "0":
            #     if not node_in_dict:
            #         nodes_by_pos[node] = (rightmostpos, 1)
            #         rightmostpos += 1
            # elif node_in_dict:
            #     if not nbr_in_dict:
            #         x, y = nodes_by_pos[node]
            #         nodes_by_pos[nbr] = (x + len(G[node]), y+1)

            for n in (node, nbr):
                if n in nodes_by_pos:
                    continue

            # add edges with data
            edata = {"name": name, comp_type[tp]: value}
            ekey = G.add_edge(node, nbr, **edata)
            edges_by_name[name] = (node, nbr, ekey)

        numcoupled = len(mutualinductor)
        indx = 0
        for [L1,L2] in mutualinductor: # Now we go over all the inductor pairs that are coupled.
            # print(L1,L2)
            assert L1 in edges_by_name
            assert L2 in edges_by_name
            # print(L1)
            # print(L2)
            L1nodes = edges_by_name[L1]
            L2nodes = edges_by_name[L2]

            L1self_inductance = G.edges[L1nodes]['L']
            L2self_inductance = G.edges[L2nodes]['L']
            # print(L1self_inductance)
            # print(L2self_inductance)
            scale = math.sqrt(L1self_inductance * L2self_inductance)
            G.mutual_inductance[L1][L2] = coupling[indx] * scale
            # print("Value", value)
            print(G.mutual_inductance[L1][L2])
            indx += 1
            # print(value*scale)
            continue
    print("Nodes:")
    for n, pos in G.nodes.data('pos'):
        print(f"{n}: {pos}")
    print("Edges:")
    for e in G.edges.data():
        print((e[0], e[1], "%s"%(e[2])))
        print("Passed")
    print("Mutual Inductance")
    for node, nbrdict in G.mutual_inductance.items():
        # print(G.mutual_inductance.items())
        for nbr, M in nbrdict.items():
        #    print(nbrdict.items())
           print(node, nbr, "%e"%M)
        #    print(node, nbr, "%e"%(G.edges[edges_by_name[node]]['L']))
#    nx.draw_networkx(G, pos)
#    nx.draw_networkx_edge_labels(G, pos, edge_labels=nx.get_edge_attributes(G,'label'))
    #plt.show()
    # print(models)
    return G


if __name__ == "__main__":
    #cir_to_networkx("L19p_Delay40p.cir")
    # cir_to_networkx("test_circuit.cir")
    # G = cir_to_networkx("jjdram_0816.cir")
    # G = cir_to_networkx("RLcir.cir")
    #G = cir_to_networkx("singleneuron_mod.cir")
    G = cir_to_networkx("singleneuron_input.cir")
    #G = cir_to_networkx("multineurons_20n/multineurons_exp_5_20n.cir")
    # cir_to_networkx("jjdram_mod.cir")
    # print(G.models["jj02"])
    # print(G.models)
    # print(G.mutual_inductance)
    # print(G.models[G.mod_of_junc["B0"]])
