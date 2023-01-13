import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np
import numpy.linalg as LA
import scipy.integrate
import time

from numpy.linalg import inv

import findtree
import construct_matrices as cm
import read_cir
import os

def ColgateJJSim(circuit,simtime):
    # Constants
    Phi0 = 2.067833848e-15


    ## construct network
    my_circuit = read_cir.cir_to_networkx(circuit)

    print("Check non-MultiGraph:", my_circuit.edges.data())
    print(my_circuit.is_multigraph())
    print(my_circuit.edges())
    # position nodes
    pos = {"0": (0,0), "1": (0,1), "2": (1,1), "3": (2,1), "4": (3,1), "5": (2,2)}

    ######################################## Done with circuit specific input data
    # inputs: my_circuit, Inductances, gamma, coupling_lambda
    # get edge names
    enames = nx.get_edge_attributes(my_circuit, 'name')
    print("enames", enames)
    # to fix bug in nx that will be fixed in v2.3
    #enames = {k[:2]: v for k, v in enames.items()}
    nodes_for_name = {v:(k if len(k)==3 else k+(0,)) for k,v in enames.items()}
    print("nodes_for_edge name:", nodes_for_name)

    ## Find the tree
    print("my_circuit.nodes",my_circuit.nodes)
    print("my_circuit.edges",my_circuit.edges)
    my_tree=findtree.circuit_tree(my_circuit)
    print('my_tree.edges:',my_tree.edges)
    print('my_tree.nodes:',my_tree.nodes)

    # draw it
    draw_me = False
    print("pos",pos)
    if draw_me:
        ecolors = [('k' if e in my_tree.edges else 'b') for e in my_circuit.edges]
        ewidth = [(5 if e in my_tree.edges else 1) for e in my_circuit.edges]
        # Now draw the circuit
        nx.draw(my_tree, pos=pos, with_labels=True)
        nx.draw_networkx_edges(my_circuit, pos=pos, edge_color=ecolors, width=ewidth)
        elabels = nx.draw_networkx_edge_labels(my_circuit, pos=pos, edge_labels=enames)


    # Create matrices (Check for errors!)
    F_mats, edge_types, num_types, external = cm.make_matrices(my_circuit, my_tree)
    F_JL, F_JB, F_JZ, F_CL, F_CB, F_CZ, F_KL, F_KB, F_KZ = F_mats
    Jedges, Cedges, Kedges, Lchords, Bchords, Zchords, Redges, RedgesC = edge_types
    numJ, numC, numK, numL, numB, numZ, numR, numRC = num_types
    iB, diB = external    # TODO add phix and dphix -- currently assumed to be zero.
    iBvalues = np.array([f(0) for f in iB])

    print('numJ:', numJ, 'numC:', numC, 'numK:', numK)
    print('numL:', numL, 'numB:', numB, 'numZ:', numZ)
    print('numR:',numR, 'numRC:', numRC)
    print('Check that the count is correct')


    # Create Inductance, resistors, and eta matrices
    L = np.zeros([numL, numL])
    for i,Li in enumerate(Lchords):
        L[i,i] = my_circuit.edges[Li]['L']
    print('L:\n', L)

    L_K = np.zeros([numK, numK])
    for i,Li in enumerate(Kedges):
        L_K[i,i] = my_circuit.edges[Li]['L']
    # print('L_K:\n', L_K)

    L_LK = np.zeros([numL, numK])

    # Add any cross inductances here
    for u, unbrs in my_circuit.mutual_inductance.items():
        for v, Muv in unbrs.items():
            uL = True
            ue = nodes_for_name[u]
            print("ue:",ue)
            print("Muv:",Muv)
            print("Lchords:",Lchords)
            print("Kedges:",Kedges)
            try:
                upos = Lchords.index(ue)
                print("upos:",upos)
            except:
                uL = False
                upos = Kedges.index(ue)
                # print("upos:",upos)
            vL = True
            ve = nodes_for_name[v]
            try:
                vpos = Lchords.index(ve)
                print("vpos:",vpos)
            except:
                vL = False
                vpos = Kedges.index(ve)
                # print("vpos:",vpos)
            if uL and vL:
                L[upos,vpos] = L[vpos,upos] = Muv
            elif not (uL or vL):
                L_K[upos,vpos] = L_K[vpos, upos] = Muv
                # print('L_K:\n', L_K)
                # print("Found!")
            elif uL:
                L_LK[upos,vpos] = Muv
            else:
                L_LK[vpos,upos] = Muv
    print('L_LK:\n', L_LK)
    print('L_K:\n', L_K)

    # Resistors # This would be an issue.
    Rz = np.zeros([numZ,numZ])
    invRz = np.zeros([numZ,numZ])
    for i, Zi in enumerate(Zchords):
        R = my_circuit.edges[Zi]['R']
        Rz[i,i] = R
        invRz[i,i] = 1/R
    print("R_z:\n", Rz)

    Rshunt = np.zeros(numR)
    invRshunt = np.zeros(numR)
    for i, Ri in enumerate(Redges):
        Rs = my_circuit.edges[Ri]['R']
        Rshunt[i] = Rs
        invRshunt[i] = 1/Rs
    print("Rshunt: \n", Rshunt)
    print("invRshunt:\n", invRshunt)

    RshuntC = np.zeros(numRC)
    invRshuntC = np.zeros(numRC)
    for i, Ri in enumerate(RedgesC):
        Rs = my_circuit.edges[Ri]['R']
        RshuntC[i] = Rs
        invRshuntC[i] = 1/Rs
    print("RshuntC: \n", RshuntC)
    print("invRshuntC:\n", invRshuntC)

    Rzs = [Rz, invRz, Rshunt, invRshunt, RshuntC, invRshuntC]

    # Junctions
    eta = np.zeros([numJ,numJ])
    inveta = np.zeros([numJ,numJ])
    area = np.zeros(numJ)
    invarea = np.zeros(numJ)
    for i, Ji in enumerate(Jedges):
        J = my_circuit.edges[Ji]['J']
        eta[i,i] = J
        inveta[i,i] = 1/J
        area[i] = J
        invarea[i] = 1/J
    print("eta:\n",eta)

    # Capacitors
    etaC = np.zeros([numC,numC])
    invetaC = np.zeros([numC,numC])
    for i, Ci in enumerate(Cedges):
        C = my_circuit.edges[Ci]['C']
        etaC[i,i] = C
        invetaC[i,i] = 1/C
    print("etaC:\n",etaC)
    etas = [eta, inveta, etaC, invetaC, area, invarea]


    L_mats = [L, L_K, L_LK]

    # Now adding the model information
    models = my_circuit.models

    rhs = cm.manipulate_matrices(L_mats, F_mats, Rzs, etas, edge_types, num_types, enames, external,models)

    # Simulate
    timerange = (0,simtime)
    yinit = np.array([0]*(2*numJ+2*numC+numZ))

    start_time = time.time()
    S = scipy.integrate.solve_ivp(rhs, timerange, yinit, atol=2.23e-14, rtol=2.23e-14,method = "DOP853")

    print("--- %s seconds ---" % (time.time() - start_time))

    # Plot

    for i in range(len(S.y[:])):
        plt.figure(i)
        if i in range(0,numJ):
            plt.plot(S.t[:],S.y[i][:],'-',color='#4258a1',label="Python")
            plt.ylabel("Phase")
        elif i in range(numJ,2*numJ):
            plt.plot(S.t[:],S.y[i][:]*1000,'-',color='#4258a1',label="Python")
            plt.ylabel(r"Voltage ($mV$)")
            plt.title(r"Junction Voltage")
        elif i in range(2*numJ,2*numJ+numC):
            plt.plot(S.t[:],S.y[i][:],'-',color='#4258a1',label="Python")
            plt.ylabel(r"Capacitor Phase")
        elif i in range(2*numJ+numC,2*numJ+2*numC):
            plt.plot(S.t[:],S.y[i][:]*1000,'-',color='#4258a1',label="Python")
            plt.ylabel(r"Capacitor Voltage ($mV$)")
        elif i in range(2*numJ+2*numC,2*numJ+2*numC+numZ):
            plt.plot(S.t[:],S.y[i][:],'-',color='#4258a1',label="Python")
            plt.ylabel(r"Resistor Flux ($Wb$)")

        plt.xlabel(r"Time $(s)$")
    #    plt.savefig("FIG/fig_"+str(i)+".png",dpi=800)
        plt.show()
    #    plt.close()

if __name__ == "__main__":
    ColgateJJSim("singleneuron_mod.cir", 4e-9)
