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

# Constants
Phi0 = 2.067833848e-15


#my_circuit = read_cir.cir_to_networkx("singleneuron_fast.cir")
# my_circuit = read_cir.cir_to_networkx("jjdram_mod.cir")
my_circuit = read_cir.cir_to_networkx("jjdram_00.cir")

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
print("iB", iB)

# rise = 1e-14
# iB0 = []
# for current in iBvalues:
#     value = read_cir.convert_value("pulse(0,"+str(current)+",0,"+str(rise)+",0,0,0)")
#     # print("VALUE:\n",value)
# for e in Bchords:
#     fb=value
#     if not hasattr(fb, '__call__'):
#         fbval = fb
#         # print("fbval:\n",fbval)
#         def fbstuff(t, value):
#             return value
#         fb = lambda t: fbstuff(t, fbval)
#         # print("fb:\n",fb)
#     iB0.append(fb)
# iB0 = np.array(iB0)
# # print(iB0)
# external0 = [iB0, None]

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
#eta = np.zeros([numJ,numJ])
#inveta = np.zeros([numJ,numJ])
area = np.zeros(numJ)
invarea = np.zeros(numJ)
for i, Ji in enumerate(Jedges):
    J = my_circuit.edges[Ji]['J']
#    eta[i,i] = J
#    inveta[i,i] = 1/J
    area[i] = J
    invarea[i] = 1/J
#print("eta:\n",eta)

# Capacitors
etaC = np.zeros([numC,numC])
invetaC = np.zeros([numC,numC])
for i, Ci in enumerate(Cedges):
    C = my_circuit.edges[Ci]['C']
    etaC[i,i] = C
    invetaC[i,i] = 1/C
print("etaC:\n",etaC)
#etas = [eta, inveta, etaC, invetaC, area, invarea]
etas = [etaC, invetaC, area, invarea]


L_mats = [L, L_K, L_LK]

# Now adding the model information
models = my_circuit.models


rhs,F_mats_tr,Fbars,Lmats_inv,Ls,model_para = cm.manipulate_matrices(L_mats, F_mats, Rzs, etas, edge_types, num_types, enames, external,models)

# Simulate
timerange = (0,2e-10)

yinit = np.array([0]*(2*numJ+2*numC+numZ))
#yinit = np.array([0]*(2*numJ+2*numC+numZ+numB))

start_time = time.time()
S = scipy.integrate.solve_ivp(lambda t,y:rhs(t,y), timerange, yinit, atol=2.23e-14, rtol=2.23e-14,method = "DOP853")
print("--- %s seconds ---" % (time.time() - start_time))

# Plot

difft = np.diff(S.t)
tperiod = S.t[1:]


for i in range(len(S.y[:])):
#for i in range(1):
    plt.figure(i)
    if i in range(0,numJ):
        plt.plot(S.t[:],S.y[i][:],'-',color='#4258a1',label="ODE Only")
        plt.ylabel("Phase")
        plt.title(r"Junction Phase "+str(i%numJ+1))
    elif i in range(numJ,2*numJ):
        plt.plot(S.t[:],S.y[i][:]*1000,'-',color='#4258a1',label="ODE Only")
        plt.ylabel(r"Voltage ($mV$)")
        plt.title(r"Junction Voltage "+str((i-numJ)%numJ+1))
    elif i in range(2*numJ,2*numJ+numC):
        plt.plot(S.t[:],S.y[i][:],'-',color='#4258a1',label="ODE Only")
        plt.ylabel(r"Phase")
        plt.title(r"Capacitor Phase "+str((i-2*numJ)%numC+1))
    elif i in range(2*numJ+numC,2*numJ+2*numC):
        plt.plot(S.t[:],S.y[i][:]*1000,'-',color='#4258a1',label="ODE Only")
        plt.ylabel(r"Voltage ($mV$)")
        plt.title(r"Capacitor Voltage "+str((i-2*numJ+numC)%numC+1))
    elif i in range(2*numJ+2*numC,2*numJ+2*numC+numZ):
        plt.plot(S.t[:],S.y[i][:],'-',color='#4258a1',label="ODE Only")
        plt.ylabel(r"Flux ($Wb$)")
        plt.title(r"Resistor Flux "+str((i-2*numJ+2*numC)%numZ+1))
    plt.legend(loc="lower right",prop={'size': 10})
    plt.xlabel(r"Time $(s)$")
    plt.savefig("FIG/fig_"+str(i)+".png",dpi=800)
    # plt.show()
    plt.close()

# Post Simulation
start_time = time.time()
transF_JL,transF_JZ,transF_CL,transF_CZ,transF_KL,transF_KZ = F_mats_tr
Fbar_KL,Fbar_JZ,Fbar_CZ,Fbar_JB,Fbar_CB = Fbars
invL,invL_K,invLbar,invLbar_K,invL_LL,invLtwidle_L,invL_D,invLtwidle_LL,invLtwidle_D = Lmats_inv
Lbar,Lbar_K,Ltwidle_K,L_LL,L_LZ,L_ZL,L_D= Ls
i0,Ca = model_para
sim_num = len(S.t)

def renameOutputs(sim_num,num_start,num_end,S):
    if num_end-num_start>0:
        var = np.array(S.y[num_start:num_end])
    else:
        var = np.zeros(sim_num)
    return var

def indexOutputs(arr,i,num):
    if num>0:
        value = arr[:,i]
    else:
        value = 0
    return value

def Ext(lst):
    return [item[i] for item in lst]

def Solvei_LK(L,L_K,L_LK,PhiL,PhiK):
    L_t = np.hstack((L,L_LK))
    L_tp = np.hstack((np.transpose(L_LK),L_K))
    L_t = np.vstack((L_t,L_tp))

    Y = np.concatenate((PhiL,PhiK),axis=None)
    Y = Y.T

    i_LK = np.linalg.solve(L_t, Y)
    iL = i_LK[:numL]
    iK = i_LK[numL:numL+numK]
    return iL, iK

phi_arr = renameOutputs(sim_num,0,numJ,S)
V_arr = renameOutputs(sim_num,numJ,2*numJ,S)
phic_arr = renameOutputs(sim_num,2*numJ,2*numJ+numC,S)
Vc_arr = renameOutputs(sim_num,2*numJ+numC,2*numJ+2*numC,S)
PhiZ_arr = renameOutputs(sim_num,2*numJ+2*numC,2*numJ+2*numC+numZ,S)

phijc_arr = phi_arr
np.append(phijc_arr,phic_arr)

transF_JCL = np.transpose(np.concatenate((F_JL, F_CL),axis=0))
transF_JCZ = np.transpose(np.concatenate((F_JZ, F_CZ),axis=0))

invL_J =  i0

iBvalue_arr = []
Vz_arr = []
PhiL_arr = []
PhiK_arr = []
iK_arr = []
iL_arr = []
ForcingJ_arr =[]
ForcingC_arr= []
sqbrack1_arr= []
sqbrack2_arr= []

for i in range(len(S.t)):
    iBvalue = [f(S.t[i]) for f in iB]

    phijc = indexOutputs(phijc_arr,i,numC+numJ)
    V = indexOutputs(V_arr,i,numJ)
    phic = indexOutputs(phic_arr,i,numC)
    Vc = indexOutputs(Vc_arr,i,numC)
    PhiZ = indexOutputs(PhiZ_arr,i,numZ)

    # NEED to check the dim of each term in the equation

    # sqbrack1
    try:
        sqbrack1a = (Phi0 /2/np.pi) * transF_JCL @ phijc
    except ValueError:
        sqbrack1a = 0
    try:
        sqbrack1b = -transF_KL @ Ltwidle_K @ F_KB @ iBvalue
    except ValueError:
        sqbrack1b = 0
    sqbrack1 = sqbrack1a+sqbrack1b

    # sqbrack2
    try:
        sqbrack2a = (Phi0/2/np.pi) * transF_JCZ @ phijc
    except ValueError:
        sqbrack2a = 0
    try:
        sqbrack2b = - transF_KZ @ Ltwidle_K @ F_KB @ iBvalue
    except ValueError:
        sqbrack2b = 0
    try:
        sqbrack2c = - PhiZ
    except ValueError:
        sqbrack2c = 0
    sqbrack2 = sqbrack2a + sqbrack2b + sqbrack2c

    # Vz
    try:
        Vza = -invL_D @ L_ZL @ invLtwidle_LL @ sqbrack1
    except ValueError:
        Vza = 0
    try:
        Vzb = invLtwidle_D @ sqbrack2
    except ValueError:
        Vzb = 0
    try:
        Vz = Rz @ (Vza+Vzb)
    except ValueError:
        Vz = 0

    # PhiL
    try:
        PhiLa = invLtwidle_LL @ sqbrack1
    except ValueError:
        PhiLa = 0
    try:
        PhiLb = invL_LL @ L_LZ @ (- invLtwidle_D) @ sqbrack2
    except ValueError:
        PhiLb = 0
    try:
        PhiL = Lbar @ (PhiLa + PhiLb)
    except ValueError:
        PhiL = 0

    # ForcingJ
    try:
        ForcingJa = F_JL @ invLtwidle_L @ PhiL
    except ValueError:
        ForcingJa = 0
    try:
        ForcingJb = Fbar_JZ @ invRz @ Vz
    except ValueError:
        ForcingJb = 0
    try:
        ForcingJc = Fbar_JB @ iBvalue
    except ValueError:
        ForcingJc = 0
    ForcingJ = ForcingJa + ForcingJb + ForcingJc

    # ForcingC
    try:
        ForcingCa = F_CL @ invLtwidle_L @ PhiL
    except ValueError:
        ForcingCa = 0
    try:
        ForcingCb = Fbar_CZ @ invRz @ Vz
    except ValueError:
        ForcingCb = 0
    try:
        ForcingCc = Fbar_CB @ iBvalue
    except ValueError:
        ForcingCc = 0
    ForcingC = ForcingCa +ForcingCb +ForcingCc

    # PhiK
    try:
        PhiKa = Fbar_KL @ invLbar @ PhiL
    except ValueError:
        PhiKa = 0
    try:
        PhiKb = F_KZ @ invRz @ Vz
    except ValueError:
        PhiKb = 0
    try:
        PhiKc = F_KB @ iBvalue
    except ValueError:
        PhiKc = 0
    try:
        PhiK = - Ltwidle_K @ (PhiKa+PhiKb+PhiKc)
    except ValueError:
        PhiK = 0

    # Solving iL+iK
    iL,iK = Solvei_LK(L,L_K,L_LK,PhiL,PhiK)

    sqbrack1_arr.append(sqbrack1)
    sqbrack2_arr.append(sqbrack2)
    Vz_arr.append(Vz)
    PhiL_arr.append(PhiL)
    iL_arr.append(iL)
    PhiK_arr.append(PhiK)
    iK_arr.append(iK)
    ForcingJ_arr.append(ForcingJ)
    ForcingC_arr.append(ForcingC)
    iBvalue_arr.append(iBvalue)


print("--- %s seconds ---" % (time.time() - start_time))

# Current sources

for i in range(numB):
    plt.plot(S.t[:],np.array(Ext(iBvalue_arr))*1e6,'-',color='#4258a1',label="ODE Only")
    plt.legend(loc="lower right",prop={'size': 10})
    plt.xlabel(r"Time $(s)$")
    plt.ylabel(r"Current $(\mu A)$")
    plt.title(r"Current Source "+str(i%numB+1))
    plt.savefig("FIG/fig_"+str(2*numJ+2*numC+numZ+i)+".png",dpi=800)
    # plt.show()
    plt.close()

# Tree inductor flux

for i in range(numK):
    plt.plot(S.t[:],np.array(Ext(PhiK_arr)),'-',color='#4258a1',label="ODE Only")
    plt.legend(loc="lower right",prop={'size': 10})
    plt.xlabel(r"Time $(s)$")
    plt.ylabel(r"Flux $(Wb)$")
    plt.title(r"Tree Inductor "+str(i%numK+1))
    plt.savefig("FIG/fig_"+str(2*numJ+2*numC+numZ+numB+i)+".png",dpi=800)
    # plt.show()
    plt.close()

# Tree inductor current

for i in range(numK):
    plt.plot(S.t[:],np.array(Ext(iK_arr))*1e6,'-',color='#4258a1',label="ODE Only")
    plt.legend(loc="lower right",prop={'size': 10})
    plt.xlabel(r"Time $(s)$")
    plt.ylabel(r"Current $(\mu A)$")
    plt.title(r"Tree Inductor "+str(i%numK+1))
    plt.savefig("FIG/fig_"+str(2*numJ+2*numC+numZ+numB+numK+i)+".png",dpi=800)
    # plt.show()
    plt.close()

# Chord inductor flux

for i in range(numL):
    plt.plot(S.t[:],np.array(Ext(PhiL_arr)),'-',color='#4258a1',label="ODE Only")
    plt.legend(loc="lower right",prop={'size': 10})
    plt.xlabel(r"Time $(s)$")
    plt.ylabel(r"Flux $(Wb)$")
    plt.title(r"Chord Inductor "+str(i%numK+1))
    plt.savefig("FIG/fig_"+str(2*numJ+2*numC+numZ+numB+2*numK+i)+".png",dpi=800)
    # plt.show()
    plt.close()

# Chord inductor current

for i in range(numL):
    plt.plot(S.t[:],np.array(Ext(iL_arr))*1e6,'-',color='#4258a1',label="ODE Only")
    plt.legend(loc="lower right",prop={'size': 10})
    plt.xlabel(r"Time $(s)$")
    plt.ylabel(r"Current $(\mu A)$")
    plt.title(r"Chord Inductor "+str(i%numK+1))
    plt.savefig("FIG/fig_"+str(2*numJ+2*numC+numZ+numB+2*numK+numL+i)+".png",dpi=800)
    # plt.show()
    plt.close()

Vz_diff = []
Vz_diff = np.diff(S.y[-1][:]) / np.diff(S.t)
Vz_diff = cm.reshape(Vz_diff)

for i in range(numZ):
    PhiZ_num = np.array(Ext(Vz_arr))[1:] * difft[:]

for i in range(numZ):
    plt.plot(S.t[:],np.array(Ext(Vz_arr))*1e3,'-',color='#4258a1',label="ODE Only")
    plt.plot(S.t[1:],np.array(Vz_diff)*1e3,'--',color='red',label="Numerical Vz")
    plt.legend(loc="lower right",prop={'size': 10})
    plt.xlabel(r"Time $(s)$")
    plt.ylabel(r"Voltage $(mV)$")
    plt.title(r"External Resistor "+str(i%numZ+1))
    plt.savefig("FIG/fig_"+str(2*numJ+2*numC+2*numZ+numB+2*numK+2*numL+i)+".png",dpi=800)
    # plt.show()
    plt.close()
