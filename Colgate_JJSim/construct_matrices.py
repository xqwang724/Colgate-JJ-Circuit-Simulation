import numpy as np
import numpy.linalg as LA
import networkx as nx
import read_cir


def makeF(tree, zedges, chords):
    """Construct a circuit matrix as in burkard, et al. (DiVincenzo)

    Our notation is heavily dependent on that paper.
    Matrices can be F_CL, F_CB, F_CZ;  F_KL, F_KB, F_KZ
    where tree edges can be
          C - junctions (or bare capacitors)
          K - inductors in the tree of the circuit (see findtree.py)
    and chord edges can be
          L - inductors not in the tree
          B - current sources not in the tree
          Z - resistors not in the tree
    Note: edges can be none of these -- they just wont impact the matrix

    Parameters
    ==========
    tree : MultiDiGraph
        The tree for the circuit found using methods in findtree.py
        The edges are directed. We treat it as undirected when finding paths.
    zedges : list of edge tuples
        A list of the tree edges of type we are using (either C or K)
    chords : list of edge tuples
        A list of the chord edges of type we are using (either L, B, or Z)

    Returns
    =======
    F : numpy.array
        A matrix with size len(zedges) by len(chords) and values 1, -1 or 0. F[i,j] = 0 -> tree edge i is not in a loop created by chord j
          F[i,j] = 1 -> tree edge i has same direction as the loop whose
                        direction is given by chord j
          F[i,j] = -1 -> tree edge i has opposite direction as loop whose
                        direction is given by chord j

    Notes
    =====
    It is possible to have no edges of some type in which case
    an empty matrix is returned. Matrices are stored as numpy arrays.
    So, matrix multiplication is done via np.dot(A, B) [or A@B in python3.6+]
    """
    m = len(zedges)
    n = len(chords)
    F = np.zeros((m,n))
    untree = nx.Graph(tree)  # undirected version of the tree

    for j, chord in enumerate(chords):
        u, v = chord[:2]
        path = nx.shortest_path(untree, u, v)
#        print("Chord:\n",chords)
        path_edges = set(nx.utils.pairwise(path))
        # print("EDGES:\n",path_edges)
        for i, e in enumerate(zedges):
            n, nbr = e[:2]
            if (n, nbr) in path_edges:
                F[i,j] = 1
            elif (nbr,n) in path_edges:
                F[i,j] = -1
            # else F[i,j] is zero
    return F

def reshape(ar):
        shape = np.shape(ar)[0]
        if len(np.shape(ar)) == 1:
            ar = np.reshape(ar,(shape,1))
        else:
            pass
        return ar

def make_matrices(G, tree):
    current_source_edges = [(u,v,k) for (u,v,k,dd) in G.edges(data=True, keys=True) if 'I' in dd]
    voltage_source_edges = [(u,v,k) for (u,v,k,dd) in G.edges(data=True, keys=True) if 'V' in dd]
    capacitor_edges = [(u,v,k) for (u,v,k,dd) in G.edges(data=True, keys=True) if 'C' in dd]
    capacitor_nodes = [(u,v) for (u,v,k,dd) in G.edges(data=True, keys=True) if 'C' in dd]
    junction_edges = [(u,v,k) for (u,v,k,dd) in G.edges(data=True, keys=True) if 'J' in dd]
    junction_nodes = [(u,v) for (u,v,k,dd) in G.edges(data=True,keys=True) if 'J' in dd]

    # print(junction_edges)
    resistor_edges = [(u,v,k) for (u,v,k,dd) in G.edges(data=True, keys=True) if ('R' in dd) and (u,v) not in junction_nodes and (u,v) not in capacitor_nodes]
#    resistor_edges = [(u,v,k) for (u,v,k,dd) in G.edges(data=True, keys=True) if ('R' in dd)]
    shunt_resistor_edges = [(u,v,k) for (u,v,k,dd) in G.edges(data=True,keys=True) if ('R' in dd) and (u,v) in junction_nodes]
#    shunt_resistor_edges = []
    shunt_resistor_edges_C = [(u,v,k) for (u,v,k,dd) in G.edges(data=True,keys=True) if ('R' in dd) and (u,v) in capacitor_nodes]
#    shunt_resistor_edges_C = []
    # print("Junction Edges", junction_edges)
    # print("Resistor Edges", resistor_edges)
#    print("shunt resistor edges", shunt_resistor_edges)
#    print("Junction Nodes", junction_nodes)
    inductor_edges = [(u,v,k) for (u,v,k,dd) in G.edges(data=True, keys=True) if 'L' in dd]

#    chords = [(u,v,k) for (u,v,k) in G.edges(keys=True) if not tree.has_edge(u,v,k)]
    Lchords = [(u,v,k) for (u,v,k) in inductor_edges if not tree.has_edge(u,v,k)]
    Bchords = [(u,v,k) for (u,v,k) in current_source_edges if not tree.has_edge(u,v,k)]
    Zchords = [(u,v,k) for (u,v,k) in resistor_edges if not tree.has_edge(u,v,k)]

    Kedges = [(u,v,k) for (u,v,k) in inductor_edges if tree.has_edge(u,v,k)]
    Jedges = [(u,v,k) for (u,v,k) in junction_edges if tree.has_edge(u,v,k)]
    Cedges = [(u,v,k) for (u,v,k) in capacitor_edges if tree.has_edge(u,v,k)]
    ##  TODO: dont we need capacitors separate from JJs at some point?
    Redges = [(u,v,k) for (u,v,k) in shunt_resistor_edges if not tree.has_edge(u,v,k)]
    RedgesC = [(u,v,k) for (u,v,k) in shunt_resistor_edges_C if not tree.has_edge(u,v,k)]


    numJ = len(Jedges)
    numC = len(Cedges)
    numK = len(Kedges)
    numL = len(Lchords)
    numB = len(Bchords)
    numZ = len(Zchords)
    numR = len(Redges)
    numRC = len(RedgesC)

#    print("Jedges",Jedges)
#    print("Cedges",Cedges)
#    print("Kedges",Kedges)
#    print("Lchords",Lchords)
#    print("Bchords",Bchords)
#    print("Zchords",Zchords)
#    print("Redges", Redges)
#    print("RedgesC", RedgesC)

    print("Find Tree:  below is a list of the chords used for each F matrix")
    if Lchords:
#        print("JL", end='')
        F_JL = makeF(tree, Jedges, Lchords)
#        print("CL", end='')
        F_CL = makeF(tree, Cedges, Lchords)
#        print("KL", end='')
        F_KL = makeF(tree, Kedges, Lchords)
    else:
        F_JL = np.zeros((numJ,0))
        F_CL = np.zeros((numC,0))
        F_KL = np.zeros((numK,0))
    if Bchords:
#        print("JB", end='')
        F_JB = makeF(tree, Jedges, Bchords)
        # print(F_JB)
#        print("CB", end='')
        F_CB = makeF(tree, Cedges, Bchords)
#        print("F_CB", F_CB)
#        print("KB", end='')
        F_KB = makeF(tree, Kedges, Bchords)
    else:
        F_JB = np.zeros((numJ,0))
        F_CB = np.zeros((numC,0))
        F_KB = np.zeros((numK,0))
    if Zchords:
#        print("JZ", end='')
        F_JZ = makeF(tree, Jedges, Zchords)
#        print("CZ", end='')
        F_CZ = makeF(tree, Cedges, Zchords)
#        print("KZ", end='')
        F_KZ = makeF(tree, Kedges, Zchords)
    else:
        F_JZ = np.zeros((numJ,0))
        F_CZ = np.zeros((numC,0))
        F_KZ = np.zeros((numK,0))


    F_mats = [F_JL, F_JB, F_JZ, F_CL, F_CB, F_CZ, F_KL, F_KB, F_KZ]
#    print("F_KB", F_KB)
#    print("F_JB", F_JB)
#    print("F_KL", F_KL)
#    print("F_KZ", F_KZ)
    edge_types = [Jedges, Cedges, Kedges, Lchords, Bchords, Zchords, Redges, RedgesC]
    num_types = [numJ, numC, numK, numL, numB, numZ, numR, numRC]

#    # external flux and bias currents
#    phix = []
#    for e in Lchords:
#        # NOT CORRECT WAY TO STORE PHI: fb=G.edges[e]['L']
#        if not hasattr(fb, '__call__'):
#            fb = fb(t)
#        phix.append(fb)
#    phix = np.array(phix)
#    dphix = 0 * phix

    iB = []
    for e in Bchords:
        fb=G.edges[e]['I']
        if not hasattr(fb, '__call__'):
            fbval = fb
            # print("fbval:\n",fbval)
            def fbstuff(t, value):
                return value
            fb = lambda t: fbstuff(t, fbval)
            print("fb:\n",fb)
        iB.append(fb)
    iB = np.array(iB)

    external = [iB, None]
    return [F_mats, edge_types, num_types, external]

def Ifinv(mat):
    try:
        return LA.inv(mat)
    except np.linalg.LinAlgError:
        # return np.divide(mat.shape[:])
        # np.divide(np.ones(mat.shape[:]), mat, matmod)
        matmod = np.divide(np.ones(np.shape(mat)), mat, out=np.zeros_like(mat), where=mat!=0)
        print("FOUND A SINGULAR MATRIX!")
        return matmod

def readRshunt(numjunc,models,Jedges,Redges,Rshunt,invRshunt,area):
    # R = np.zeros([numjunc,numjunc])
    # invR =np.zeros([numjunc,numjunc])
    R = np.zeros((numjunc,1))
    invR = np.zeros((numjunc,1))
    # print("R:", R[1,0])
    # print("Rshunt:\n",Rshunt)
    # print("invRshunt:\n",invRshunt)

    for i in range(numjunc):
        jnode = (Jedges[i][0],Jedges[i][1])
        rnodes = [(u,v) for (u,v,k) in Redges]
        # print("RNODES", rnodes)
        # print("JNODE", jnode)
        # print("Rshunt",Rshunt[i])
        # print(models["B"+str(i)]["rtype"])
        if jnode in rnodes and models["B"+str(i)]["rtype"] == 1:
            # R[i][i] = 1/(1/np.sqrt(area)*models["B"+str(i)]["rn"] + invRshunt)
            # R.append(1/(1/np.sqrt(area)*models["B"+str(i)]["rn"] + invRshunt[i]))
            R[i,0] = (1/(1/np.sqrt(area)*models["B"+str(i)]["rn"] + invRshunt[i]))
            # invR[i][i] = (1/np.sqrt(area)*models["B"+str(i)]["rn"] + invRshunt)
            invR[i,0] = ((1/np.sqrt(area)*models["B"+str(i)]["rn"] + invRshunt[i]))
        elif jnode in rnodes and models["B"+str(i)]["rtype"] == 0:
            # R[i][i] =Rshunt
            R[i,0] = Rshunt[i]
            # invR[i][i] =1/Rshunt
            invR[i,0]=invRshunt[i]
#        elif jnode not in rnodes and models["B"+str(i)]["rtype"] == 1:
#            print("PRINT:",type(1/np.sqrt(area)*models["B"+str(i)]["rn"]))
#            print(1/np.sqrt(area)*models["B"+str(i)]["rn"])
#            R[i,0] = 1/np.sqrt(area)*models["B"+str(i)]["rn"]
        else:
            print("No shunt resistor found.")
            continue


    R = reshape(np.transpose(np.array(R)))[0]
    invR = reshape(np.transpose(np.array(invR)))[0]
    print("R: \n",R)
    return R,invR


def readRshuntC(numc,Cedges,RedgesC,Rshunt,invRshunt):
    # R = np.zeros([numjunc,numjunc])
    # invR =np.zeros([numjunc,numjunc])
    R = np.zeros((numc,1))
    invR = np.zeros((numc,1))
    # print("R:", R[1,0])
    # print("Rshunt:\n",Rshunt)
    # print("invRshunt:\n",invRshunt)

    for i in range(numc):
        cnode = (Cedges[i][0],Cedges[i][1])
        rnodes = [(u,v) for (u,v,k) in RedgesC]
        # print("RNODES", rnodes)
        # print("JNODE", jnode)
        # print("Rshunt",Rshunt[i])
        # print(models["B"+str(i)]["rtype"])
        if cnode in rnodes:
            R[i,0] = Rshunt[i]
            invR[i,0] = invRshunt[i]
        else:
            R[i,0] = 0
            invR[i,0] = 0
            print("No shunt for capacitor found.")
            continue

    R = reshape(np.transpose(np.array(R)))[0]
    invR = reshape(np.transpose(np.array(invR)))[0]
#    print("R: \n",R)
    return R,invR


# def manipulate_matrices(L_mats, F_mats, Rzs, etas, edge_types, num_types, enames, external,models):
def manipulate_matrices(L_mats, F_mats, Rzs, etas, edge_types, num_types, enames, external,models):
    L, L_K, L_LK = L_mats
    F_JL, F_JB, F_JZ, F_CL, F_CB, F_CZ, F_KL, F_KB, F_KZ = F_mats
    Jedges, Cedges, Kedges, Lchords, Bchords, Zchords, Redges, RedgesC = edge_types
    numJ, numC, numK, numL, numB, numZ, numR, numRC = num_types
    iB, diB = external

#    print("iB, diB:", iB, diB)
    iBvalues = [f(0) for f in iB]
    # print("iBvalues:",iBvalues)
    # print("function?",iB[1](0))

    # iBint = []
    # num_steps = 1000
    # for j in range(numB):
    #     iBint.append([])
    #     for i in range(num_steps):
    #         # print("iBint:",iBint[j])
    #         # print("iBvalues:",iBvalues[j])
    #         iBint[j].append(i*iBvalues[j]/num_steps)
    # iBint = np.array(iBint)
    # iB = np.concatenate(iB,iBint)

    # for i in range(len(iBvalues)):
    #     iBint.append(read_cir.make_pulse("pulse(0,"+str(iBvalues[i])+",0,"+str(rise)+",0,0,0)"))
        # print("i:",i)
    # print("IBint: \n",IBint)


    Rz, invRz, Rshunt, invRshunt , RshuntC, invRshuntC = Rzs
    eta, inveta, etaC, invetaC, area, invarea= etas
    #invetaC = reshape(np.transpose(np.array(invetaC)))[0]
    #Rz = reshape(np.transpose(np.array(Rz)))[0]
    #invRz = reshape(np.transpose(np.array(invRz)))[0]
    # gamma, coupling_lambda = params
    # print("invRz:",invRz)

    nume=len(next(iter(enames)))
#    print('Edges that are Chords:')
    for i,e in enumerate(Lchords):
        print('self-inductance for ', enames[e[:nume]],': ', L[i,i])
    for i,ei in enumerate(Lchords):
        for j,ej in enumerate(Lchords):
            if i!=j and L[i,j] != 0:
                print('inductance between ', enames[ei[:nume]], ' and ', enames[ej[:3]], ": ", L[i,j])

    print('Edges that are in the tree:')
    for i,e in enumerate(Kedges):
        print('self-inductance for ', enames[e[:nume]],': ', L_K[i,i])
    for i,ei in enumerate(Kedges):
        for j,ej in enumerate(Kedges):
            if i!=j and L_K[i,j] != 0:
                print('inductance between ', enames[ei[:nume]], ' and ', enames[ej[:nume]], ": ", L_K[i,j])

    print('Inductance from chords to tree inductors')
    for i,ei in enumerate(Kedges):
        for j,ej in enumerate(Lchords):
            print(j,i,ej,ei,L_LK)
            if L_LK[j,i] != 0:
                print('inductance between ', enames[ei[:nume]], ' and ', enames[ej[:nume]], ": ", L_LK[i,j])


    #1
    transF_JL = np.transpose(F_JL)
#    print("transF_JL",transF_JL)
    transF_JZ = np.transpose(F_JZ)
#    print("transF_JZ", transF_JZ)
    transF_CL = np.transpose(F_CL)
#    print("transF_CL",transF_CL)
    transF_CZ = np.transpose(F_CZ)
#    print("transF_CZ", transF_CZ)
    transF_KL = np.transpose(F_KL)
#    print("transF_KL",transF_KL)
    transF_KZ = np.transpose(F_KZ)
#    print("transF_KZ", transF_KZ)

#    invL = LA.inv(L)
    invL = Ifinv(L)
#    print("invL", invL)
#    invL_K = LA.inv(L_K)
    invL_K = Ifinv(L_K)
#    print("invL_K", invL_K)
    transL_LK = np.transpose(L_LK)
#    print("transL_LK", transL_LK)

    #3
    Lbar = L - (L_LK @ invL_K @ transL_LK)
#    print("Lbar", Lbar)
#    invLbar = LA.inv(Lbar)
    invLbar = Ifinv(Lbar)
#    print("invLbar", invLbar)

    #4
    Fbar_KL = F_KL - (invL_K @ transL_LK)
#    print("Fbar_KL", Fbar_KL)
#    transFbar_KL = np.transpose(Fbar_KL)
#    print("transnFbar_KL", transFbar_KL)

    #5
    Lbar_K = L_K - (transL_LK @ invL @ L_LK)
#    print("Lbar_K", Lbar_K)
#    invLbar_K = LA.inv(Lbar_K)
    invLbar_K = Ifinv(Lbar_K)
#    print("invLbar_K", invLbar_K)

    #6
    L_Kshape = list(L_K.shape)
#    Ltwidle_K = (LA.inv(np.eye(L_Kshape[0]) - (((L_K)@(Fbar_KL)@(invL))@(L_LK)@(invLbar_K)) )) @ L_K
    Ltwidle_K = (Ifinv(np.eye(L_Kshape[0]) - (((L_K)@(Fbar_KL)@(invL))@(L_LK)@(invLbar_K)) )) @ L_K
#    print("Ltwidle_K", Ltwidle_K)
#    transLtwidle_K = np.transpose(Ltwidle_K)
#    print("transLtwidle_K", transLtwidle_K)

    #7
    L_LKshape = list(L_LK.shape)
    invLtwidle_L = (np.eye(L_LKshape[0])+invL@(L_LK@invLbar_K@Ltwidle_K@Fbar_KL))@(invLbar)
#    print("invLtwidle_L", invLtwidle_L)

    #12
    L_LL = Lbar + transF_KL @ Ltwidle_K @ Fbar_KL
#    print("L_LL", L_LL)

    #8
#    invL_LL = LA.inv(L_LL)
    invL_LL = Ifinv(L_LL)
#    print("invL_LL", invL_LL)
#    transInvL_LL = np.transpose(invL_LL)
#    print("transInvL_LL", transInvL_LL)

    #9
#    print(F_JZ, F_JL, F_KZ)
    Fbar_JZ = F_JZ + F_JL @ invL @ L_LK @ invLbar_K @ Ltwidle_K @ F_KZ
#    print("Fbar_JZ", Fbar_JZ)
    Fbar_CZ = F_CZ + F_CL @ invL @ L_LK @ invLbar_K @ Ltwidle_K @ F_KZ
#    print("Fbar_CZ", Fbar_CZ)

    #10
    Fbar_JB = F_JB + F_JL @ invL @ L_LK @ invLbar_K @ Ltwidle_K @ F_KB
#    print("Fbar_JB", Fbar_JB)
    Fbar_CB = F_CB + F_CL @ invL @ L_LK @ invLbar_K @ Ltwidle_K @ F_KB
#    print("Fbar_CB", Fbar_CB)

    #14
    L_LZ = transF_KL @ Ltwidle_K @ F_KZ
#    print("L_LZ", L_LZ)

    #15
    L_ZL = transF_KZ @ Ltwidle_K @ Fbar_KL
#    print("L_ZL", L_ZL)

    #extra eqn: D = transF_KZ @ Ltwidle_K @ F_KZ
    L_Z = Rz/1e15
    L_D = transF_KZ @ Ltwidle_K @ F_KZ  #+ L_Z
#    L_D = transF_KZ @ Ltwidle_K @ F_KZ  + L_Z
    #print("L_D", L_D)
#    invL_D = LA.inv(L_D)
    invL_D = Ifinv(L_D)
#    print("invL_D", invL_D)


    # print("Print Equations?:   dv First equation")
    # print("F_JL @ invLtwidle_L @ phiL",F_JL[0,:],"@",invLtwidle_L,"@phiL  = ",F_JL@invLtwidle_L,"@phiL")
    # print("Fbar_JB @ iBvalue",Fbar_JB[0,:],"@iBvalue")
    # print("second equation")
    # print("F_JL @ invLtwidle_L @ phiL",F_JL[1,:],"@",invLtwidle_L,"@phiL  = ",F_JL@invLtwidle_L,"@phiL")
    # print("Fbar_JB @ iBvalue",Fbar_JB[1,:],"@iBvalue")

#    print("F_JL @ invLtwidle_L @ phiL",F_JL,"@",invLtwidle_L,"@phiL  = ",F_JL@invLtwidle_L,"@phiL")
#    print("Fbar_JB @ iBvalue",Fbar_JB,"@iBvalue")

#    invLtwidle_LL = LA.inv(L_LL - L_LZ @ invL_D @ L_ZL)
    invLtwidle_LL = Ifinv(L_LL - L_LZ @ invL_D @ L_ZL)
    print("invLtwidle_LL",invLtwidle_LL)
#    invLtwidle_D = LA.inv(L_D - L_ZL @ invL_LL @ L_LZ)
    invLtwidle_D = Ifinv(L_D - L_ZL @ invL_LL @ L_LZ)
#    print("invLtwidle_D",invLtwidle_D)
#    print("Lbar",Lbar)

#    print("Ltwidle_K",Ltwidle_K)

    Phi0 = 2.067833848e-15
    Req,invReq = readRshunt(numJ,models,Jedges,Redges,Rshunt,invRshunt,area)
    ReqC, invReqC = readRshuntC(numC,Cedges,RedgesC,RshuntC,invRshuntC)
    i0 = np.zeros((numJ,1))
    Ca = np.zeros((numJ,1))

    for i in range(numJ):
#        print(models["B"+str(i)]["icrit"])
        i0[i,0] = area[i]*models["B"+str(i)]["icrit"]
#    for i in range(numJ,numJ+numC):
#        i0[i,0] = 0
    i0 = reshape(np.transpose(np.array(i0)))[0]
    #print("i0",i0)

    for i in range(numJ):
        Ca[i,0] = area[i]*models["B"+str(i)]["cap"]
    Ca = reshape(np.transpose(np.array(Ca)))[0]
    #print("Ca", Ca)

    # def rhs(t,y):
    #     phi = y[:numJ]
    #     v = y[numJ:(2*numJ)]
    #     vcap = y[(2*numJ):(2*numJ+numC)]
    #     phiZ = y[(2*numJ+numC):(2*numJ+numC+numZ)]
    #     #
    #     iBvalue = np.array([f(t) for f in iB])
    #     #
    #     # Note: sqbrack would have phi_x in there is we set up a way to store it.
    #     lam = coupling_lambda
    #     sqbrack1 = transF_JL @ phi - transF_KL @ Ltwidle_K @ F_KB @ iBvalue/lam
    #     sqbrack2 = transF_JZ @ v - transF_KZ @ Ltwidle_K @ F_KB @ iBvalue/lam - phiZ
    #     vz = lam/gamma * Rz @ (- invLtwidle_LL @ sqbrack1 + invLtwidle_D @ sqbrack2)
    #     phiL = Lbar @ (invLtwidle_LL @ sqbrack1 + invL_LL @ L_LZ @ (- invLtwidle_D) @ sqbrack2)
    #     #
    #     ForcingJ = inveta @ (lam * F_JL @ invLtwidle_L @ phiL + gamma * Fbar_JZ @ invRz @ vz + Fbar_JB @ iBvalue)
    #     ForcingC = invetaC @ (lam * F_CL @ invLtwidle_L @ phiL + gamma * Fbar_CZ @ invRz @ vz + Fbar_CB @ iBvalue)
    #     #
    #     dphi = v
    #     dv = -np.sin(phi) - gamma * v - ForcingJ
    #     dvc = -ForcingC
    #     dphiz = vz
    #     result_rhs = np.concatenate((dphi, dv, dvc, dphiz))
    #     return result_rhs
    # return rhs

    def rhs(t,y):
        phi = y[:numJ]
        # v = y[numJ:(2*numJ)]
        V = y[numJ:(2*numJ)]
        phic = y[(2*numJ):(2*numJ+numC)]
        Vc = y[(2*numJ+numC):(2*numJ+2*numC)]
        PhiZ = y[(2*numJ+2*numC):(2*numJ+2*numC+numZ)]
#        PhiZ = y[(2*numJ):(2*numJ+numZ)]

        iBvalue = np.array([f(t) for f in iB])
        invL_J =  i0
    #Correct but partial
#        sqbrack1 = (Phi0 /2/np.pi) * transF_JL @ phi - transF_KL @ Ltwidle_K @ F_KB @ iBvalue
#        sqbrack2 = (Phi0/2/np.pi) * transF_JZ @ phi - transF_KZ @ Ltwidle_K @ F_KB @ iBvalue - PhiZ

    #V2 (Wrong)
#        sqbrack1 = (Phi0 /2/np.pi) * (transF_JL @ phi + transF_CL @ phic) - transF_KL @ Ltwidle_K @ F_KB @ iBvalue
#        sqbrack2 = (Phi0/2/np.pi) * (transF_JZ @ phi + transF_CZ @ phic) - transF_KZ @ Ltwidle_K @ F_KB @ iBvalue - PhiZ
        #print(transF_CZ @ phic)
    #V1 Merge vector both
        phijc = np.concatenate((phi,phic))
        transF_JCL = np.transpose(np.concatenate((F_JL, F_CL),axis=0))
        transF_JCZ = np.transpose(np.concatenate((F_JZ, F_CZ),axis=0))
        sqbrack1 = (Phi0 /2/np.pi) * transF_JCL @ phijc - transF_KL @ Ltwidle_K @ F_KB @ iBvalue
        sqbrack2 = (Phi0/2/np.pi) * transF_JCZ @ phijc - transF_KZ @ Ltwidle_K @ F_KB @ iBvalue - PhiZ
#        print(sqbrack2)
#        print(phijc)

        Vz = Rz @ (- invL_D @ L_ZL @ invLtwidle_LL @ sqbrack1 + invLtwidle_D @ sqbrack2)
        #print(Vz)
        # WE CAN'T SIMULATE A CIRCUIT WITHOUT AN INDUCTOR?? That's Because
        # L_Zis taken out


        PhiL = Lbar @ (invLtwidle_LL @ sqbrack1 + invL_LL @ L_LZ @ (- invLtwidle_D) @ sqbrack2)
        ForcingJ = (F_JL @ invLtwidle_L @ PhiL + Fbar_JZ @ invRz @ Vz + Fbar_JB @ iBvalue)

    #Wrong Equations Right Result for RC
        #ForcingC = (F_CL @ invLtwidle_L @ PhiL + Fbar_CZ @ invRz @ Vc + Fbar_CB @ iBvalue)

        ForcingC = (F_CL @ invLtwidle_L @ PhiL + Fbar_CZ @ invRz @ Vz + Fbar_CB @ iBvalue)
#        print(ForcingC)
#        print(np.shape(Fbar_JB))
#        print(np.shape(Fbar_CB))
#        print(np.shape(Vc))
#        print(np.shape(Vz))
        #ForcingC =  Fbar_CZ @ invRz @ Vc

        dphi = 2*np.pi/Phi0*V
        dV = 1/Ca * (-invL_J * np.sin(phi) - invReq * V - ForcingJ)
        dphic = 2*np.pi/Phi0*Vc
        #print(invetaC)
#        ForcingC = (F_CL @ invLtwidle_L @ PhiL + Fbar_CZ @ invRz @ Vz + Fbar_CB @ iBvalue)
        dVc = - invetaC @ (invReqC * Vc + ForcingC)
        dPhiz = Vz
#        result_rhs = np.concatenate((dphi, dV,dPhiz))
        result_rhs = np.concatenate((dphi, dV, dphic, dVc, dPhiz))
        return result_rhs
    return rhs
