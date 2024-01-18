import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import igraph as ig  # may have to pip install
import itertools
from tqdm import tqdm


def type_mat_rearrange(input_mat):
    typee = type(input_mat[0,0])
    rearranged = np.zeros(input_mat.shape).astype(typee)
    if input_mat.shape[0]==4:
        rearranged[3,3]=input_mat[0,0]; rearranged[3,0]=input_mat[0,1]
        rearranged[3,1]=input_mat[0,2]; rearranged[3,2]=input_mat[0,3]
        rearranged[0,3]=input_mat[1,0]; rearranged[0,0]=input_mat[1,1]
        rearranged[0,1]=input_mat[1,2]; rearranged[0,2]=input_mat[1,3]
        rearranged[1,3]=input_mat[2,0]; rearranged[1,0]=input_mat[2,1]
        rearranged[1,1]=input_mat[2,2]; rearranged[1,2]=input_mat[2,3]
        rearranged[2,3]=input_mat[3,0]; rearranged[2,0]=input_mat[3,1]
        rearranged[2,1]=input_mat[3,2]; rearranged[2,2]=input_mat[3,3]
    else:
        rearranged[4,4]=input_mat[0,0]; rearranged[4,0]=input_mat[0,1]; rearranged[4,2]=input_mat[0,2]
        rearranged[4,3]=input_mat[0,3]; rearranged[4,1]=input_mat[0,4]
        rearranged[0,4]=input_mat[1,0]; rearranged[0,0]=input_mat[1,1]; rearranged[0,2]=input_mat[1,2]
        rearranged[0,3]=input_mat[1,3]; rearranged[0,1]=input_mat[1,4]
        rearranged[2,4]=input_mat[2,0]; rearranged[2,0]=input_mat[2,1]; rearranged[2,2]=input_mat[2,2]
        rearranged[2,3]=input_mat[2,3]; rearranged[2,1]=input_mat[2,4]
        rearranged[3,4]=input_mat[3,0]; rearranged[3,0]=input_mat[3,1]; rearranged[3,2]=input_mat[3,2]
        rearranged[3,3]=input_mat[3,3]; rearranged[3,1]=input_mat[3,4]
        rearranged[1,4]=input_mat[4,0]; rearranged[1,0]=input_mat[4,1]; rearranged[1,2]=input_mat[4,2]
        rearranged[1,3]=input_mat[4,3]; rearranged[1,1]=input_mat[4,4]
    return rearranged


def all_nonzero(arr1, arr2, arr3, arr4):  # has to be same length
    arr1_bi = np.copy(arr1); arr1_bi[arr1_bi!=0] = 1
    arr2_bi = np.copy(arr2); arr2_bi[arr2_bi!=0] = 1
    arr3_bi = np.copy(arr3); arr3_bi[arr3_bi!=0] = 1
    arr4_bi = np.copy(arr4); arr4_bi[arr4_bi!=0] = 1
    arr1_new = []; arr2_new = []; arr3_new = []; arr4_new = []
    
    for ii in range(len(arr1)):
        if arr1_bi[ii]+arr2_bi[ii]+arr3_bi[ii]+arr4_bi[ii]>=4:  # change
            arr1_new.append(arr1[ii]); arr2_new.append(arr2[ii])
            arr3_new.append(arr3[ii]); arr4_new.append(arr4[ii])
        else:
            continue
    
    arr1_new = np.array(arr1_new); arr2_new = np.array(arr2_new)
    arr3_new = np.array(arr3_new); arr4_new = np.array(arr4_new)
    return arr1_new, arr2_new, arr3_new, arr4_new


def sort_network_dict(network_dict):
    sorted_list = sorted(network_dict.items(), key = lambda item: item[1], reverse = True)
    return sorted_list


def to_strid(list_of_tuple, com_id):
    if isinstance(list_of_tuple[0][0],int)==True:
        tuplelen = 1
    else:
        tuplelen = 2
    strid = np.zeros((len(list_of_tuple),tuplelen)).astype(str)
    vals = np.zeros((len(list_of_tuple),1))
    if tuplelen == 1:
        for ii in range(len(list_of_tuple)):
            temp = com_id[ np.where(com_id[:,1] == (list_of_tuple[ii][0])+1)[0], 0 ]
            strid[ii] = temp
            vals[ii] = list_of_tuple[ii][1]
        new_list = np.hstack((strid,vals))
    else:
        for ii in range(com_id.shape[0]):
            strid[ii,0] = com_id[ int(np.where(com_id[:,1] == (list_of_tuple[ii][0][0])+1)[0]), 0 ]
            strid[ii,1] = com_id[ int(np.where(com_id[:,1] == (list_of_tuple[ii][0][1])+1)[0]), 0 ]
            vals[ii] = list_of_tuple[ii][1]
        new_list = np.hstack((strid,vals))
    return new_list


def module_sim(list1, list2):
    if len(list1)==0 or len(list2)==0:
        raise("empty array in the input.")
    else:
        numer = np.intersect1d(list1, list2, assume_unique=False)
        denom = np.union1d(list1, list2)
        Similarity = numer.size / denom.size
    
    return Similarity


def module_sim2(list1, list2):
    if len(list1)==0 or len(list2)==0:
        raise("empty array in the input.")
    else:
        numer = np.intersect1d(list1, list2, assume_unique=False)
        denom1 = len(list1); denom2 = len(list2)
        Sim1 = numer.size / denom1; Sim2 = numer.size / denom2
    
    return Sim1, Sim2


def config_model(original_matrix, rewire_num):
    import random
    cfg_model = np.copy(original_matrix)
    
    for ii in range(rewire_num):
        nz_list = list(zip(np.nonzero(cfg_model)[0],np.nonzero(cfg_model)[1]))  # list of tuple
        rand_pick = random.sample(range(len(nz_list)),2)
        pick1 = nz_list[rand_pick[0]]
        pick2 = nz_list[rand_pick[1]]
        
        # check
        check = True
        while check:
            if (cfg_model[pick1[0],pick2[1]]==0 and cfg_model[pick2[0],pick1[1]]==0):
                check = False
            else:
                rand_pick = random.sample(range(len(nz_list)),2)
                pick1 = nz_list[rand_pick[0]]
                pick2 = nz_list[rand_pick[1]]

        # rewire
        temp1=cfg_model[pick1]; temp2=cfg_model[pick2]
        cfg_model[pick1]=0; cfg_model[pick2]=0

        cfg_model[pick1[0],pick2[1]]=temp1; cfg_model[pick2[0],pick1[1]]=temp2
        
    return cfg_model


def twocell_motif(Adj_matrix):    # assumes row = pre // col = post
    mat_nodes = np.arange(Adj_matrix.shape[0])
    pick2cell = list(itertools.combinations(mat_nodes, 2))
    uncon = 0; uni = 0; bi = 0
    
    for ii in range(len(pick2cell)):
        pick = list(pick2cell[ii])
        submat = Adj_matrix[np.ix_(pick,pick)]
        if submat[0,1]==0 and submat[1,0]==0:
            uncon += 1
        elif (submat[0,1]==0 and submat[1,0]!=0) or (submat[0,1]!=0 and submat[1,0]==0):
            uni += 1
        elif submat[0,1]!=0 and submat[1,0]!=0:
            bi += 1
    return (uncon, uni, bi)


def find_twocell_motif(Adj_matrix, pattern_number):    # refer to Song et al. 2005
    mat_nodes = np.arange(Adj_matrix.shape[0])
    pick2cell = list(itertools.combinations(mat_nodes, 2))
    mask = np.array([[0,1],[1,0]])
    motif_list = []
    weights = []
    for ii in range(len(pick2cell)):
        pick = list(pick2cell[ii])
        submat = Adj_matrix[np.ix_(pick,pick)]
        submat = submat*mask
        w1 = submat[0,1]; w2 = submat[1,0]
        if w1==0 and w2==0:
            pat = 0    # unconnected
        elif (w1==0 and w2!=0):
            pat = 1    # uni-
            pick = copy.deepcopy(pick[::-1])
        elif (w1!=0 and w2==0):
            pat = 1
        else:
            pat = 2    # bi-
        
        if pat == pattern_number:
            motif_list.append(pick)
            ww = list(submat[np.nonzero(submat)])
            weights.append(ww)
        else:
            continue

    return motif_list, weights


def difftype_er(Adj_matrix, type1_idx, type2_idx):
    submat1 = Adj_matrix[np.ix_(type1_idx, type2_idx)]
    submat2 = Adj_matrix[np.ix_(type2_idx, type1_idx)]
    linknum = np.nonzero(submat1)[0].size + np.nonzero(submat2)[0].size 
    ermat = np.zeros((type1_idx.size,type2_idx.size,2))  # 3d array
    idx_list = list(itertools.product(np.arange(type1_idx.size), np.arange(type2_idx.size), np.arange(2)))
        # all possible indices
    pick_idx = np.random.permutation(len(idx_list))[:linknum].astype(int)

    for ii in range(pick_idx.size):
        ermat[idx_list[pick_idx[ii]]] = 1

    motif_cnt = np.zeros(3).astype(int)  # unc-, uni-, bi-
    for ii in range(type1_idx.size):
        for jj in range(type2_idx.size):
            if (ermat[ii,jj,0]==0) and (ermat[ii,jj,1]==0):
                motif_cnt[0] += 1
            elif (ermat[ii,jj,0]==0 and ermat[ii,jj,1]!=0) or (ermat[ii,jj,0]!=0 and ermat[ii,jj,1]==0):
                motif_cnt[1] += 1
            elif (ermat[ii,jj,0]!=0) and (ermat[ii,jj,1]!=0):
                motif_cnt[2] += 1
 
    return motif_cnt


def connection_prob(Adj_matrix, type1_idx, type2_idx):
    if np.array_equal(type1_idx, type2_idx)==True:
        submat = Adj_matrix[np.ix_(type1_idx, type2_idx)]
        linknum = np.nonzero(submat)[0].size
        denom = type1_idx.size*(type2_idx.size-1)
    else:
        submat1 = Adj_matrix[np.ix_(type1_idx, type2_idx)]; submat2 = Adj_matrix[np.ix_(type2_idx, type1_idx)]
        linknum = np.nonzero(submat1)[0].size + np.nonzero(submat2)[0].size
        denom = type1_idx.size*type2_idx.size*2
    prob = linknum / denom
    
    theoretic_cnt = np.zeros(3)
    theoretic_cnt[0] = (denom/2) * (1-prob)**2  # unc-
    theoretic_cnt[1] = (denom/2) * 2 * prob * (1-prob)  # uni-
    theoretic_cnt[2] = (denom/2) * prob**2  # bi-
    
    return prob, theoretic_cnt


def twop_to_threep(two_p):
    three_p = np.zeros(16)
    three_p[0] = 1*two_p[0]*two_p[0]*two_p[0]; three_p[1] = 6*two_p[0]*two_p[0]*two_p[1]
    three_p[2] = 3*two_p[0]*two_p[0]*two_p[2]; three_p[3] = 3*two_p[0]*two_p[1]*two_p[1]
    three_p[4] = 3*two_p[0]*two_p[1]*two_p[1]; three_p[5] = 6*two_p[0]*two_p[1]*two_p[1]
    three_p[6] = 6*two_p[0]*two_p[1]*two_p[2]; three_p[7] = 6*two_p[0]*two_p[1]*two_p[2]
    three_p[8] = 3*two_p[0]*two_p[2]*two_p[2]; three_p[9] = 6*two_p[1]*two_p[1]*two_p[1]
    three_p[10] = 2*two_p[1]*two_p[1]*two_p[1];three_p[11] = 3*two_p[1]*two_p[1]*two_p[2]
    three_p[12] = 6*two_p[1]*two_p[1]*two_p[2];three_p[13] = 3*two_p[1]*two_p[1]*two_p[2]
    three_p[14] = 6*two_p[1]*two_p[2]*two_p[2];three_p[15] = 1*two_p[2]*two_p[2]*two_p[2]
    return three_p


def threecell_motif(Adj_matrix):
    mat_nodes = np.arange(Adj_matrix.shape[0])
    pick3cell = list(itertools.combinations(mat_nodes, 3))
    patterns = np.zeros(16).astype(int)    # 16 patterns 
    
    from tqdm import tqdm
    for ii in tqdm(range(len(pick3cell))):
        pick = list(pick3cell[ii])
        submat = Adj_matrix[np.ix_(pick,pick)]
        mask = np.array([[0,1,1],[1,0,1],[1,1,0]])
        submat = submat*mask    # diagonal -> 0
        nz1, nz2 = np.nonzero(submat)
        if nz1.size==0:
            patterns[0] += 1
        elif nz1.size==1:
            patterns[1] += 1
        elif nz1.size==2:
            if (nz1[0]==nz2[1]) and (nz1[1]==nz2[0]):
                patterns[2] += 1
            elif nz1[0]==nz1[1]:
                patterns[3] += 1
            elif nz2[0]==nz2[1]:
                patterns[4] += 1
            elif (nz2[0]==nz1[1]) or (nz2[1]==nz1[0]):
                patterns[5] += 1
        elif nz1.size==3:
            if (np.unique(nz1).size==3) and (np.unique(nz2).size==2):
                patterns[6] += 1
            elif (np.unique(nz1).size==2) and (np.unique(nz2).size==3):
                patterns[7] += 1
            elif (np.unique(nz1).size==2) and (np.unique(nz2).size==2):
                patterns[9] += 1
            elif (np.unique(nz1).size==3) and (np.unique(nz2).size==3):
                patterns[10] += 1
        elif nz1.size==4:
            if np.array_equal(np.sort(nz1),np.sort(nz2)):
                patterns[8] += 1
            elif (np.unique(nz2).size==2) and (np.unique(nz1).size==3):
                patterns[11] += 1
            elif (np.unique(nz1).size==2) and (np.unique(nz2).size==3):
                patterns[13] += 1
            else:
                patterns[12] += 1
        elif nz1.size==5:
            patterns[14] += 1
        elif nz1.size==6:
            patterns[15] += 1
            
    return patterns


def fast_3motif_count(Adj_matrix):     
    g = ig.Graph.Weighted_Adjacency(Adj_matrix)
    temp = np.copy(np.array(g.triad_census()))
    order_change = np.copy(temp)
    order_change[9] = temp[8]; order_change[10] = temp[9]; order_change[8] = temp[10]
    order_change[12] = temp[13]; order_change[13] = temp[12]
    
    patterns = np.copy(order_change).astype(int)
    
    return patterns


def constrained_ER(node_number, motif_count):
    c_ER = np.zeros((node_number, node_number))
    unc = [0]*motif_count[0]; uni = [1]*motif_count[1]; bi = [2]*motif_count[2]
    assign = unc + uni + bi
    comb = list(itertools.combinations(np.arange(node_number), 2))
    
    import random
    random.shuffle(assign)
    
    for ii in range(len(comb)):
        idx = comb[ii]
        if assign[ii]==0:
            continue
        elif assign[ii]==1:
            nn = random.random()
            if nn < 0.5:
                c_ER[idx] = 1
            else:
                c_ER[idx[::-1]] = 1
        elif assign[ii]==2:
            c_ER[idx] = 1; c_ER[idx[::-1]] = 1
            
    return c_ER


def find_threecell_motif(Adj_matrix, pattern_number):
    mat_nodes = np.arange(Adj_matrix.shape[0])
    pick3cell = list(itertools.combinations(mat_nodes, 3))
    mask = np.array([[0,1,1],[1,0,1],[1,1,0]])
    motif_list = []
    weights = []
    for ii in tqdm(range(len(pick3cell))):
        pick = list(pick3cell[ii])
        submat = Adj_matrix[np.ix_(pick,pick)]
        submat = submat*mask
        nz1, nz2 = np.nonzero(submat)
        
        if nz1.size==0:
            pat = 0
        elif nz1.size==1:
            pat = 1
        elif nz1.size==2:
            if (nz1[0]==nz2[1]) and (nz1[1]==nz2[0]):
                pat = 2
            elif nz1[0]==nz1[1]:
                pat = 3
            elif nz2[0]==nz2[1]:
                pat = 4
            elif (nz2[0]==nz1[1]) or (nz2[1]==nz1[0]):
                pat = 5
        elif nz1.size==3:
            if (np.unique(nz1).size==3) and (np.unique(nz2).size==2):
                pat = 6
            elif (np.unique(nz1).size==2) and (np.unique(nz2).size==3):
                pat = 7
            elif (np.unique(nz1).size==2) and (np.unique(nz2).size==2):
                pat = 9
            elif (np.unique(nz1).size==3) and (np.unique(nz2).size==3):
                pat = 10
        elif nz1.size==4:
            if np.array_equal(np.sort(nz1),np.sort(nz2)):
                pat = 8
            elif (np.unique(nz2).size==2) and (np.unique(nz1).size==3):
                pat = 11
            elif (np.unique(nz1).size==2) and (np.unique(nz2).size==3):
                pat = 13
            else:
                pat = 12
        elif nz1.size==5:
            pat = 14
        elif nz1.size==6:
            pat = 15
        
        if pat == pattern_number:
            motif_list.append(pick)
            ww = list(submat[np.nonzero(submat)])
            weights.append(ww)
        else:
            continue

    return motif_list, weights

