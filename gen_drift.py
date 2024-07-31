import math

def wf(pc):
    tp_dim = pc + 1
    trans_ps = []
    for i in range(tp_dim):
        trans_ps_row = []
        for j in range(pc+1): trans_ps_row += [0]
        trans_ps += [trans_ps_row]
    for i in range(tp_dim):
        for j in range(tp_dim):
            trans_ps[i][j] = math.comb(pc, j)*((i/pc)**j)*(1 - i/pc)**(pc-j)
    return trans_ps