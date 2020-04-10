import numpy as np
import timeit
from scipy.spatial.distance import cdist
 
# define a dot product function used for the rotate operation
def v_dot(a):return lambda b: np.dot(a,b)
 
class lattice_SAW:
    def __init__(self,N,l0):
        self.N = N
        self.l0 = l0
        # initial configuration. Usually we just use a straight chain as inital configuration
        self.init_state = np.dstack((np.arange(N),np.zeros(N),np.zeros(N)))[0]
        self.state = self.init_state.copy()
 
        # define a rotation matrix
        # 9 possible rotations: 3 axes * 3 possible rotate angles(90,180,270)
        self.rotate_matrix = np.array([[[1,0,0],[0,0,-1],[0,1,0]],[[1,0,0],[0,-1,0],[0,0,-1]]
        ,[[1,0,0],[0,0,1],[0,-1,0]],[[0,0,1],[0,1,0],[-1,0,0]]
        ,[[-1,0,0],[0,1,0],[0,0,-1]],[[0,0,-1],[0,1,0],[-1,0,0]]
        ,[[0,-1,0],[1,0,0],[0,0,1]],[[-1,0,0],[0,-1,0],[0,0,1]]
        ,[[0,1,0],[-1,0,0],[0,0,1]]])
 
    # define pivot algorithm process where t is the number of successful steps
    def walk(self,t):
        acpt = 0
        # while loop until the number of successful step up to t
        while acpt <= t:
            pick_pivot = np.random.randint(1,self.N-1) # pick a pivot site
            pick_side = np.random.choice([-1,1]) # pick a side
 
            if pick_side == 1:
                old_chain = self.state[0:pick_pivot+1]
                temp_chain = self.state[pick_pivot+1:]
            else:
                old_chain = self.state[pick_pivot:]
                temp_chain = self.state[0:pick_pivot]
 
            # pick a symmetry operator
            symtry_oprtr = self.rotate_matrix[np.random.randint(len(self.rotate_matrix))]
            # new chain after symmetry operator
            new_chain = np.apply_along_axis(v_dot(symtry_oprtr),1,temp_chain - self.state[pick_pivot]) + self.state[pick_pivot]
 
            # use cdist function of scipy package to calculate the pair-pair distance between old_chain and new_chain
            overlap = cdist(new_chain,old_chain)
            overlap = overlap.flatten()
 
            # determinte whether the new state is accepted or rejected
            if len(np.nonzero(overlap)[0]) != len(overlap):
                continue
            else:
                if pick_side == 1:
                    self.state = np.concatenate((old_chain,new_chain),axis=0)
                elif pick_side == -1:
                    self.state = np.concatenate((new_chain,old_chain),axis=0)
                acpt += 1
 
        # place the center of mass of the chain on the origin
        self.state = self.l0*(self.state - np.int_(np.mean(self.state,axis=0)))
 
N = 100 # number of monomers(number of steps)
l0 = 1 # bond length(step length)
t = 1000 # number of pivot steps
 
chain = lattice_SAW(N,l0)

 
#timeit chain.walk(t)
#1
#1 loops, best of 3: 2.61