import sys
def NeighborJoining(D,N):
    """
    Function that Recursively Constructs the evolutionary tree using Neighbor-Joining Algorithm
    Input: 
        D: Distance Matrix (Dictionary of dictionaries)
        N: Number of Species
    Returns:
        T: Evolutionary Tree (Dictionary of dictionaries)
    """
    ##Base case, when there are only 2 species, then connect them directly
    if  N==2:
        #Tree is stored as a dictionary of dictionary
        T = {}
        n = max(D.keys())
        m = min(D.keys())
        nodes = list(D.keys())
        T[m] = {n:round(D[nodes[0]][nodes[1]],3)}
        T[n] = {m:round(D[nodes[0]][nodes[1]],3)}
    
        return T
    else:
        T = []
        N = len(D)
        #Constructing Dstar and finding the indices (min_pos(i,j)) of the minimum value
        Dstar,(i,j) = DstarConstruction(D,N)
        #Distance between i and j
        delta_ij = (sum(D[i].values())-sum(D[j].values()))/(N-2)    
        
        #Limb length each cherry
        limb_i = round((D[i][j]+delta_ij)/2,3)
        limb_j = round((D[i][j]-delta_ij)/2,3)
        
        #Modify distance matrix by removing columns i and j and replacing with another column
        D = modifyD(D,N,i,j)     
        #Recursively calling the function with modified distance matrix
        T = NeighborJoining(D,N-1)
        
        ##Adding the cherries to the existing tree
        n = max(D.keys())
        
        if i in T:
            T[i][n] = limb_i
        else:
            T[i] = {n:limb_i}
            
        if j in T:
            T[j][n] = limb_j
        else:
            T[j] = {n:limb_j}
            
        T[n][i] = limb_i
        T[n][j] = limb_j
        
    return T

def DstarConstruction(D,N):
    """
    Function that constructs Dstar matrix and finds the position corresponding to minimum element of 
    Dstar
    """
    Dstar = {}
    MIN = sys.maxsize

    for key1 in sorted(D.keys()):
        Dstar[key1] = {key1:0}
        for key2 in sorted(D.keys()):
            if key1!=key2:       #As we do not change the diagonal elements         
                #Calculating entries of Dstar
                Dstar[key1][key2] = (N-2)*D[key1][key2] - sum(D[key1].values()) - sum(D[key2].values())
                #The command below ensures only the first ocurring minimum value is taken 
                #into consideration
                if Dstar[key1][key2]<MIN:
                    MIN = Dstar[key1][key2]
                    #Storing the minimum index
                    min_pos = (key1,key2)
    return Dstar,min_pos

def modifyD(D,N,i,j):
    """
    Function to remove entries corresponding to species i and j, and replacing with a new species 
    that is equivalent to both the cherries combined.
    """
    #As the new species should have index starting from n
    n = max(D.keys())+1
    
    D_mod = {n:{n:0}}
    dij = D[i][j] #ij distance
    
    for key in D.keys():
        if key!=i and key!=j:            
            D_mod[key] = D[key].copy()
            ##Deleting i and jth column
            del D_mod[key][i]
            del D_mod[key][j]
            
            #Distance between other species and the new entry
            D_mod[n][key] = ((D[key][i]+D[key][j]-dij)/2)        
            D_mod[key][n] = D_mod[n][key]
        
    return D_mod

#Getting the number of species
N = int(input())
#Getting the distance matrix, stored as dictionary of dictionary
D = {} 
for i in range(N):
    d_vec = list(map(lambda x: int(x),input().split("\t")[:-1]))
    D[i] = {}
    for j in range(N):
        D[i][j]  = d_vec[j]
Tree = NeighborJoining(D,N)
#Printing Output
for node in sorted(list(Tree.keys())):    
    neigh = sorted(list(Tree[node].keys()))
    for near in neigh:
        print("%d"%node+"->"+"%d"%near+":%.3f"%(Tree[node][near]))