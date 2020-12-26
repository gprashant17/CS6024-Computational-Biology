def FittingAlignment(v,w,n1,n2):
    """
    Computes the fitting alignment score of two strings and returns the alignment
    v: Larger String 
    w: Smaller String
    n1: Length of v
    n2: Length of w
    """
    path = {} #To store the annotation corresponding to each entry in the DAG matrix, as to whether it was a            mismatch/indel
    
    DAG = {i:[0]*(n1+1) for i in range(n2+1)} #DAG Matrix -> stored as dictionary
    
    for i in range(1,n2+1):
        DAG[i][0] = -i
        for j in range(1,n1+1):
            #if a match occurs
            if w[i-1] == v[j-1]:
                opt = [DAG[i-1][j-1]+1,DAG[i-1][j]-1,DAG[i][j-1]-1]
            else:
                opt = [DAG[i-1][j-1]-1,DAG[i-1][j]-1,DAG[i][j-1]-1]
            DAG[i][j] = max(opt)
            #As we are backtracking, the order of preferences are reversed (i.e., d>i>m)
            opt.reverse()
            k = opt.index(max(opt))
            if k==0:
                path[(i,j)] = "d"
            elif k==1:
                path[(i,j)] = "i"
            elif k==2:
                path[(i,j)] = "m"
                
    
    #Obtain the last row of the DAG matrix
    last_row = DAG[n2]
    Score = max(last_row)
    last_row.reverse()
    ind =  n1 - last_row.index(Score)
    pos = (n2,ind)
    str1 = ""
    str2 = ""
    
    #Backtracking... stop when first row is reached
    while pos[0]!=0:
        if path[pos] == "m":
            str1 = v[pos[1]-1] + str1
            str2 = w[pos[0]-1] + str2
            pos = (pos[0]-1,pos[1]-1)
        elif path[pos] == "i":
            str1 = "-" + str1
            str2 = w[pos[0]-1] + str2
            pos = (pos[0]-1,pos[1])
        elif path[pos] == "d":
            str1 = v[pos[1]-1] + str1
            str2 = "-" + str2
            pos = (pos[0],pos[1]-1)
            
    return Score, str1,str2           

v = input()
w = input()
n1 = len(v)
n2 = len(w)
Score,str1,str2 = FittingAlignment(v,w,n1,n2)
print(Score)
print(str1)
print(str2)