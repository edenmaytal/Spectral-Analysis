import mykmeanssp as myk
import numpy as np
import pandas as pd

def kmeans_pp(k, input, max_iter=300):
    n = len(input) #num of vecs in input 
    vecs = input
    np.random.seed(0)
    ind=int(np.random.choice(n))#sample a random index
    miu_1= vecs[ind] #get miu1 acoording to the index we sampled
    res_ind= [ind] #list of chosen miu indicies
    Z=1
    D= [0 for i in range(n)]
    miu= [miu_1] #list of chosen miu vectors
    while Z<k: 
        for i,vec in enumerate(vecs): #go over xi (vecs)
            if vec!=0:
                min_dist= np.linalg.norm(np.array(miu_1)-np.array(vec))**2 
                for j in range(1,Z): #go over the centroids we found so far
                    dist= np.linalg.norm(np.array(vec)-np.array(miu[j]))**2
                    if dist < min_dist:
                        min_dist = dist
                D[i]= min_dist
        Z+=1
        S= sum(D)
        prob=[(D[i]/S) for i in range(n)] #list of probability for each xi
        ind=int(np.random.choice(a=n, p=prob))
        miu.append(vecs[ind])
        res_ind.append(ind)
    
    inds_to_print=""
    for i in range(len(res_ind)-1):
        inds_to_print += str(res_ind[i]) + ","
    print(inds_to_print, end = "")
    print(res_ind[-1])

    vecs_new = [0 for i in range(k)] # change the output list vecs - k initial centroids first
    for i in range(k):
        vecs_new[i] = vecs[res_ind[i]]
        vecs[res_ind[i]] = None
    for i in range(n):
        if vecs[i] != None:
            vecs_new.append(vecs[i])
    return vecs_new,res_ind


if __name__ == '__main__':
    import sys
    args= sys.argv[1:]
    k = int(args[0])
    goal = args[1]
    file_name = args[2]
    input = pd.read_csv(file_name, header = None)
    vecs_input = input.values.tolist()
    # fit prints the result from the C code
    if goal == "wam":
        myk.fit(0, vecs_input, len(vecs_input[0]))
    elif goal == "ddg":
        myk.fit(1, vecs_input, len(vecs_input[0])) 
    elif goal == "lnorm":
        myk.fit(2, vecs_input, len(vecs_input[0]))
    elif goal == "jacobi":
        myk.fit(3, vecs_input, len(vecs_input[0]))   
    else: # goal = 'spk'
        vecs = myk.fit_spk(k, vecs_input, len(vecs_input[0])) # needs to return updated k
        k = len(vecs[0])
        vecs, result_ind = kmeans_pp(k, vecs) 
        myk.fit_kmeans(vecs, k ,len(vecs[0]))


    
