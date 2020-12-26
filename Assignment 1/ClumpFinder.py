def ClumpFinder(k,L,t,Genome):
    #Length of Genome
    N = len(Genome)
    
    #Storing Frequent patterns
    freq_patterns = []
    
    for i in range(N-L+1):
        #choosing region of length L in the Genome
        region = Genome[i:i+L]
        #Calculating the Frequency array in the first iteration
        if i==0:
            freq_dict = {}
            for j in range(L-k+1):
                #for each kmer in this region
                kmer = region[j:j+k]
                #Update count
                freq_dict[kmer] = freq_dict.get(kmer,0) + 1
            #Append to freq_patterns if the kmer occurs at least t times
            for pattern in freq_dict:
                if freq_dict[pattern] >= t:
                    freq_patterns.append(pattern)
        else:
            #Reduce count of the first kmer of previous region by 1
            first = Genome[i-1:k+i-1]
            freq_dict[first] -= 1
            #Increase count of the last kmer of current region by 1
            last = Genome[i+L-k:i+L]
            freq_dict[last] = freq_dict.get(last,0) + 1
            
            #If the last kmer's count is at least t and not present already in freq_patterns, append it
            if freq_dict[last] >=t and last not in freq_patterns:
                freq_patterns.append(last)
                
    #Sort patterns in lexicological order
    freq_patterns = sorted(freq_patterns)
    return freq_patterns

Genome = input()
inputs = input().split()
k = int(inputs[0])
L = int(inputs[1])
t = int(inputs[2])

start_time = time.time()
print(" ".join(ClumpFinder(k,L,t,Genome)))
print(time.time()-start_time)