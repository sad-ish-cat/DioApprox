# generate all lattices in nxn box that can appear as a state
def generate_basises(n): 
	candidates = [[[a,b],[c,d]] for a in range(n+1) for b in range(n+1) for c in range(n+1) for d in range(n+1)] #makes all possible bases in nxn box
	candidates = [B for B in candidates if (B[0][0]*B[1][1]-B[0][1]*B[1][0] == n)] #filters out lattices with index not equal to n
	candidates = [B for B in candidates if (B[0][0]> B[0][1] and B[1][0] < B[1][1])]
	
	return candidates