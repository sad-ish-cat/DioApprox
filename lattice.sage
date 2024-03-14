# generate all lattices in nxn box that can appear as a state
def generate_lattices(n): 
	candidates = [[[a,b],[c,d]] for a in range(n) for b in range(n) for c in range(n) for d in range(n)] #makes all possible bases in nxn box
	candidates = [B for B in candidates if (B[0][0]*B[1][1]-B[0][1]*B[1][0] == n)] #filters out lattices with index not equal to n
	#candidates = [B in candidates if]
	
	return candidates