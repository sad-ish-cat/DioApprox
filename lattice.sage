


def generate_bases(n): 

	candidates = [[[a,b],[c,d]] for a in range(n+1) for b in range(a) for d in range(n+1) for c in range(d)] 
	candidates = [B for B in candidates if (B[0][0]*B[1][1]-B[0][1]*B[1][0] == n)] 
	
	return candidates
	
	


def count_by_left_right(index, max_steps=100):

	count = 0
	for i in range(index):
		count += len(find_left_right_relative([[i,index],[2*i+1,2*index],[i+1,index]],max_steps = max_steps))
	return count
	
	

def test_cases(indices = []):

	for i in indices:
		if not (count_by_left_right(i) == len(generate_bases(i))):
			return False
	return True
	
	
	
def change_basis(v,b1,b2):
		
		det = b1[0]*b2[1] - b1[1]*b2[0]
		a = (v[0]*b2[1] - v[1]*b2[0])/det
		b = (b1[0]*v[1] - b1[1]*v[0])/det
		
		return [a,b]
	
def next_lattice(v1, v2, move, bases):
	
	
	if move == "L":
		v1_new = [v1[0]+v2[0],v1[1]+v2[1]]
		v2_new = v2
		
	if move == "R":
		v1_new = v1
		v2_new = [v1[0]+v2[0],v1[1]+v2[1]]
		
		
	is_w_best = False
	left_right = ""
	
	sequence = find_left_right_relative([v2_new,[],v1_new], max_steps = 100)
	
	w1 = create_data_from_left_right(sequence)[2]
	w2 = create_data_from_left_right(sequence)[0]
	
	
	return [[change_basis(v1_new, w1, w2), change_basis(v2_new, w1, w2)], sequence]
	
