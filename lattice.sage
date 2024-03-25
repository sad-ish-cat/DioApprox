load("sternbrocot.sage")

@cached_function
def generate_bases(n): 
	candidates = []
	for a in range(1, n + 1):
		for d in range(1, n + 1): 
			if a * d > n: 
				for b in divisors(a * d - n): 
					c = (a * d - n) // b
					if a > b and c < d: 
						candidates.append(((a, b), (c, d)))
			else:
				if a * d == n: 
					# first let b = 0 
					candidates.extend([((a, 0), (c, d)) for c in range(0, d)])
					# then let c = 0
					candidates.extend([((a, b), (0, d)) for b in range(1, a)])
	return [basis for basis in candidates if gcd(flatten(basis)) == 1]
	
	
@cached_function 
def generate_bases_old(n): 
	candidates = [((a,b),(c,d)) for a in range(n+1) for b in range(a) for d in range(n+1) for c in range(d)] 
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
	
	return (a,b)
	
def next_basis(old, move, max_steps = 100):
	
	v1, v2 = old;
	
	if move == "L":
		v1_new = [v1[0]+v2[0],v1[1]+v2[1]]
		v2_new = v2
		
	if move == "R":
		v1_new = v1
		v2_new = [v1[0]+v2[0],v1[1]+v2[1]]
		
		
	sequence = find_left_right_relative([v2_new,[],v1_new], max_steps = max_steps)
	
	w1 = create_data_from_left_right(sequence)[2]
	w2 = create_data_from_left_right(sequence)[0]
	
	
	return [(change_basis(v1_new, w1, w2), change_basis(v2_new, w1, w2)), sequence]

	
@cached_function
def state_machine(n):
	bases = generate_bases(n);
	return [[old, move] + next_basis(old, move, n+3) for old in bases for move in ["L", "R"]]
	
@cached_function
def get_graph(n):
	return DiGraph(
		{basis : [next_basis(basis, "L", max_steps = n+3)[0],
		   next_basis(basis, "R", max_steps = n+3)[0]]
			for basis in generate_bases(n)},
		format='dict_of_lists', loops=True, multiedges=True,
		immutable=False)
		
		
def cycle_to_lr(cycle):

	if len(cycle) <= 1:
		return ""
	
	moves = ""
	sequence = ""
	
	for i in range(0, len(cycle) - 1):
		if next_basis(cycle[i], "L")[0] == cycle[i+1]:
			move = "L"
		else:
			move = "R"
		moves += move
		sequence += next_basis(cycle[i],move)[1]
	return (moves, sequence)

def list_cycles_lr(n, max_length = 5):

	g = get_graph(n);
	g.remove_loops();
	return [cycle_to_lr(c) for c in g.all_cycles_iterator(max_length = max_length)]
	
	
def count_lr_consecutive(lrseq):
	
	seq = list(lrseq)
	if (len(seq) == 1):
		return 1
	elif (seq[0] == seq[1]):
		return 1 + count_lr_consecutive(lrseq[1:])
	else:
		return 1
	
	
def to_continued_fraction(lrseq):

	if lrseq == "":
		return []
		
	seq = lrseq
	continued_frac = []
	
	while(len(seq)>0):
	
		num = count_lr_consecutive(seq)
		continued_frac.append(num)
		
		seq = seq[num:]
	
	if len(continued_frac) == 1:
		return [Infinity]
	if len(continued_frac) % 2 == 0:
	    return continued_frac
	return continued_frac[1:-1] + [continued_frac[0] + continued_frac[-1]]
	
	
	
def approx_value(ctdfrac):
	
	biggest = max(ctdfrac)
	
	if biggest == Infinity:
		return Infinity
	
	rotations = []
	for i in range(len(ctdfrac)):
		if ctdfrac[i] >= biggest - 1:
			rotations.append(ctdfrac[i:] + ctdfrac[:i])
		
	largest = 0
	
	for rot in rotations:
		approx = continued_fraction([[rot[0]],rot[1:] + [rot[0]]]).value() +\
			continued_fraction([[0], rot[::-1]]).value();
		if approx > largest:
			largest = approx
			
	return largest
	
	
def list_cycle_approx(n, max_length = 5, sort = True):

	cycles = list_cycles_lr(n, max_length);
	ret = [
		[N(max(approx_value(to_continued_fraction(string)) for string in cycle))] + list(cycle)
	    for cycle in cycles
	]
	if sort:
		ret.sort()
	return ret
