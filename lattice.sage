load("sternbrocot.sage")
from sage.rings.continued_fraction import last_two_convergents

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
	
# Gives a conjectured formula for the number of good bases of index n.
def count_by_left_right(index, max_steps=100):
	count = 0
	for i in range(index):
		count += len(find_left_right_relative(((i,index),(),(i+1,index)),max_steps = max_steps))
	return count
		
		
		
def test_cases(indices = []):
	
	for i in indices:
		if not (count_by_left_right(i,max_steps=(max(indices)+10)) == len(generate_bases(i))):
			print(i)
			return False
	return True
			
			
			
def change_basis(v,b1,b2):
	
	det = b1[0]*b2[1] - b1[1]*b2[0]
	a = (v[0]*b2[1] - v[1]*b2[0])/det
	b = (b1[0]*v[1] - b1[1]*v[0])/det
	
	return (a,b)
	
@cached_function
def next_basis(old, move, max_steps = 100):
	
	v1, v2 = old;
	
	
	
	if move == "L":
		v1_new = (v1[0]+v2[0],v1[1]+v2[1])
		v2_new = v2
		
	if move == "R":
		v1_new = v1
		v2_new = (v1[0]+v2[0],v1[1]+v2[1])
		
	sequence = find_left_right_relative((v2_new,(),v1_new), max_steps = max_steps)
	
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
		{basis : {next_basis(basis, "L", max_steps = n+3)[0]:["L/" + next_basis(basis, "L", max_steps = n+3)[1]],
			 next_basis(basis, "R", max_steps = n+3)[0]:["R/" + next_basis(basis, "R", max_steps = n+3)[1]]}
			for basis in generate_bases(n)},
		format='dict_of_dicts', loops=True, multiedges=True,
		immutable=False)
		
def path_to_lr(path):

	if len(path) == 0:
		return ""
	
	moves = ""
	sequence = ""
	
	
	
	for i in range(0, len(path) - 1):
		next_basis_l = next_basis(path[i], "L")
		if next_basis_l[0] == path[i+1]:
			moves += "L"
			sequence += next_basis_l[1]
		else:
			moves += "R"
			sequence += next_basis(path[i], "R")[1]
		
	return (moves, sequence)

def list_cycles_lr(n, max_length = 5):

	g = get_graph(n);
	g.remove_loops();
	return [path_to_lr(c) for c in g.all_cycles_iterator(max_length = max_length)]
	
# Convert LR sequence to continued fraction
def to_continued_fraction(lrseq, cyclic=False):

	if lrseq == "":
		return []
		
	prev = ""
	count = 0
	continued_frac = []
	
	for ch in lrseq:
		if ch == prev:
			count += 1;
		else:
			prev = ch;
			if count >= 1:
				continued_frac.append(count)
			count = 1
	continued_frac.append(count)
	
	if not(cyclic):
		return continued_frac
	if len(continued_frac) == 1:
		return [Infinity]
	if len(continued_frac) % 2 == 0:
			return continued_frac
	return continued_frac[1:-1] + [continued_frac[0] + continued_frac[-1]]

# Convert eventually periodic LR sequence to continued fraction
def to_eventually_periodic_cf(transient, period):
	cf1 = to_continued_fraction(transient)
	if len(transient) == 0 or transient[0] == "L":
		cf1 = [0] + cf1
	cf2 = to_continued_fraction(period, cyclic=True)
	
	last_transient = transient[-1] if len(transient) > 0 else "R"
	
	if last_transient != period[0] and period[-1] != period[0]:
		return (cf1, cf2)
	elif last_transient == period[0] and period[-1] != period[0]:
		return (cf1[:-1] + [cf1[-1] + cf2[0]], cf2[1:] + [cf2[0]])
	elif last_transient != period[0] and period[-1] == period[0]:
		return (cf1 + [to_continued_fraction(period, cyclic=False)[0]], cf2)
	elif last_transient == period[0] and period[-1] == period[0]:
		return (cf1[:-1] + [cf1[-1] + to_continued_fraction(period, cyclic=False)[0]], cf2)

def test_to_eventually_periodic_cf():
	assert to_eventually_periodic_cf("RR", "LLLRLRR") == ([2], [3, 1, 1, 2])
	assert to_eventually_periodic_cf("RR", "RRRLRLL") == ([5], [1, 1, 2, 3])
	assert to_eventually_periodic_cf("RR", "LLLRLL")  == ([2, 3], [1, 5])
	assert to_eventually_periodic_cf("RR", "RRLLRRR") == ([4], [2, 5])
	assert to_eventually_periodic_cf("LL", "LLLRLRR") == ([0, 5], [1, 1, 2, 3])
	assert to_eventually_periodic_cf("LL", "RRRLRLL") == ([0, 2], [3, 1, 1, 2])
	assert to_eventually_periodic_cf("LL", "LLLRLLL") == ([0, 5], [1, 6])
	assert to_eventually_periodic_cf("LL", "RRLLRRR") == ([0, 2, 2], [2, 5])
	assert to_eventually_periodic_cf("", "RRLR") == ([2], [1, 3])
	assert to_eventually_periodic_cf("", "RRLLL") in [([2], [3, 2]), ([], [2, 3])]
	assert to_eventually_periodic_cf("", "LLRRR") == ([0], [2, 3])
	assert to_eventually_periodic_cf("", "LLRL") == ([0, 2], [1, 3])

# Convert continued fraction to LR sequence
@cached_function
def to_lr_seq(ctdfrac):
	
	lr_seq = ""
	
	for i in range(len(ctdfrac)):
		if i % 2 == 0:
			lr_seq += "R" * ctdfrac[i]
		else:
			lr_seq += "L" * ctdfrac[i]
			
	return lr_seq[:len(lr_seq)-1]
	
def approx_value(ctdfrac, fast=False):
	
	max_term = max(ctdfrac)
	
	if max_term == Infinity:
		return Infinity
	
	rotations = []
	for i in range(len(ctdfrac)):
		if ctdfrac[i] >= max_term - 1:
			rotations.append(ctdfrac[i:] + ctdfrac[:i])
		
	largest = 0
	
	for rot in rotations:
		if fast:
			approx = fast_value(continued_fraction([[rot[0]],rot[1:] + [rot[0]]])) +\
				fast_value(continued_fraction([[0], rot[::-1]]));
		else:
			approx = continued_fraction([[rot[0]],rot[1:] + [rot[0]]]).value() +\
				continued_fraction([[0], rot[::-1]]).value();
		if approx > largest:
			largest = approx
			
	return largest
	
	
def list_cycle_approx(n, max_length = 5, sort = True):

	cycles = list_cycles_lr(n, max_length);
	ret = [
		[N(max(approx_value(to_continued_fraction(string, cyclic = True))
			 for string in cycle))] + list(cycle)
			for cycle in cycles
	]
	if sort:
		ret.sort()
	return ret
	
# def min_lambda(ctdfrac):
#		n = len(ctdfrac)
#		
#		if n == 0:
#			return 0
#		
#		max_term = max(ctdfrac)
#		
#		if max_term == Infinity:
#			return Infinity
#		
#		largest = 0
#		for i in range(len(ctdfrac)):
#			if ctdfrac[i] >= max_term - 1:
#				print(ctdfrac[i:n-((n-i-1)%2)], [0] + ctdfrac[i-1:-((i-1)%2):-1])
#				approx = continued_fraction(ctdfrac[i:n-((n-i-1)%2)]).value() +\
#				continued_fraction([0] + ctdfrac[i-1:-((i-1)%2):-1]).value();
#				if approx > largest:
#					largest = approx
#				
#		return largest

@cached_function
def min_lambda_new(ctdfrac):

	ctdfrac = tuple(ctdfrac)
	lr_seq = to_lr_seq(ctdfrac)
	
	n = len(ctdfrac)
	
	if n == 0:
		return 0
	
	max_term = max(ctdfrac)
	
	if max_term == Infinity:
		return Infinity
	
	largest = 0
	for i in range(len(ctdfrac)):
		if ctdfrac[i] >= max_term - 1: 
			approx = (create_number_from_left_right(to_lr_seq(tuple(
				ctdfrac[i:n-((n-i-1)%2)]))) +
				create_number_from_left_right(to_lr_seq(tuple(
					(0,) + ctdfrac[-n+i-1:-n-((i-1)%2):-1])))
				)
			if approx > largest:
				largest = approx
			
	return largest
	

	
# def find_good_paths(n, lambda_limit, max_length = 5):
#		
#		graph = get_graph(n)
#		prev = set([(node, ) for node in graph.vertices()]) # k = 0
#		current = set()
#		for k in [1..max_length]:
#			for oldpath in prev:
#				vtx = oldpath[-1];
#				for edge in graph.outgoing_edge_iterator(vtx,labels=False):
#					newpath = oldpath + (edge[1],);
#					
#					# Check head segment
#					if newpath[1:] in prev:
#						xi_moves, nxi_moves = path_to_lr(newpath)
#						if (min_lambda_new(tuple(to_continued_fraction(xi_moves))) <= lambda_limit
#							and min_lambda_new(tuple(to_continued_fraction(nxi_moves))) <= lambda_limit):
#								current.add(newpath)
#			prev = current
#			print(len(prev), "good paths of length", k)
#			current = set()
#		return prev


# Adapted from the usual value() method for CF's, but returns a numeric value instead.
# Cf. /usr/lib/python3/dist-packages/sage/rings/continued_fraction.py
def fast_value(cf): 
	if cf._x1 and cf._x1[0] < 0:
		return -fast_value(-cf)

	if cf._x2[0] is Infinity:
			return cf._rational_()
	
	# determine the equation for the purely periodic cont. frac. determined
	# by cf._x2
	p0, q0, p1, q1 = last_two_convergents(cf._x2)
	
	# now x is one of the root of the equation
	#		q1 x^2 + (q0 - p1) x - p0 = 0
	D = (q0-p1)**2 + 4*q1*p0
	x = ((p1 - q0) + sqrt(RR(D))) / (2*q1)
	
	# we add the preperiod
	p0, q0, p1, q1 = last_two_convergents(cf._x1)
	return (p1*x + p0) / (q1*x + q0)
		 
