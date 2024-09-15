load("lattice.sage")

# raney.sage: Compute transducers by the methods of Raney

matrix_L = matrix([[1,0], [1,1]])
matrix_R = matrix([[1,1], [0,1]])
matrix_L_inv = matrix_L^-1
matrix_R_inv = matrix_R^-1

def PROD(word):
	ret = identity_matrix(2);
	for ch in word:
		if ch == "L":
			ret *= matrix_L
		elif ch == "R":
			ret *= matrix_R
		else:
			raise ValueError("Invalid character " + repr(ch) + " in input " +
				repr(word))
	return ret
	
def to_CB(mat_orig):
	mat = mat_orig;
	word = "";
	while True:
		if mat[1,0] >= mat[1,1]:
			word = "L" + word;
			mat = mat * matrix_L_inv;
		elif mat[0,1] >= mat[0,0]:
			word = "R" + word;
			mat = mat * matrix_R_inv;
		else:
			break;
	assert mat * PROD(word) == mat_orig;
	return mat, word;

def to_RB(mat_orig):
	mat = mat_orig;
	word = "";
	while True:
		if mat[0,1] >= mat[1,1]:
			word = word + "R";
			mat = matrix_R_inv * mat;
		elif mat[1,0] >= mat[0,0]:
			word = word + "L";
			mat = matrix_L_inv * mat;
		else:
			break;
	assert PROD(word) * mat == mat_orig;
	return mat, word;

# Get the s'th DB matrix, 0 <= s < n, n prime
# (for composite n this generates some of the DB matrices)
def get_DB_prime(n, s):
	mat_orig = matrix([[n - s, s], [n - s - 1, s + 1]]);
	mat, word = to_CB(mat_orig);
	# assert get_index(mat, n) == s;
	return mat, word;

def get_index(mat, n):
	if n == 1:
		return 0;
	try:
		return (GF(n)(mat[0,0]) / GF(n)(mat[1,0] - mat[0,0])).lift();
	except ZeroDivisionError:
		return (GF(n)(mat[0,1]) / GF(n)(mat[1,1] - mat[0,1])).lift();
	
@cached_function
def get_DB_prime_all(n):
	# assert is_prime(n)
	return [get_DB_prime(n, s) for s in range(n)]

def immed_offshoots(word):
	ret = [];
	for i in range(len(word)):
		if word[i] == "L":
			ret.append(word[0:i] + "R");
		elif word[i] == "R":
			ret.append(word[0:i] + "L");
		else:
			raise ValueError("Invalid character " + repr(word[i]) + " in input " +
				repr(word))
	ret.append(word + "L");
	ret.append(word + "R");
	return ret
	
@cached_function
def get_edges(n):
	nodes = get_DB_prime_all(n);
	edges = [];
	for s in range(n):
		M1, word1 = nodes[s];
		for V1 in immed_offshoots(word1):
			M = M1 * PROD(V1);
			M2, V2 = to_RB(M);
			assert M == PROD(V2) * M2;
			edges.append((s, get_index(M2, n), (V1, V2)))
	return edges
	
# Overrides corr. method in lattice.sage
@cached_function
def get_graph(n):
	return DiGraph(get_edges(n), loops=True, multiedges=True);

# Overrides corr. method in lattice.sage
def find_good_paths(n, lambda_limit, max_length=30):
	
	graph = get_graph(n)
	
	# k = 1: Check edges
	prev = set()
	for edge in graph.edge_iterator(): # Edge is triple (v0, v1, (w_xi, w_nxi))
		xi_moves, nxi_moves = edge[2];
		if (min_lambda_new(tuple(to_continued_fraction(xi_moves))) <= lambda_limit
			and min_lambda_new(tuple(to_continued_fraction(nxi_moves))) <= lambda_limit
		):
			prev.add((edge, ))
	
	print(len(prev), "good paths of length", 1)
	current = set()
	for k in [2..max_length]:
		for oldpath in prev:
			vtx = oldpath[-1][1];
			for edge in graph.outgoing_edge_iterator(vtx):
				newpath = oldpath + (edge,);
				
				# Check head segment
				if newpath[1:] in prev:
					# Check overall approximability
					xi_moves = "".join(edge[2][0] for edge in newpath);
					nxi_moves = "".join(edge[2][1] for edge in newpath);
					if (min_lambda_new(tuple(to_continued_fraction(xi_moves)))
						<= lambda_limit
						and min_lambda_new(tuple(to_continued_fraction(nxi_moves)))
						<= lambda_limit):
							# print(xi_moves, nxi_moves);
							current.add(newpath)
		prev = current
		print(len(prev), "good paths of length", k)
		current = set()
	return prev

# Find the lowest point of the Lagrange spectrum Ln.
# Inputs:
# - n: a prime number (n = 1 is also permitted)
# - lambda_limit: a predetermined upper limit for the approximability values to
#			check (default Infinity)
# - max_length: the lengths of the paths in the graph to check (default 1000,
#			to keep looking until it stabilizes).
# - verbosity: a level of verbosity from 0 to 4:
#		= 0: print nothing
#		= 1: print final findings
#		= 2: print timely updates (number of paths of each length, running record of
#				 lowest point of Ln found)
#		= 3: print info about all cycles
#		= 4: print info about all paths
# Returns (lambda, alpha), where lambda is the lowest point of Ln, and alpha >
# lambda is a lower bound for all other points. In particular, if this method
# terminates and returns a value, it verifies the following conjectures for that
# value of n:
# - min(Ln) is an algebraic number of degree <= 2.
# - min(Ln) = min(Mn).
# - min(Ln) is an isolated point in Ln and Mn.


def find_min_Ln(n, lambda_limit = 5, max_length=1000, verbosity=1):
	
	if verbosity >= 1:
		print("p =", n)
	
	graph = get_graph(n)
	if verbosity >= 2:
		print("Graph found")
	cycle_found = False
	
	# k = 1: Start with all edges within the limit
	current = set()
	for edge in get_edges(n):
		if verbosity >= 4:
			print("Path", edge[0:2], end=" ")
		lxi = 0; lnxi = 0;
		xi_moves = "".join(edge[2][0]);
		nxi_moves = "".join(edge[2][1]);
		if (
			(lxi := max(cf := to_continued_fraction(xi_moves))) <= lambda_limit and
			(lnxi := max(ncf := to_continued_fraction(nxi_moves))) <= lambda_limit and
			(lxi := min_lambda_new(tuple(cf))) <= lambda_limit
			and ((lnxi := min_lambda_new(tuple(ncf))) <= lambda_limit)
		):
			current.add((edge,));
			if verbosity >= 4:
				print("Good", max(lxi, lnxi))
		elif verbosity >= 4:
			print("Bad", max(lxi, lnxi))
	
	if verbosity >= 2:
		print(len(current), "good paths of length", 1)
	path_counts = [len(current)]
	prev = current
	current = set()
	excluded = Infinity;
	
	
	for k in [2..max_length]:
		for oldpath in prev:
			vtx = oldpath[-1][1];
			for edge in graph.outgoing_edge_iterator(vtx):
				newpath = oldpath + (edge,);
				
				# Check head segment
				if newpath[1:] in prev:
					if verbosity >= 4:
						print("Path", [edge[0] for edge in newpath] + [newpath[-1][1]], end=" ")
					# Check overall approximability
					lxi = 0; lnxi = 0;
					xi_moves = "".join(edge[2][0] for edge in newpath);
					nxi_moves = "".join(edge[2][1] for edge in newpath);
					if (
						(lxi := max(cf := to_continued_fraction(xi_moves))) <= lambda_limit and
						(lnxi := max(ncf := to_continued_fraction(nxi_moves))) <= lambda_limit and
						(lxi := min_lambda_new(tuple(cf))) <= lambda_limit
						and ((lnxi := min_lambda_new(tuple(ncf))) <= lambda_limit)
					):
						current.add(newpath);
						if verbosity >= 4:
							print("Good", max(lxi, lnxi))
						
						# If cycle, decrease the lambda-limit.
						if newpath[0][0] == newpath[-1][1]:
							
							if verbosity >= 3:
								print("Cycle", xi_moves);
							# First compute approximately to save time for long cycles.
							lambda_new_N = max(
								approx_value(to_continued_fraction(moves, cyclic=True),
									fast = True
								) for moves in (xi_moves, nxi_moves)
							);
							if lambda_new_N < N(lambda_limit):
								# Now compute exactly.
								lambda_new = max(
									approx_value(
										to_continued_fraction(moves, cyclic=True)
									) for moves in (xi_moves, nxi_moves)
								);
								if lambda_new != lambda_limit:
									excluded = min(excluded, lambda_limit)
									lambda_limit = lambda_new
									cycle_found = True;
									if verbosity >= 2:
										print("Found new best cycle of type", [
											to_continued_fraction(xi_moves, cyclic=True),
											to_continued_fraction(nxi_moves, cyclic=True)],
											"with approximability", lambda_limit, "=", N(lambda_limit)
										);
					else:
						if verbosity >= 4:
							print("Bad", max(lxi, lnxi));
						excluded = min(excluded, max(lxi, lnxi))
		if verbosity >= 2:
			print(len(current), "good paths of length", k)
		if len(current) == 0:
			return None
		path_counts.append([len(current)])
		
		# Check for stability
		if cycle_found:
			continuations = {}
			stable = True
			tail = floor((k - sqrt(k))/2) # How many leading and trailing moves to clip.
			for path in current:
				part = path[tail:-tail-1]
				cnt1 = path[tail+1:k-tail];
				cnt2 = continuations.get(part)
				if cnt2 is None:
					continuations[part] = cnt1
				elif cnt2 != cnt1:
					stable = False
					if verbosity >= 3:
						print("Unstable")
					break
			if stable:
				if verbosity >= 3:
					print("Stable")
				# Find all cycles. Make sure no long cycles were missed.
				new_graph = DiGraph(continuations.items(), format='list_of_edges')
				cycles = [[path[0] for path in cycle]
					for cycle in new_graph.all_simple_cycles()
				]
				cycle_qualities = {}
				for cycle in cycles:
					xi_moves	= "".join(edge[2][0] for edge in cycle[:-1]);
					nxi_moves = "".join(edge[2][1] for edge in cycle[:-1]);
					lambda_value = max(
						approx_value(
							to_continued_fraction(moves, cyclic=True)
						) for moves in (xi_moves, nxi_moves)
					)
					lambda_new_N = max(
						approx_value(to_continued_fraction(moves, cyclic=True),
							fast = True
						) for moves in (xi_moves, nxi_moves)
					)
					if lambda_new_N < lambda_limit:
						if verbosity >= 3/2:
							print("In endgame, found new best cycle of type", [
								to_continued_fraction(xi_moves, cyclic=True),
								to_continued_fraction(nxi_moves, cyclic=True)],
								"with approximability", lambda_value, "=", lambda_new_N
							);
						if lambda_limit != lambda_value:
							excluded = min(excluded, lambda_limit)
						lambda_limit = lambda_value
					cycle_qualities[tuple(cycle)] = lambda_value
				if verbosity >= 1:
					print("Lowest point of Ln is:" + " "*13, N(lambda_limit), "=", lambda_limit)
					print("Next lowest point of Ln is at least", N(excluded), "=", excluded)
					# print(cycle_qualities)
					for entry in cycle_qualities.items():
						if entry[1] == lambda_limit:
							cycle = entry[0]
							cycle_start = cycle[0][0]
							if cycle_start == 0:
								xi_transient = ""
								nxi_transient = ""
							else:
								transient = graph.edge_label(0, cycle_start)[0]
								xi_transient = transient[0]
								nxi_transient = transient[1]
							xi_period = "".join(edge[2][0] for edge in cycle[:-1]);
							nxi_period = "".join(edge[2][1] for edge in cycle[:-1]);
							
							xi_cf = continued_fraction(
								to_eventually_periodic_cf(xi_transient, xi_period))
							nxi_cf = continued_fraction(
								to_eventually_periodic_cf(nxi_transient, nxi_period))
							xi = xi_cf.value(); nxi = nxi_cf.value()
							assert nxi == n * xi
							
							print("Minimum point at cycle", [edge[0] for edge in entry[0]])
							print("\txi =", xi, "=", xi_cf)
							print("\tp*xi =", nxi, "=", nxi_cf)
				return lambda_limit, excluded
		
		prev = current
		current = set()
	if verbosity >= 1:
		print("Lowest point found is", lambda_limit, "=", N(lambda_limit))
	return None

markoff_data = {
	(5,) : sqrt(5),
	(2,) : sqrt(8),
	(13,17) : sqrt(221)/5,
	(37,41) : sqrt(1517)/13
}
	
def is_first_6(p):
	return any(all(kronecker_symbol(p, m) >= 0
		for m in ent) for ent in markoff_data)

def min_Ln_easy(p):
	for (ent, val) in markoff_data.items():
		if all(kronecker_symbol(p, m) >= 0 for m in ent):
			return val

def test_all_primes(verbosity = 1):
	for p in Primes():
			print("-"*40)
			find_min_Ln(p, verbosity = verbosity);

def test_interesting_primes(verbosity = 2):
	for p in Primes():
		if not is_first_6(p):
			print("-"*40)
			find_min_Ln(p, verbosity = verbosity);

def print_boring_primes():
	for p in Primes():
		if p > 2000:
			break
		val = min_Ln_easy(p)
		if val is not None:
			print(p, N(val))
			
def compact(cf):
	return "".join([str(a) for a in cf])

#### Testing:
# show(get_graph(7))
# find_good_paths(67, 3.67821976, 30)
# find_min_Ln(67, verbosity=2)
# test_interesting_primes()

# An:
# for p in Primes():
#		if not is_first_6(p):
#			print("p =", p)
#			find_min_Ln(p, verbosity = 1) 