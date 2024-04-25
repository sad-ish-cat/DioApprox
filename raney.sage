load("lattice.sage")

matrix_L = matrix([[1,0], [1,1]])
matrix_R = matrix([[1,1], [0,1]])
matrix_L_inv = matrix([[1,0], [-1,1]])
matrix_R_inv = matrix([[1,-1], [0,1]])

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
		if mat[1][0] >= mat[1][1]:
			word = "L" + word;
			mat = mat * matrix_L_inv;
		elif mat[0][1] >= mat[0][0]:
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
		if mat[0][1] >= mat[1][1]:
			word = word + "R";
			mat = matrix_R_inv * mat;
		elif mat[1][0] >= mat[0][0]:
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
		return (GF(n)(mat[0][0]) / GF(n)(mat[1][0] - mat[0][0])).lift();
	except ZeroDivisionError:
		return (GF(n)(mat[0][1]) / GF(n)(mat[1][1] - mat[0][1])).lift();
	
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
def find_good_paths(n, lambda_limit, max_length = 300):
	
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
					if (min_lambda_new(tuple(to_continued_fraction(xi_moves))) <= lambda_limit
						and min_lambda_new(tuple(to_continued_fraction(nxi_moves))) <= lambda_limit):
							# print(xi_moves, nxi_moves);
							current.add(newpath)
		prev = current
		print(len(prev), "good paths of length", k)
		current = set()
	return prev
	
#### Testing:
# show(get_graph(7))
# find_good_paths(67, 3.67821976, 30)