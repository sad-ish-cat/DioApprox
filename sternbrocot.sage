# Class for recursive binary trees
# left - Node - left descendant of tree
# right - Node - right descendant of tree
# data - data stored in node
# height - height of tree
class Node:

	def __init__(self, left = None, right = None, data = None, height = 0):
		self.left = left
		self.right = right
		self.data = data
		self.height = 0


# Makes a Stern-Brocot tree of given height and initial root data
# height - int - height of tree
# root_data - [[a, b], [c, d], [e, f]] a,b,c,d,e,f int - node data for the start of the tree
def make_tree(height, root_data = ((0,1),(1,1),(1,0))):

	if height == 1:
		root = Node(data = root_data)
	else:
		root = Node(data = root_data)
		root.height = height
		root.left = make_tree(height-1, root_data = (root_data[0],(root_data[0][0]+root_data[1][0],root_data[0][1]+root_data[1][1]),root_data[1]))
		root.right = make_tree(height-1, root_data = (root_data[1],(root_data[2][0]+root_data[1][0],root_data[2][1]+root_data[1][1]),root_data[2]))
		
	return root
		
# Helper function for display_tree
def convert_to_labelled_binary_tree(tree):

	if tree.left == None:
		return LabelledBinaryTree([None, None], label=str(tree.data[1][0]) + "/" + str(tree.data[1][1]))
		
		
	binary_tree = LabelledBinaryTree([convert_to_labelled_binary_tree(tree.left), convert_to_labelled_binary_tree(tree.right)], label=str(tree.data[1][0]) + "/" + str(tree.data[1][1]))

	return binary_tree
	
# Given a tree, draws the tree with ASCII art
# tree - Node - tree to be drawn
def display_tree(tree):

	return ascii_art(convert_to_labelled_binary_tree(tree))
	

# Given a numberator, denominator, and reference tree this will attempt to find the LR sequence in this tree within the maximum number of steps
# num - int - numerator of fraction of number to be found
# den - int - denominator of fraction of number to be found
# root_data - [[a, b], [c, d], [e, f]] a,b,c,d,e,f int - node data for the start of the tree
# max_steps - int - maximum number of left-right moves the function will take before terminating and outputting the LR sequence up to that point
def find_left_right (num, den, root_data = [[0,1],[1,1],[1,0]], max_steps = 10):

	if max_steps == 0:
		return " rest of sequence not found within number of max_steps"

	if (root_data[1][0]/root_data[1][1]) == (num/den):
		return ""
		
	if (root_data[1][0]/root_data[1][1]) > (num/den):
		return "L" + find_left_right(num, den, [root_data[0],[root_data[0][0]+root_data[1][0],root_data[0][1]+root_data[1][1]],root_data[1]], max_steps - 1)
		
	if (root_data[1][0]/root_data[1][1]) < (num/den):
		return "R" + find_left_right(num, den, [root_data[1],[root_data[2][0]+root_data[1][0],root_data[2][1]+root_data[1][1]],root_data[2]], max_steps - 1)
		
# Given a decimal number and reference tree this will attempt to find the LR sequence in this tree within the maximum number of steps
# dec - float - number to be found
# root_data - [[a, b], [c, d], [e, f]] a,b,c,d,e,f int - node data for the start of the tree
# max_steps - int - maximum number of left-right moves the function will take before terminating and outputting the LR sequence up to that point
def find_left_right_dec (dec, root_data = [[0,1],[1,1],[1,0]], max_steps = 10):

	if max_steps == 0:
		return ""

	if (root_data[1][0]/root_data[1][1]) == dec:
		return ""
		
	if (root_data[1][0]/root_data[1][1]) > dec:
		return "L" + find_left_right_dec(dec, [root_data[0],[root_data[0][0]+root_data[1][0],root_data[0][1]+root_data[1][1]],root_data[1]], max_steps - 1)
		
	if (root_data[1][0]/root_data[1][1]) < dec:
		return "R" + find_left_right_dec(dec, [root_data[1],[root_data[2][0]+root_data[1][0],root_data[2][1]+root_data[1][1]],root_data[2]], max_steps - 1)
		
		
# This will return the data found by traversing a reference tree by the given LR sequence
# left_right_sequence - str - String containing the letters L and R
# root_data - [[a, b], [c, d], [e, f]] a,b,c,d,e,f int - node data for the start of the tree
def create_data_from_left_right_old (left_right_sequence, root_data = [[0,1],[1,1],[1,0]]):

	letters = list(left_right_sequence)
	if len(letters) == 0:
		return root_data
	
	if letters[0] == "L": 
		return create_data_from_left_right(left_right_sequence[1:], [root_data[0],[root_data[0][0]+root_data[1][0],root_data[0][1]+root_data[1][1]],root_data[1]])
	
	elif letters[0] == "R":
		return create_data_from_left_right(left_right_sequence[1:], [root_data[1],[root_data[2][0]+root_data[1][0],root_data[2][1]+root_data[1][1]],root_data[2]])
		
@cached_function
def create_data_from_left_right (left_right_sequence, root_data = ((0,1),(1,1),(1,0))):

	letters = list(left_right_sequence)
	
	output_data = root_data
	
	for letter in letters:
		if letter == "L":
			output_data = (output_data[0],(output_data[0][0]+output_data[1][0],output_data[0][1]+output_data[1][1]),output_data[1])
		elif letter == "R":
			output_data = (output_data[1],(output_data[2][0]+output_data[1][0],output_data[2][1]+output_data[1][1]),output_data[2])
	
	return output_data
	
@cached_function
def create_number_from_left_right (left_right_sequence, root_data = ((0,1),(1,1),(1,0))):
	
	if left_right_sequence == "":
		return 0
	
	output_data = create_data_from_left_right(left_right_sequence, root_data = root_data)
	return output_data[1][0]/output_data[1][1]
		
# Helper function for find_left_right_relative, returns true if the interval corresponding to A is a subset of the one corresponding to B

def is_subset(A, B):

	if(B[2][1] == 0) and (A[0][0]/A[0][1] >= B[0][0]/B[0][1]):
		return True
	
	elif(B[2][1] == 0):
		return False
		
	elif(A[2][1] == 0):
		return False
		
	elif (A[0][0]/A[0][1] >= B[0][0]/B[0][1]) and (A[2][0]/A[2][1] <= B[2][0]/B[2][1]):
		return True
		
	return False
	
# Find the sequence of left-right moves in a tree with respect to a reference tree
# data - [[a, b], [c, d], [e, f]], a,b,c,d,e,f int - node that the function will attempt to find an LR sequence in the tree that leads to this node
# root_data - [[a, b], [c, d], [e, f]] a,b,c,d,e,f int - node data for the start of the tree
# max_steps - int - maximum number of left-right moves the function will take before terminating and outputting the LR sequence up to that point
@cached_function
def find_left_right_relative (data, root_data = ((0,1),(1,1),(1,0)), max_steps = 10):
	
	T = make_tree(2, root_data)
 
	if max_steps == 0:
		return " rest of sequence not found within number of max_steps"
		
	if(is_subset(data, T.right.data)):
		return "R" + find_left_right_relative(data, T.right.data, max_steps - 1)
		
	if(is_subset(data, T.left.data)):
		return "L" + find_left_right_relative(data, T.left.data, max_steps - 1)
		
	return ""
		
		
		
		

		
		






		
