load("lattice.sage")
import time
import random

	
	
def test_find_good_paths(n, lambda_limit, max_length = 5):
		
	start = time.time()
	
	find_good_paths(n, lambda_limit, max_length)
	
	end = time.time()
	
	return end - start
	
def test_brocot(iterations):

	timer = 0
	
	for i in range(iterations):
		lrseq = ""
		for j in range(20):
			lrseq += random.choice(["L", "R"])
		

		start = time.time()
		data = create_data_from_left_right_old(lrseq)
		a = data[1][0]/data[1][1]
		end = time.time()
		timer += end - start
		
	return timer
	
	
def test_ctd(iterations):

	timer = 0
	
	for i in range(iterations):
		lrseq = ""
		for j in range(20):
			lrseq += random.choice(["L", "R"])
		

		start = time.time()
		
		ctdfrac = to_continued_fraction(lrseq)
		
		a = continued_fraction(ctdfrac).value()
		
		end = time.time()
		timer += end - start
		
	return timer
	

	