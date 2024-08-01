# load("approx.sage")
from sage.misc.search import search

def words():
  for b in [1..5]:
    for a in [0..b]:
      if gcd(a,b) > 1:
        continue
      word = "".join(["22" if floor((i+1)*a/b) > floor(i*a/b) else "11"
        for i in range(b)]);
      print("")
      print(word);
      # print(d_quality(word, verb=True))
    
    
def flipz(sol):
  x,y,z = sol
  return 3*x*y - z, x, y
def flipy(sol):
  x,y,z = sol
  return 3*x*z - y, x, z
  
def find_Markoff_in_order(max_num):
  sols = [(1,1,1)];
  for i in (0..):
    if len(sols) > i:
      sol = sols[i];
    else:
      return sols;
    for new in {flipz(sol), flipy(sol)}:
      if new[0] <= max_num:
        found, j = search(sols, new);
        assert not found;
        sols.insert(j, new)
        
def check_strong_Markoff_uniqueness(max_num):
  nums = [(9*sol[0]^2 - 4).squarefree_part()
    for sol in find_Markoff_in_order(max_num)];
  nums.sort();
  for i in [0..len(nums) - 2]:
    if nums[i] == nums[i+1]:
      print("Duplicate value", nums[i]);
      return False
  return True
  
def compute_class_groups(max_num):
  for sol in find_Markoff_in_order(max_num):
    # print(sol[0])
    D = (9*sol[0]^2 - 4).squarefree_part()
    # print("D =", D)
    print(D, end=" ", flush=True);
    print(list(QuadraticField(D).class_group().elementary_divisors()))
    
    
def real_qf_class_groups(min_rad, max_rad):
  for D in [min_rad..max_rad]:
    if D.is_squarefree():
      print(D, end=" ", flush=True);
      print(list(QuadraticField(D).class_group().elementary_divisors()))