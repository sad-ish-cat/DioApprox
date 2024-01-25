load("approx.sage")

def words():
  for b in [1..5]:
    for a in [0..b]:
      if gcd(a,b) > 1:
        continue
      word = "".join(["22" if floor((i+1)*a/b) > floor(i*a/b) else "11"
        for i in range(b)]);
      print("")
      print(word);
      print(d_quality(word, verb=True))
    
    
def flipx(sol):
  x,y,z = sol
  return 2*y*z - x, y, z
def flipy(sol):
  x,y,z = sol
  return x, 3*x*z - y, z
def flipz(sol):
  x,y,z = sol
  return x, y, 6*x*y - z
  
def find_Markoff_6(iterns = 10):
  sols = {(1,1,1)}
  for i in range(iterns):
    old = list(sols);
    for sol in old:
      sols.update({flipx(sol), flipy(sol), flipz(sol)})
  return sols;