# modular_group.sage: A program for computing examples of modular groups
# from Vulakh's paper.

RRR.<x> = PolynomialRing(ZZ)
G = Gamma0(13)
Farey = FareySymbol(G);

# for n in [0..8]:
#   print(Farey.word_problem(G.gen()^n))

# Try fixed points of random short words.
gens = G.gens()

def print_form_data(M):
  F = (matrix([1,x])*M*matrix([[x],[-1]]))[0,0]
  print(F)
  period = NumberField(F, names = "rtF").0.continued_fraction().period()
  print(period)
  return period

# Get a tree from three involutions.
def tree(S1, S, S2, num_iter=7):
  S1 = GL(2,ZZ)(S1)
  S  = GL(2,ZZ)(S )
  S2 = GL(2,ZZ)(S2)
  
  assert (S1^4).is_one(); 
  assert (S ^4).is_one(); 
  assert (S2^4).is_one(); 
  
  print("T =\n", S1*S*S2, sep="")
  
  found_ctd_fracs = []
  found_ctd_fracs += [print_form_data(S1*S), print_form_data(S*S2),
    print_form_data(S1*S2)]
  old_triples = [(S1, S, S2)]
  for i in range(num_iter):
    print("----------")
    new_triples = [];
    for U1, U, U2 in old_triples:
      UU1 = U1*U*U1
      new_triples += [(UU1, U1, U2)];
      found_ctd_fracs += [print_form_data(UU1*U2)]
      UU2 = U2*U*U2
      new_triples += [(U1, U2, UU2)];
      found_ctd_fracs += [print_form_data(U1*UU2)]
    old_triples = new_triples
  return found_ctd_fracs