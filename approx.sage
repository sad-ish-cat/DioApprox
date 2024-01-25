# approx.sage: A program for computing badly approximable points on conics,
# leading to low-lying elements of the Lagrange and Markoff n-spectra.
load("ComparableMixin.sage")

from sage.rings.continued_fraction import last_two_convergents

# First version.
def find_badly_approximable_numbers(nmax, limit = 3, verb=False):
  result = [];
  for c in [-nmax .. -1]:
    print(-c)
    for a in [1 .. -c]:
      for b in [-a .. 0]:
        d = b^2 - 4*a*c;
        if d.is_square() or gcd([a,b,c]) > 1:
          continue;
        form = a*x^2 + b*x + c;  
        if verb:
          print("d =", d)
          print("form is", form)
        K.<sqrtd> = QuadraticField(d);
        xi_raw = (-b + sqrtd)/(2*a);
        xibar_raw = (-b - sqrtd)/(2*a);
        xi = xi_raw - ceil(xibar_raw) # to get purely periodic
        cfrac = continued_fraction(xi);
        if verb:
          print(cfrac)
        assert cfrac.preperiod() == (); 
        period = cfrac.period()
        l = len(period);
        # Find the biggest partial quotient(s).
        m = max(period);
        quality = 0;
        for i in range(l):
          if period[i] == m:
            conv = cfrac.convergent((i-1) % l);
            Zquality = a * ((xi - conv)*conv.denominator()).norm();
            if verb:
              print("off by", Zquality);
            quality = max(quality, sqrtd/abs(Zquality));
            if quality > limit:
              break;
        if quality <= limit:
          result.append((RR(quality), d, cfrac, form));
  
  result.sort();
  return result;

def is_nonstandard(period):
  new_period = period;
  l = len(period);
  for it in range(l-1):
    new_period = new_period[1:l] + new_period[0:1];
    if new_period <= period:
      return True;
  return False;
  
# Second version.
def examine_ctd_fracs(lmax, max_term = 2, limit = 3, verb=True):
  result = [];
  for l in [1..lmax]:
    if verb:
      print("length", l)
    for period in Tuples([1..max_term], l):
      if is_nonstandard(period):
        continue;
      
      cfrac = continued_fraction([(), period])
      xi = cfrac.value();
      form = xi.minpoly();
      form = form * form.denominator();
      a = form.lc();
      d = form.discriminant();
      sqrtd = d.sqrt();
      
      # Find the biggest partial quotient(s).
      m = max(period);
      quality = 0;
      for i in range(l):
        if period[i] == m:
          conv = cfrac.convergent((i-1) % l);
          Zquality = a * ((xi - conv)*conv.denominator()).norm();
          if verb:
            print("off by", Zquality);
          quality = max(quality, sqrtd/abs(Zquality));
          if quality > limit:
            break;
      if quality <= limit:
        result.append((RR(quality), d, cfrac));
  
  result.sort();
  return result;

# Computes the disc & approx'ty of a periodic ctd frac.
def d_quality(period, limit=Infinity, verb=False):
  if isinstance(period, str):
    period = tuple(ZZ(i) for i in period);
  
  Z_poly.<xx> = PolynomialRing(ZZ)
  cfrac = continued_fraction(period)
  p0, q0, p1, q1 = last_two_convergents(period)
  form = q1 * xx^2 + (q0 - p1) * xx - p0
  if verb:
    print(form);
  a = form.lc();
  d = form.discriminant();
  
  # Find the biggest partial quotient(s).
  l = len(period);
  m = max(period);
  quality = 0;
  for i in range(l):
    if period[i] == m:
      conv = cfrac.convergent((i-1) % l);
      Zquality = a * conv.denominator()^2 * form(xx = conv) / form.lc();
      if verb:
        print("off by", Zquality);
      quality = max(quality, d/(Zquality)^2);
      if quality >= limit^2:
        # It's bad.
        return (d, limit^2);
        
  return (d, quality);


# Third version
def n_conic(n, lmax=1, max_term=Infinity, limit=Infinity, choices=None,
  verb=False):
  
  Z_poly.<xx> = PolynomialRing(ZZ)
  result = [];
  if n == 1:
    largest_prime_power = 1
  else:
    largest_prime_power = max(d for d in divisors(n) if d.is_prime_power());
  Z2 = FreeModule(ZZ, 2)
  
  for l in [1..lmax]:
    if verb and not(choices):
      print("length", l)
    for period1 in choices if choices else Tuples([1..max_term], l):
      if is_nonstandard(period1) and not(choices):
        continue;
      if verb and choices:
        print(period1)
      p0, q0, p1, q1 = last_two_convergents(period1)
      form = q1 * xx^2 + (q0 - p1) * xx - p0
      c, b, a = form.coefficients()
      assert b^2 - 4*a*c == form.disc();
      
      d1, quality1 = d_quality(period1);
      
      if quality1 < limit^2:
        # Now apply the various n-transformations.
        for dd in divisors(n):
          if not(largest_prime_power.divides(a*dd^2)):
            continue;
          d = ZZ(n/dd);
          for i in range(dd):
            if gcd([d, dd, i]) > 1:
              continue;
            if not(largest_prime_power.divides(2*a*dd*i) and
                largest_prime_power.divides(a*i^2 - b*d*i + c*d^2)):
              continue;
          
            # Test the resulting k*xi companions
            quality = quality1;
            ds = [MarkoffNumber(d1)];
            periods = [CompactSequence(period1)]
            for k in divisors(n)[1:]:
              M = Z2.submodule([(k,0),(0,k),(dd,0),(i,d)]);
              u2, u1 = M.0;
              u4, u3 = M.1;
              form2 = Z_poly(form(xx = (u4*xx - u2)/(-u3*xx + u1))
                * (-u3*xx + u1)^2)
              form2 = Z_poly(form2 / form2.content())
              
              cfrac2 = continued_fraction(
                NumberField(form2,names='gen',check=False,
                  maximize_at_primes = []).0);
              period2 = cfrac2.period();
              d2, quality2 = d_quality(period2);
              
              if quality2 < limit^2:
                ds += [MarkoffNumber(d2)]
                periods += [CompactSequence(period2)]
                quality = max(quality, quality2)
              else:
                quality = limit^2;
                break;
            if quality < limit^2:
              if len(set(ds)) == 1:
                ds = [ds[0]];
              result += [(RR(sqrt(quality)), ds, periods)]
  try:
    result.sort(key = lambda x : x[0]); # Sort by quality only
  except:
    if verb:
      print("Could not sort")
  return result;
  
class MarkoffNumber(ComparableMixin):
  def __init__(self, d):
    self.d = d
    
  def __repr__(self):
    d = self.d;
    if d.is_square():
      return str(d.sqrt()) + "^2";
    elif (d + 1).is_square():
      return str((d + 1).sqrt()) + "^2 - 1";
    elif (d + 4).is_square():
      return str((d + 4).sqrt()) + "^2 - 4";
    else:
      return str(d);
      
  def _cmpkey(self):
    return (self.d, self.d);
    
class CompactSequence(ComparableMixin):
  def __init__(self, seq):
    self.seq = seq
  def __repr__(self):
    seq = self.seq;
    if all(len(str(item)) == 1 for item in seq):
      return "(" + "".join(str(item) for item in seq) + ")";
    else:
      return repr(seq);