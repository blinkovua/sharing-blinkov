# -*- coding: utf-8 -*-

class Monom(tuple):
  variables = None
  zero = None

  def __new__ (cls, *args):
    if len(args) == 0:
      return Monom.zero
    else:
      return super().__new__(cls, *args)

  def __init__(self, *args):
    assert len(self) == Monom.variables
    assert all(v >= 0 for v in self)

  def position(self):
    return self[0]

  def degree(self):
    return sum(self[1:])

  def lcm(self, other):
    assert self[0] == other[0] or other[0] == 0
    return Monom(max(self[i], other[i]) for i in range(Monom.variables))

  def gcd(self, other):
    assert self[0] == other[0] or other[0] == 0
    return Monom(min(self[i], other[i]) for i in range(Monom.variables))

  def __nonzero__(self):
    return sum(self[1:]) > 0

  def prolong(self, var):
    assert 1 <= var < Monom.variables
    return Monom(self[i] + (1 if i == var else 0) for i in range(Monom.variables))

  def __mul__(self, other):
    assert self[0] == 0 or other[0] == 0
    return Monom(self[i] + other[i] for i in range(Monom.variables))

  def divisible(self, other):
    if self[0] != other[0] and other[0] != 0:
      return False
    for i in range(1, Monom.variables):
      if self[i] < other[i]:
        return False
    return True

  def divisibleTrue(self, other):
    if self[0] != other[0] and other[0] != 0:
      return False
    i = 1
    while i < Monom.variables and self[i] == other[i]:
      i += 1
    if i == Monom.variables or self[i] < other[i]:
      return False
    else:
      i += 1
      while i < Monom.variables and self[i] >= other[i]:
        i += 1
      return i == Monom.variables

  def __truediv__(self, other):
    assert self.divisible(other)
    return Monom(self[i] - other[i] for i in range(Monom.variables))

  def __pow__(self, other):
    assert other >= 0
    return Monom(self[i]*other if i > 0 else self[i] for i in range(Monom.variables))

  def __lt__(self, other):
    return self.cmp(other) < 0

  def __eq__(self, other):
    return self.cmp(other) == 0

  def __hash__(self):
    return tuple.__hash__(self)

  def cmp(self, other):
    assert False

  def __lex(self, other):
    for i in range(1, Monom.variables):
      if self[i] > other[i]:
        return 1
      elif self[i] < other[i]:
        return -1
    return 0

  def __deglex(self, other):
    d1, d2 = self.degree(), other.degree()
    if d1 > d2:
      return 1
    elif d1 < d2:
      return -1
    else:
      for i in range(Monom.variables-1, 0, -1):
        if self[i] < other[i]:
          return 1
        elif self[i] > other[i]:
          return -1
      return 0

  def __alex(self, other):
    d1, d2 = self.degree(), other.degree()
    if d1 > d2:
      return -1
    elif d1 < d2:
      return 1
    else:
      for i in range(1, Monom.variables):
        if self[i] > other[i]:
          return 1
        elif self[i] < other[i]:
          return -1
      return 0

  def POTlex(self, other):
    return (lambda a, b: (a > b) - (a < b))(self[0], other[0]) \
            or self.__lex(other)

  def TOPlex(self, other):
    return self.__lex(other) or \
           (lambda a, b: (a > b) - (a < b))(self[0], other[0])

  def POTdeglex(self, other):
    return (lambda a, b: (a > b) - (a < b))(self[0], other[0]) \
           or self.__deglex(other)

  def TOPdeglex(self, other):
    return self.__deglex(other) or \
      (lambda a, b: (a > b) - (a < b))(self[0], other[0])

  def POTalex(self, other):
    return (lambda a, b: (a > b) - (a < b))(self[0], other[0]) \
           or self.__alex(other)

  def TOPalex(self, other):
    return self.__alex(other) or \
      (lambda a, b: (a > b) - (a < b))(self[0], other[0])

if __name__ == '__main__':
  variables = ['a', 'b', 'c', 'd', 'e', 'f']
  Monom.variables = len(variables)+1
  # POTlex TOPlex POTdeglex TOPdeglex POTalex TOPalex
  Monom.cmp = Monom.POTlex
  Monom.zero = Monom(0 for v in range(Monom.variables))
  for i in range(len(variables)):
    globals()[variables[i]] = Monom(0 if l-1 != i else 1 for l in range(Monom.variables))

  print(Monom())
  m1 = b*c**2
  print(m1)
  print(m1/m1)
  print(repr(m1))

  m2 = a*e*f**2 * Monom()
  print(m2, m2.prolong(3))
  print(m1.POTlex(m2))
  print(m1.TOPdeglex(m2))
  print(m1.TOPalex(m2))
  print(m1 == m2)
  print(m1 < m2 <= m2)
  a = set([m1, m2, a, b, c])
  print(a)

