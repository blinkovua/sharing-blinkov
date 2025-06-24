from sympy import *

from monom import *

def init(f, v, h):
    global pda_f, pda_v, pda_h, pda_fun, pda_var
    # pda_f, pda_v, pda_h = f, v, h
    # Monom.variables = 2
    # Monom.cmp = Monom.TOPdeglex
    # Monom.zero = Monom(0 for v in range(Monom.variables))
    # pda_fun = dict(zip(pda_f,\
    #             (Monom(0 if v else i for v in range(Monom.variables))\
    #             for i in range(1, 1 + len(pda_f)))))
    # pda_var = {pda_v: Monom((0, 1))}

    pda_f, pda_v, pda_h = f, v, h
    Monom.cmp = Monom.TOPdeglex
    Monom.init([str(pda_v)], [str(f) for f in pda_f])
    pda_var = Monom(0)
    pda_fun = dict(zip(pda_f, (Monom(pos=i) for i in range(len(pda_f)))))

def set_clip(n, clp, p):
    global pda_n, pda_clp, pda_p
    pda_n, pda_clp, pda_p = n, clp, p

# def T(f, j1):
#     global pda_n, pda_p
#     return sum(\
#         diff(f, pda_v, j)*(pda_h*(j1 - pda_p))**j/\
#                  (factorial(j))\
#         for j in range(pda_n))
def T(f):
    global pda_n, pda_v, pda_h, pda_p
    g = f.subs(pda_v, pda_v - pda_p*pda_h)
    return sum(\
        diff(g, pda_h, j).subs(pda_h, 0).doit()*pda_h**j/\
                 (factorial(j))
        for j in range(pda_n))

def clip(f):
    f = f.expand()
    return [f.coeff(pda_h, i) for i in range(pda_clp)]

def clip_n(f):
    f = f.expand()
    return sum(f.coeff(pda_h, i)*pda_h**i for i in range(pda_n))

def df2m(a):
    if a.func != Derivative:
        m = pda_fun[a]
    else:
        m = pda_fun[a.args[0]]
        for xi in a.args[1:]:
            if isinstance(xi, Symbol):
                m = m*pda_var
            else:
                m = m*pda_var**xi[1]
    return m

def m2df(m):
    return pda_f[m.position()].diff(pda_v, m[0])

def findDiv(a, d):
    r = None
    def find(a, r):
        if a.args:
            if (a in pda_fun) or \
                    (a.func == Derivative and a.args[0] in pda_fun):
                m = df2m(a)
                if m.divisible(d) and (not r or m.cmp(r) > 0):
                    r = m
            else:
                for s in a.args:
                    r = find(s, r)
        return r
    return find(a, r)

def reduction(f1, f2, m, c, shift):
    assert shift < pda_clp
    r = [f1[i] for i in range(shift)]
    if not m:
        for i in range(shift, pda_clp):
            r.append(expand(f1[i] - f2[i-shift]*c))
    else:
        for i in range(shift, pda_clp):
            r.append(expand(f1[i] - f2[i-shift].diff(*m)*c))
    return r

def NF(f, df, G, head=False):
    assert len(df) == len(G)
    ms = [df2m(d) for d in df]
    for i in range(0 if head else 1, pda_clp):
        t = 0
        if f[i]:
            while True:
                r = None
                for l in range(len(ms)):
                    r = findDiv(f[i], ms[l])
                    if r:
                        break
                if not r:
                    break
                m = m2df(r)
                deg = f[i].as_poly(m).degree()
                c = f[i].coeff(m, deg)
                if deg > 1:
                    c = f[i].coeff(m, deg)*m**(deg - 1)
                m = r/ms[l]
                d = []
                if m[0] > 0:
                    d.append(pda_v)
                    if m[0] > 1:
                        d.append(m[0])

                f = reduction(f, G[l], tuple(d), c/G[l][0].coeff(df[l]), i)
    return f

def compact(f):
    def cmpct(a):
        if not a.args:
            return a
        else:
            if a in pda_fun:
                return Symbol("%s" % a.func, real=True)
            elif a.func == Derivative and a.args[0] in pda_fun:
                m = []
                for xi in a.args[1:]:
                    if isinstance(xi, Symbol):
                        m.append(str(xi))
                    else:
                        m.append(str(xi[0])*xi[1])
                return Symbol("%s_{%s}" % (a.args[0].func, "".join(m)), real=True)
            else:
                return a.func(*tuple(cmpct(s) for s in a.args))
    return cmpct(f)

def prn(a, p=None):
    if p:
        for i in range(pda_clp):
            print("%s^%d =>" % (pda_h,  i))
            display(compact(a[i]).collect(p, factor))
    else:
        for i in range(pda_clp):
            print("%s^%d =>" % (pda_h,  i))
            display(compact(a[i]).factor())


def prnlatex(a, p=None):
    if p:
        print(latex(compact(a[0]).collect(p, factor)))
        print(r"+%s\left(" % pda_h)
        print(latex(compact(a[1]).collect(p, factor)))
        print(r"\right)")
        for i in range(2, pda_clp):
            print(r"+%s^%d\left(" % (pda_h, i))
            print(latex(compact(a[i]).collect(p, factor)))
            print(r"\right)")
    else:
        print(latex(compact(a[0]).factor()))
        print(r"+%s\left(" % pda_h)
        print(latex(compact(a[1]).factor()))
        print(r"\right)")
        for i in range(2, pda_clp):
            print(r"+%s^%d\left(" % (pda_h, i))
            print(latex(compact(a[i]).factor()))
            print(r"\right)")


if __name__ == '__main__':
    eps, t, h = symbols(r'\varepsilon, t, h', real=True)
    u, u1 = symbols('u, u1', cls=Function)

    eq = u(t).diff(t, 2)  + u(t) + eps*u(t)**3
    print(eq)
    print()

    eq1 = (u(t).diff(t)**2/2 + u(t)**2/2 + eps*u(t)**4/4) - h
    print(expand(eq1.diff(t) - eq*u(t).diff(t)))
    print()

    init((u(t), u1(t)), t, h)
    set_clip(8, 7, Rational(0, 1))

    print(expand(T(u(t+h))-T(u(t)))/(2*h))
    print()

    def RungeKutta4(f, y):
        k1 = f(y)
        k2 = f(y + h*k1/2)
        k3 = f(y + h*k2/2)
        k4 = f(y + h*k3)
        return expand((k1 + 2*k2 + 2*k3 + k4)/6)

    def f(y):
        return Matrix([\
            clip_n(y[1]),\
            clip_n(-y[0] - eps*y[0]**3),\
        ])

    set_clip(7, 6, Rational(0, 1))
    r = RungeKutta4(f, Matrix([u(t), u1(t)]))

    F1 = clip((T(u(t+h))-T(u(t)))/h - T(r[0]))
    prnlatex(F1, eps)
    print()

    F2 = clip((T(u1(t+h))-T(u1(t)))/h - T(r[1]))
    prnlatex(F2, eps)
    print()

    f1 = NF(F1, [u(t).diff(t), u1(t).diff(t)], [F1, F2], head=False)
    prnlatex(f1, eps)
    print()

    f2 = NF(F2, [u(t).diff(t), u1(t).diff(t)], [F1, F2], head=False)
    prnlatex(f2, eps)
    print()

    F3 = clip((T(u1(t+h))**2/2 + T(u(t+h))**2/2 + eps*T(u(t+h))**4/4) -\
          (T(u1(t))**2/2 + T(u(t))**2/2 + eps*T(u(t))**4/4))
    prnlatex(F3, eps)
    print()

    f3 = NF(F3, [u(t).diff(t), u1(t).diff(t)], [f1, f2], head=False)
    prnlatex(f3, eps)
    print()
