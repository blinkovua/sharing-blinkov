#!/usr/bin/env python
# coding: utf-8

# In[1]:


from IPython.display import display, Image
from pprint import pprint
import json

from ginv import *
sympy.init_printing()


# The database of polynomial systems.
# 
# The demonstration database is an essential part of PHCpack
# 
# http://homepages.math.uic.edu/~jan/

# In[2]:


with open('poly_exmpl1.json', 'r') as f:
      poly_exmpl = json.loads(f.read())


# In[3]:


" ".join(sorted(poly_exmpl.keys()))


# In[4]:


" ".join(sorted(k for k, v in poly_exmpl.items()    if 'timing' in v and v['timing']['GinvBlockLow'] > 1.0))


# In[5]:


" ".join(sorted(k for k, v in poly_exmpl.items()    if len(v['eqs']) < 6))


# In[6]:


ex = poly_exmpl['cyclic7']
ex


# In[7]:


Monom.cmp = Monom.TOPdeglex
var = ex['var']
fun = []
Monom.init(var, fun)
for var_i, var_g in enumerate(var):
    globals()[var_g] = Poly(Monom(var_i))
# invdiv = JanetCache()
invdiv = Janet()
res = ginvBlockLow([eval(eq) for eq in ex['eqs']], invdiv, level=3)
print(f"crit1: {res[1]}")
print(f"crit2: {res[2]}")
print(f" time: {res[0]:.2f} sec")
print(f"    count: {invdiv.count}")
print(f"reduction: {invdiv.reduction}")
print(f"       HP: {invdiv.HP()}")
invdiv.saveImage('invdiv.pdf', level=1)
invdiv.saveImage('invdiv.png', level=1)
Image('./invdiv.png')


# In[7]:


# for w in invdiv:
#     print(w.poly)


# In[8]:


for eq in ex['eqs']:
    print(eq)
    print(eval(eq))


# In[9]:


eval("-a**2*c+a**2")


# In[ ]:




