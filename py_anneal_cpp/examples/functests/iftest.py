"""
Test how the sympy ccode converter handles a conditional statement.
"""

problem_name = 'iftest'

stims = ()

syms = []
X_list = ['x']
NX = len(X_list)
for X in X_list:
    symbol_ptr = sym.symbols(X)
    tmp = '{0} = symbol_ptr'.format(X)
    exec tmp
    syms.append(symbol_ptr)
