import sys,itertools,json
from pathlib import Path
import sympy as S

I=S.I;t=S.symbols('t',real=True);c=(1-t*t)/(1+t*t);d=2*t/(1+t*t)
def G(v):
 p,q,z=v;return S.Matrix([[0,p,q,z],[p,0,I*z,-I*q],[q,-I*z,0,I*p],[z,I*q,-I*p,0]])
x=S.symbols('x:12',real=True);T=S.zeros(4)
for k,(i,j) in enumerate([(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]):T[i,j]=x[2*k]+I*x[2*k+1];T[j,i]=S.conjugate(T[i,j])
V=G([0,-d,c]);eq=[]
for z in V*T-T:
 z=S.expand(z);eq.extend([S.together(S.re(z)).as_numer_denom()[0],S.together(S.im(z)).as_numer_denom()[0]])
sol=S.linsolve(eq,x);print('supported',sol,flush=True)
T=T.subs(dict(zip(x,next(iter(sol)))))
for axis in ([1,0,0],[0,1,0],[0,0,1]):print('diag Jordan',axis,[S.factor(z) for z in (T*G(axis)+G(axis)*T).diagonal()],flush=True)
print(T.applyfunc(S.factor))
r,h=S.symbols('r h',real=True)
for s in (1,-1):
 T=S.Matrix([[0,-1,-I*s*c,-I*s*d],[-1,0,s*d,-s*c],[I*s*c,s*d,0,I],[I*s*d,-s*c,-I,0]])
 assert all(S.cancel(z)==0 for z in s*V*T-T)
 assert list((T*G([r,h*c,h*d])+G([r,h*c,h*d])*T).diagonal())==[-2*r,-2*r,2*r,2*r]
print('both-sign exceptional extra line killed',flush=True)
M=S.linear_eq_to_matrix(eq,x)[0]
rowidx=M.T.rref()[1];colidx=M[list(rowidx),:].rref()[1]
minor=S.factor(M[list(rowidx),list(colidx)].det())
assert len(rowidx)==len(colidx)==11
print('rank-eleven minor',minor,flush=True)
# Its only real zeros are chart endpoints; t=d/(1+c) lies strictly between zero and one.
num=S.fraction(minor)[0]
assert S.factor(num)==-4*t*(t-1)**2*(t+1)**2*(t**2+1)**8
Path(__file__).with_name('global_exceptional_complex_term_verification.json').write_text(json.dumps({'supported_complex_kernel_dimension':1,'rank_minor':str(minor),'chart':'0<t=d/(1+c)<1','both_signs_extra_line_killed_by_r':True},indent=2)+'\n')
