"""Symbolic two-direction basis of supported traceless matrices on a full-support involution (used in the pair chart)."""
import sys,itertools,json
from pathlib import Path
import sympy as sp
import helpers;ns=vars(helpers)
qm,qc,mm,ma,ms,identity,qa=[ns[n] for n in ('qm','qc','mm','ma','ms','identity','qa')]
r,s,t=sp.symbols('r s t',real=True);ii=sp.I
M=sp.Matrix([[0,r,s,t],[r,0,ii*t,-ii*s],[s,-ii*t,0,ii*r],[t,ii*s,-ii*r,0]])
X=[[(sp.re(M[a,b]),sp.im(M[a,b]),0,0) for b in range(4)] for a in range(4)]
rel=r*r+s*s+t*t-1
def red(z):return sp.Poly(sp.expand(z),r,s,t).rem(sp.Poly(rel,r,s,t)).as_expr()
def equal(A,B):return all(red(x-y)==0 for a,b in zip(A,B) for q,p in zip(a,b) for x,y in zip(q,p))
reports=[]
for sign in (1,-1):
 e=(sp.eye(4)+sign*M)/2;v=e[:,0];w=e[:,1];F=v*w.T-w*v.T
 Y=[]
 for scalar in (1,ii):
  A=scalar*F;B=[[(0,0,sp.re(A[a,b]).expand(),sp.im(A[a,b]).expand()) for b in range(4)] for a in range(4)]
  E=[[(sp.re(e[a,b]),sp.im(e[a,b]),0,0) for b in range(4)] for a in range(4)]
  assert all(red(x)==0 for a in range(4) for x in B[a][a])
  assert equal(mm(X,B),ms(sign,B));assert equal(mm(B,X),ms(sign,B))
  assert equal(mm(B,B),ms((s*s+t*t)/4,E));Y.append(B)
 assert equal(ma(mm(*Y),mm(*Y[::-1])),ms(0,identity))
 reports.append({'eigenspace_sign':sign,'basis':'(e col0 wedge e col1) j and i times that matrix','squared_basis_norm_on_support':'(s^2+t^2)/4','mixed_anticommutator_zero':True})
report={'arithmetic':'symbolic rational polynomials modulo r^2+s^2+t^2=1','full_support_requires':'r,s,t nonzero','two_real_basis_directions_per_support':reports,'rank_test':'For each sign, 12x2 matrix Tr_H((P_a Q_b+Q_b P_a)Y_sign,k), a=0..3,b=0..2. L_sign nonzero iff rank <=1.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
