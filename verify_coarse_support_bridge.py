"""Symbolic identities used in the coarse-support lemma (matching, four-cycle and full-support cases)."""
import sys,itertools,json
from pathlib import Path
import sympy as sp
import helpers;ns=vars(helpers)
qm,qc,qa,qs,mm,ma,ms,O,Z,identity=[ns[n] for n in ('qm','qc','qa','qs','mm','ma','ms','O','Z','identity')]
def simplify(A):return [[tuple(sp.factor(z) for z in q) for q in row] for row in A]
def equal(A,B):return all(sp.factor(x-y)==0 for a,b in zip(A,B) for q,p in zip(a,b) for x,y in zip(q,p))
t=sp.symbols('t',real=True);r=(1-t*t)/(1+t*t);s=2*t/(1+t*t)
a1,a2,a3,f1,f2,f3=sp.symbols('a1 a2 a3 f1 f2 f3',real=True)
a=(0,a1,a2,a3);f=(0,f1,f2,f3)
X=[[Z,Z,qs(r,O),qs(s,O)],[Z,Z,qs(s,O),qs(-r,O)],[qs(r,O),qs(s,O),Z,Z],[qs(s,O),qs(-r,O),Z,Z]]
Y=[[Z,a,qs(s,f),qs(-r,f)],[qc(a),Z,qs(-r,f),qs(-s,f)],[qc(qs(s,f)),qc(qs(-r,f)),Z,qs(-1,a)],[qc(qs(-r,f)),qc(qs(-s,f)),a,Z]]
Zxy=simplify(mm(X,Y));assert equal(mm(X,X),identity);assert equal(mm(X,Y),mm(Y,X));assert all(sp.factor(x)==0 for k in range(4) for x in Zxy[k][k])
norm=a1*a1+a2*a2+a3*a3+f1*f1+f2*f2+f3*f3;dot=a1*f1+a2*f2+a3*f3
D=simplify(ma(mm(Y,Y),ms(-norm,identity)))
# Every off-diagonal square defect is a real multiple of dot(a,f).
for row in D:
 for q in row:
  for v in q:
   if v:assert sp.factor(v/dot).free_symbols <= {t}
# Full-support involution canonical form and compressed Gram eigenvalues.
u1,u2,u3=sp.symbols('u1 u2 u3',real=True);i=ns['I']
F=[[Z,qs(u1,O),qs(u2,O),qs(u3,O)],
 [qs(u1,O),Z,qs(u3,i),qs(-u2,i)],
 [qs(u2,O),qs(-u3,i),Z,qs(u1,i)],
 [qs(u3,O),qs(u2,i),qs(-u1,i),Z]]
assert equal(mm(F,F),ms(u1*u1+u2*u2+u3*u3,identity))
G=sp.Matrix(4,4,lambda a,b:sp.Rational(1,8) if a==b else sum(v*v for v in F[a][b])/4-sp.Rational(1,8))
for signs,u in [((1,1,-1,-1),u1),((1,-1,1,-1),u2),((1,-1,-1,1),u3)]:
 v=sp.Matrix(signs);err=G*v-u*u/2*v
 assert all(sp.factor(z).subs(u3*u3,1-u1*u1-u2*u2)==0 for z in err)
report={'symbolic_C4_commuting_zero_diagonal_form_checked':True,'C4_remaining_square_conditions':['norm(a)^2+norm(f)^2=1','dot(a,f)=0'],'if_second_coarse_support_not_full':'a=0 or f=0; one of Y or XY has matching support','full_support_normal_form_square_checked':True,'compressed_P_Bloch_Gram_eigenvalues':['0','u1^2/2','u2^2/2','u3^2/2'],'scope':'Checks the matrix identities used in the coarse-support lemma. The case analysis itself, and its use in the global commuting-pair theorem, is in the manuscript.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
