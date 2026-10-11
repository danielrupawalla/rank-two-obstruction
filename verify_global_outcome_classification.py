"""Exact all-parameter candidate and complete-setting identities for the global proof. Requires SymPy."""
import sys,itertools,json
from pathlib import Path
import sympy as S
from fractions import Fraction as F
from pathlib import Path
import itertools,json
O=(F(1),F(0),F(0),F(0));Z=(F(0),)*4;I=(0,1,0,0);J=(0,0,1,0);K=(0,0,0,1)
def qm(x,y):
 a,b,c,d=x;e,f,g,h=y
 return (a*e-b*f-c*g-d*h,a*f+b*e+c*h-d*g,a*g-b*h+c*e+d*f,a*h+b*g-c*f+d*e)
def qc(x):return (x[0],-x[1],-x[2],-x[3])
def qs(c,x):return tuple(c*a for a in x)
def qa(*xs):return tuple(sum(x[a] for x in xs) for a in range(4))
def dot(x,y):return sum(a*b for a,b in zip(x,y))
def mm(A,B):return [[qa(*(qm(A[i][k],B[k][j]) for k in range(4))) for j in range(4)] for i in range(4)]
def ma(A,B):return [[qa(a,b) for a,b in zip(x,y)] for x,y in zip(A,B)]
def ms(c,A):return [[qs(c,z) for z in row] for row in A]
identity=[[O if i==j else Z for j in range(4)] for i in range(4)]
def projections(M):return [[[qs(F(1,4),qm(M[a][b],qc(M[d][b]))) for d in range(4)] for a in range(4)] for b in range(4)]
r,h,c,d,a,b=S.symbols('r h c d a b',real=True);ii=S.I
Xc=S.Matrix([[0,r,h*c,h*d],[r,0,ii*h*d,-ii*h*c],[h*c,-ii*h*d,0,ii*r],[h*d,ii*h*c,-ii*r,0]])
X=[[(S.re(Xc[x,y]),S.im(Xc[x,y]),0,0) for y in range(4)] for x in range(4)]
rels=[(h,h*h+r*r-1),(d,d*d+c*c-1),(b,b*b+a*a-1)]
def red(z):
 for v,rel in rels:z=S.Poly(S.expand(z),v).rem(S.Poly(rel,v)).as_expr()
 return S.expand(z)
def equal(A,B):return all(red(x-y)==0 for a,b in zip(A,B) for q,p in zip(a,b) for x,y in zip(q,p))
def Gamma(v):
 p,q,t=v;return S.Matrix([[0,p,q,t],[p,0,ii*t,-ii*q],[q,-ii*t,0,ii*p],[t,ii*q,-ii*p,0]])
def fromcd(C,D):return [[tuple(S.expand(z) for z in (S.re(C[i,j]),S.im(C[i,j]),S.re(D[i,j]),S.im(D[i,j]))) for j in range(4)] for i in range(4)]
def project(q):return [[qm(q[i],qc(q[j])) for j in range(4)] for i in range(4)]
def tau(M):return red(sum(M[k][k][0] for k in range(4)))
u=S.Matrix([a*h,-a*r*c+b*d,-a*r*d-b*c]);v=S.Matrix([b*h,-b*r*c-a*d,-b*r*d+a*c]);w=S.Matrix([-r,-h*c,-h*d])
classes=[]
for prod in (1,-1):
 rows=[]
 for eps in itertools.product((1,-1),repeat=3):
  if eps[0]*eps[1]*eps[2]!=prod:continue
  rows.append(project([(S.Rational(1,2),0,0,0)]+[(0,S.Rational(e,2),0,0) for e in eps]))
 classes.append(rows)
for name,n0,x0,y0 in [('u',u,v,w),('v',v,w,u)]:
 rows=[]
 for s,e in itertools.product((1,-1),repeat=2):
  xi=e*(x0-ii*s*y0)/4;D=S.zeros(4)
  for (i,j),z in zip([(0,1),(0,2),(0,3),(2,3),(1,3),(1,2)],list(xi)+[ii*xi[0],-ii*xi[1],ii*xi[2]]):D[i,j]=z;D[j,i]=-z
  rows.append(fromcd((S.eye(4)+s*Gamma(n0))/4,D))
 classes.append(rows)
for k,rows in enumerate(classes):
 for Q in rows:assert equal(mm(Q,Q),Q)
 assert equal(ma(ma(rows[0],rows[1]),ma(rows[2],rows[3])),identity)
print('generic projectors passed',flush=True)
# All distinct-class overlaps are 1/4, within-class distinct are zero.
for i,j in itertools.combinations(range(16),2):
 Q=classes[i//4][i%4];T=classes[j//4][j%4];assert tau(mm(Q,T))==(0 if i//4==j//4 else S.Rational(1,4)),(i,j,tau(mm(Q,T)))
print('generic orthogonality exhaustion passed',flush=True)
def fixed_R(phasea,phaseb):
 out=[]
 for sign in (1,-1):
  F=S.zeros(4);phase=phasea+sign*ii*phaseb
  for (i,j),z in {(0,1):h/2,(0,2):(-r*c-sign*ii*d)/2,(0,3):(-r*d+sign*ii*c)/2,(1,2):(-sign*c-ii*r*d)/2,(1,3):(-sign*d+ii*r*c)/2,(2,3):ii*h/2}.items():F[i,j]=phase*z;F[j,i]=-phase*z
  for eta in (1,-1):out.append(fromcd((S.eye(4)+sign*Gamma([r,h*c,h*d]))/4,eta*F/2))
 return out

for rows in classes:
 for Q in rows:
  for T in fixed_R(a,b):
   for k0 in range(4):
    assert red(sum(qm(T[k0][j],Q[j][k0])[0] for j in range(4)))==S.Rational(1,16)
print('all generic strong moment identities passed',flush=True)
classes=[[[[tuple(z.subs({a:1,b:0}) for z in q) for q in row] for row in Q] for Q in rows] for rows in classes]
rows=[]
for e,eta in itertools.product((1,-1),repeat=2):
 rows.append(project([(S.Rational(1,2),0,0,0),(0,0,0,S.Rational(e,2)),(0,-eta*d/2,0,eta*e*c/2),(0,eta*c/2,0,eta*e*d/2)]))
classes.append(rows)
for Q in rows:assert equal(mm(Q,Q),Q)
for i,j in itertools.combinations(range(20),2):
 value=tau(mm(classes[i//4][i%4],classes[j//4][j%4]))
 
 if i//4==j//4:assert value==0,(i,j,value)
 elif j>=16:
  allowed=[(1+sgn*z)/4 for sgn in (1,-1) for z in (r,h,c,d)]
  assert value in allowed,(i,j,value)
 else:assert value==S.Rational(1,4)
print('exceptional five K4 classes passed',flush=True)

# Independently check all candidate outcomes against the defining strong-triple linear moments.
def fixed_R(phasea,phaseb):
 out=[]
 for sign in (1,-1):
  F=S.zeros(4);phase=phasea+sign*ii*phaseb
  for (i,j),z in {(0,1):h/2,(0,2):(-r*c-sign*ii*d)/2,(0,3):(-r*d+sign*ii*c)/2,(1,2):(-sign*c-ii*r*d)/2,(1,3):(-sign*d+ii*r*c)/2,(2,3):ii*h/2}.items():F[i,j]=phase*z;F[j,i]=-phase*z
  for eta in (1,-1):out.append(fromcd((S.eye(4)+sign*Gamma([r,h*c,h*d]))/4,eta*F/2))
 return out
R=fixed_R(1,0)
for rows in classes:
 for Q in rows:
  for T in R:
   for k0 in range(4):
    moment=red(sum(qm(T[k0][j],Q[j][k0])[0] for j in range(4)))
    assert moment==S.Rational(1,16),(k0,moment)
print('all exceptional strong moment identities passed',flush=True)

for axis,sgn in [(u,-1),(v,1)]:
 T=fromcd(Gamma(axis),S.zeros(4));R=fixed_R(a,b)
 Ap=ma(R[0],ms(-1,R[1]));Am=ma(R[2],ms(-1,R[3]))
 assert equal(mm(mm(T,Ap),T),ms(sgn,Am))
 assert equal(ma(mm(T,X),mm(X,T)),ms(0,identity))
print('all-parameter matching intertwiners passed',flush=True)
# Canonical blocks of the exceptional frame lie in one complex subfield.
U=[[tuple(4*z for z in classes[4][b][a0][0]) for b in range(4)] for a0 in range(4)]
canonical=[[qm(qc(U[a0][0]),U[a0][b]) for b in range(4)] for a0 in range(4)]
assert all(red(q[1])==red(q[3])==0 for row in canonical for q in row)
for q1,q2 in itertools.combinations(sum(canonical,[]),2):
 assert all(red(x-y)==0 for x,y in zip(qm(q1,q2),qm(q2,q1)))
print('exceptional canonical commutativity passed',flush=True)
report={'arithmetic':'Exact polynomial quaternion arithmetic modulo circle identities',
 'generic_candidate_outcomes':16,'exceptional_candidate_outcomes':20,
 'all_candidate_projectors_and_strong_moments':True,
 'generic_orthogonality_graph':'four disjoint K4',
 'exceptional_orthogonality_graph':'five disjoint K4; cross overlaps strictly positive for r,h,c,d in (0,1)',
 'matching_intertwiners':True,
 'limits':'Exhaustion and chart coverage are analytic arguments in the paper, not consequences of candidate enumeration.'}
Path(__file__).with_name('global_outcome_classification_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
