"""Symbolic projector chart covering ALL noncommuting quaternionic order-four pairs.
Exhaustion follows from the analytic coarse-support proof, not symbolic examples.
"""
import sys,itertools,json
from pathlib import Path
import sympy as sp
import helpers;ns=vars(helpers)
qm,qc,mm,ma,ms,qa,identity,O,Z=[ns[n] for n in ('qm','qc','mm','ma','ms','qa','identity','O','Z')]
r,h,c,d,a,b=sp.symbols('r h c d a b',real=True);ii=sp.I
Xc=sp.Matrix([[0,r,h*c,h*d],[r,0,ii*h*d,-ii*h*c],[h*c,-ii*h*d,0,ii*r],[h*d,ii*h*c,-ii*r,0]])
X=[[(sp.re(Xc[x,y]),sp.im(Xc[x,y]),0,0) for y in range(4)] for x in range(4)]
rels=[(h,h*h+r*r-1),(d,d*d+c*c-1),(b,b*b+a*a-1)]
def red(z):
 for v,rel in rels:z=sp.Poly(sp.expand(z),v).rem(sp.Poly(rel,v)).as_expr()
 return sp.expand(z)
def equal(A,B):return all(red(x-y)==0 for a,b in zip(A,B) for q,p in zip(a,b) for x,y in zip(q,p))
projections=[]
for sign in (1,-1):
 F=sp.zeros(4)
 for (x,y),v in {(0,1):h/2,(0,2):(-sign*ii*d-r*c)/2,(0,3):(sign*ii*c-r*d)/2,(1,2):(-ii*r*d-sign*c)/2,(1,3):(ii*r*c-sign*d)/2,(2,3):ii*h/2}.items():F[x,y]=v;F[y,x]=-v
 F=(a+sign*ii*b)*F
 Y=[[(0,0,sp.re(F[x,y]).expand(),sp.im(F[x,y]).expand()) for y in range(4)] for x in range(4)]
 e=ms(sp.Rational(1,2),ma(identity,ms(sign,X)))
 assert equal(mm(Y,Y),e) and equal(mm(X,Y),ms(sign,Y))
 for sgn in (1,-1):
  Q=ms(sp.Rational(1,2),ma(e,ms(sgn,Y)))
  assert equal(mm(Q,Q),Q)
  assert all(Q[k][k]==(sp.Rational(1,4),0,0,0) for k in range(4))
  assert sum(Q[k][k][0] for k in range(4))==1
  assert all(qc(Q[y][x])==Q[x][y] for x,y in itertools.product(range(4),repeat=2))
  projections.append(Q)
assert equal(mm(X,X),identity)
assert equal(ma(ma(*projections[:2]),ma(*projections[2:])),identity)
for x,y in itertools.combinations(range(4),2):assert equal(mm(projections[x],projections[y]),ms(0,identity))
P=[[[O if x==y==k else Z for y in range(4)] for x in range(4)] for k in range(4)]
for p,q in itertools.product(P,projections):
 assert equal(mm(mm(p,q),p),ms(sp.Rational(1,4),p))
 assert equal(mm(mm(q,p),q),ms(sp.Rational(1,4),q))
report={'arithmetic':'exact symbolic rational polynomials modulo three circle identities','parameters':['r,h,c,d,a,b'],'identities':['r^2+h^2=1','c^2+d^2=1','a^2+b^2=1'],'positive_chart':'r,h,c,d>0','degrees_of_freedom':3,'outcome_projectors_polynomial':True,'PVM_relations_and_all_pair_MUM_identities':True,'coverage':'Analytic coarse-support dichotomy plus supported two-plane exhausts every noncommuting canonical pair. No special-family exclusions or parameter rank assumptions.','scope':'Chart for first pairs. Third outcomes over this chart are checked in verify_global_outcome_classification.py.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
