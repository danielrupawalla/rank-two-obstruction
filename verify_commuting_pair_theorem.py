"""Exact polynomial identities supporting the uniform commuting-pair theorem.

Analytic exhaustiveness and normal-form arguments are in the manuscript.
"""
import itertools,json
from pathlib import Path
from fractions import Fraction as F
def add(A,B):
 C=A.copy()
 for k,v in B.items():C[k]=C.get(k,0)+v
 return {k:v for k,v in C.items() if v}
def sc(x,A):return {k:x*v for k,v in A.items() if x*v}
def mul(A,B):
 C={}
 for a,x in A.items():
  for b,y in B.items():e=tuple(i+j for i,j in zip(a,b));C[e]=C.get(e,0)+x*y
 return {k:v for k,v in C.items() if v}
one={(0,0,0,0):F(1)};var=[]
for i in range(4):e=[0]*4;e[i]=1;var.append({tuple(e):F(1)})
c,s,d,e=var
def reduce(A):
 A=A.copy()
 for axis,other in [(1,0),(3,2)]:
  while any(k[axis]>=2 for k in A):
   C={}
   for k,v in A.items():
    if k[axis]<2:C=add(C,{k:v});continue
    f=list(k);f[axis]-=2;C=add(C,{tuple(f):v});f[other]+=2;C=add(C,{tuple(f):-v})
   A=C
 return A
plus=lambda x:sc(F(1,2),add(one,x))
minus=lambda x:sc(F(1,2),add(one,sc(-1,x)))
single={('X','Y'):('Z',1),('Y','Z'):('X',1),('Z','X'):('Y',1),('Y','X'):('Z',3),('Z','Y'):('X',3),('X','Z'):('Y',3)}
def pm(a,b):
 word='';phase=0
 for x,y in zip(a,b):
  if x=='I':word+=y
  elif y=='I':word+=x
  elif x==y:word+='I'
  else:z,p=single[x,y];word+=z;phase+=p
 return word,phase%4
def sym(a,b,c):
 total=0;label=None
 for x,y,z in itertools.permutations([a,b,c]):
  t,p=pm(x,y);t,q=pm(t,z);phase=(p+q)%4
  label=t;total +=[1,1j,-1,-1j][phase]
 assert complex(total).imag==0
 return label,int(complex(total).real)
A=[{'III':one},{'ZII':one},{'IZI':one},{'ZZI':one}]
B=[{'III':one},{'XII':one},{'IXI':plus(c),'XXI':minus(c),'IYZ':sc(F(1,2),s),'XYZ':sc(-F(1,2),s)},{'IXI':minus(c),'XXI':plus(c),'IYZ':sc(-F(1,2),s),'XYZ':sc(F(1,2),s)}]
C=[{'III':one},{'YIZ':one},{'ZYY':plus(d),'XYX':minus(d),'XYY':sc(F(1,2),e),'ZYX':sc(F(1,2),e)},{'ZYY':minus(d),'XYX':plus(d),'XYY':sc(-F(1,2),e),'ZYX':sc(-F(1,2),e)}]
rows=[]
for ai,bi,ci in itertools.product(range(4),repeat=3):
 row={}
 for a,b,z in itertools.product(A[ai],B[bi],C[ci]):
  w,t=sym(a,b,z);p=sc(t,mul(mul(A[ai][a],B[bi][b]),C[ci][z]));row[w]=add(row.get(w,{}),p)
 row={w:p for w,p in row.items() if p};rows.append(row)
important=['YIX','YIY','YZX','YZY']
G=[[{} for j in range(4)] for i in range(4)]
for i,x in enumerate(important):
 for j,y in enumerate(important):
  out={}
  for row in rows:out=add(out,mul(row.get(x,{}),row.get(y,{})))
  G[i][j]=reduce(out)
W=[[s,mul(s,d),{},{}],[{},sc(-1,mul(s,e)),{},{}],[{},sc(-1,mul(c,e)),{},e],[sc(-1,c),sc(-1,mul(c,d)),one,d]]
for i,j in itertools.product(range(4),repeat=2):
 value={}
 for k in range(4):value=add(value,mul(W[i][k],W[j][k]))
 target=add(sc(8,value),sc(4,one) if i==j else {})
 assert reduce(add(G[i][j],sc(-1,target)))=={}
comm=['III','IZI','XII','XZI','YIX','YIY','YIZ','YZX','YZY','YZZ','ZII','ZZI']
othercomm=[x for x in comm if x not in important]
alllabels=[''.join(w) for w in itertools.product('IXYZ',repeat=3)]
quatlabels=[w for w in alllabels if ((w[0]=='Y')+(w[1]=='Y'))%2==(w[2]!='I')]
anti=[w for w in quatlabels if w not in comm]
for row in rows:
 assert not (any(row.get(x) for x in comm) and any(row.get(x) for x in anti))
for x,y in itertools.product(comm,repeat=2):
 if x in important and y in important:continue
 value={}
 for row in rows:value=add(value,mul(row.get(x,{}),row.get(y,{})))
 assert reduce(value)==(sc(36,one) if x==y else {}),(x,y)
def const(x):return sc(F(x),one)
qzero=({},)*4
qone=(one,{},{},{})
qi=({},one,{},{})
qj=({},{},one,{})
qk=({},{},{},one)
def qa(a,b):return tuple(add(x,y) for x,y in zip(a,b))
def qs(a,t):return tuple(mul(x,t) for x in a)
def qneg(a):return tuple(sc(-1,x) for x in a)
def qc(a):return a[0],sc(-1,a[1]),sc(-1,a[2]),sc(-1,a[3])
def qm(a,b):
 out=[]
 for terms in [[(0,0,1),(1,1,-1),(2,2,-1),(3,3,-1)],[(0,1,1),(1,0,1),(2,3,1),(3,2,-1)],[(0,2,1),(2,0,1),(3,1,1),(1,3,-1)],[(0,3,1),(3,0,1),(1,2,1),(2,1,-1)]]:
  x={}
  for i,j,t in terms:x=add(x,sc(t,mul(a[i],b[j])))
  out.append(x)
 return tuple(out)
def qnorm(a):return qm(a,qc(a))[0]
def qeq(a,b):return all(reduce(add(x,sc(-1,y)))=={} for x,y in zip(a,b))
def outcome(kind,tags,angle):
 x,y=angle;v=({}, {},x,y);h=({}, {},sc(-1,y),x)
 if kind=='A':
  eps,delta=tags;return [qone,v,qs(qi,const(eps)),qs(h,const(delta))]
 eta,eps=tags;a=sc(F(eta,2),s)
 u=qa(qs(qi,sc(F(1,2),add(const(eps),sc(-eta,c)))),qa(qs(h,sc(F(1,2),add(const(eps),sc(eta,c)))),qs(v,sc(-F(eta,2),s))))
 v3=qa(qs(qi,sc(F(1,2),add(const(eps),sc(eta,c)))),qa(qs(h,sc(F(1,2),add(const(eps),sc(-eta,c)))),qs(v,sc(F(eta,2),s))))
 return [qone,qa((a,{},{},{}),u),v,qa((sc(-1,a),{},{},{}),v3)]
def inner(U,V):
 z=qzero
 for a,b in zip(U,V):z=qa(z,qm(qc(a),b))
 return z
q=(c,s,{},{});U=[[qone,qone,qone,qone],[qone,qneg(qone),q,qneg(q)],[qone,qone,qneg(qone),qneg(qone)],[qone,qneg(qone),qneg(q),q]]
signs=list(itertools.product([1,-1],repeat=2))
checks=0
for kind,tags in itertools.product(['A','B'],signs):
 T=outcome(kind,tags,(d,e));S=[]
 for b in range(4):
  z=qzero
  for a in range(4):z=qa(z,qm(qc(U[a][b]),T[a]))
  S.append(qs(z,const(F(1,2))))
 for t in T+S:assert reduce(qnorm(t))==one
 for a,b in itertools.product(range(4),repeat=2):assert reduce(qm(qm(U[a][b],S[b]),qc(T[a]))[0])==const(F(1,2));checks+=1
 if kind=='B':
  transformed=[qm(S[b],qc(S[0])) for b in [0,2,1,3]]
  assert reduce(transformed[1][0])=={}
  assert qeq(transformed[2],qs(qi,const(-tags[1])))
for kind,a,b in itertools.product(['A','B'],signs,signs):
 ov=inner(outcome(kind,a,(one,{})),outcome(kind,b,(d,e)))
 expected=const(0) if a==tuple(-x for x in b) else sc(8,add(one,d)) if a==b else const(4)
 assert reduce(add(qnorm(ov),sc(-1,expected)))=={}
for a,b in itertools.product(signs,repeat=2):
 ov=inner(outcome('A',a,(one,{})),outcome('B',b,(d,e)));p=reduce(qnorm(ov))
 L={};M={}
 for power,value in p.items():
  if power[2:]==(1,0):L[power[:2]+(0,0)]=value
  elif power[2:]==(0,1):M[power[:2]+(0,0)]=value
  else:assert power==(0,0,0,0) and value==4,(a,b,p)
 parity=a[0]*a[1]*b[0]*b[1]
 assert reduce(add(add(mul(L,L),mul(M,M)),sc(-8,add(one,sc(parity,c)))))=={}
# Check observable formulas directly from quaternion outcome frames.
def qpauli(label):
 P={'I':[[1,0],[0,1]],'X':[[0,1],[1,0]],'Y':[[0,-1j],[1j,0]],'Z':[[1,0],[0,-1]]};M=[[1]]
 for t in label:M=[[a*b for a in row for b in r] for row in M for r in P[t]]
 out=[]
 for a in range(4):
  row=[]
  for b in range(4):
   z=[[M[2*a+i][2*b+j] for j in range(2)] for i in range(2)]
   vals=[(complex(z[0][0]).real+complex(z[1][1]).real)/2,(complex(z[0][0]).imag-complex(z[1][1]).imag)/2,(complex(z[0][1]).real-complex(z[1][0]).real)/2,(complex(z[0][1]).imag+complex(z[1][0]).imag)/2]
   row.append(tuple(const(F(v)) for v in vals))
  out.append(row)
 return out
H=[[1,1,1,1],[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]]
Rframes=[outcome('A',tags,angle) for tags,angle in [((1,1),(one,{})),((1,1),(sc(-1,one),{})),((-1,-1),(d,e)),((-1,-1),(sc(-1,d),sc(-1,e)))]]
for index,frames,target in [(b,[ [U[a][x] for a in range(4)] for x in range(4)],B[b]) for b in range(4)]+[(b,Rframes,C[b]) for b in range(4)]:
 for a,v in itertools.product(range(4),repeat=2):
  lhs=qzero
  for f,sign in zip(frames,H[index]):lhs=qa(lhs,qs(qm(f[a],qc(f[v])),const(F(sign,4))))
  rhs=qzero
  for label,coef in target.items():rhs=qa(rhs,qs(qpauli(label)[a][v],coef))
  assert qeq(lhs,rhs),(index,a,v)
report={'arithmetic':'Rational polynomial coefficients modulo c^2+s^2=1 and d^2+e^2=1',
        'outcome_classes_verified':8,'triangle_identities_checked':checks,
        'within_class_overlap_norms':'0 for both sign reversals, 8(1+cos angle) for identical signs, 4 for one sign reversal',
        'cross_class_overlap_norm':'4+L cos angle+M sin angle, L^2+M^2=8(1+sign_product c)',
        'second_class_exchange_verified':True,'observable_formulas_verified':True,
        'grading':'IZI','commuting_sector_dimension':12,'gram_blocks':'36 I_8 and 4 I_4 + 8 W W^T','polynomial_gram_factorization_verified':True,
        'conclusion':'Uniform positive commuting-sector Gram forces all fourth-candidate deviations to anticommute with IZI in the normalized first class. Cubic trace then excludes rank-two projector outcomes. Exhaustiveness and normalization steps are proved analytically in the manuscript.',
        'scope':'Checks the identities of the commuting-pair nonextension theorem. That every strong triple contains a commuting pair is the global commuting-pair theorem, checked by the verify_global_* scripts.'}
Path(__file__).with_name('commuting_pair_theorem_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
