"""Exact check of the supported-kernel description of complete third settings (Proposition on supported kernels) on a fixed exact triple."""
from fractions import Fraction as F
from pathlib import Path
import itertools,json
import helpers;ns=vars(helpers)
qm,qc,qa,qs,mm,ma,ms,O,Z,U,W,identity=[ns[n] for n in ('qm','qc','qa','qs','mm','ma','ms','O','Z','U','W','identity')]
def tr(A):return sum(A[a][a][0] for a in range(4))
def sub(A,B):return ma(A,ms(-1,B))
def flat(A):return list(itertools.chain.from_iterable(A[a][b] for a,b in itertools.product(range(4),repeat=2)))
def lin(cs,B):
 A=ms(0,identity)
 for c,b in zip(cs,B):A=ma(A,ms(c,b))
 return A

def rref(A):
 A=[list(map(F,row)) for row in A];pivs=[];r=0
 for c in range(len(A[0])):
  p=next((i for i in range(r,len(A)) if A[i][c]),None)
  if p is None:continue
  A[r],A[p]=A[p],A[r];v=A[r][c];A[r]=[x/v for x in A[r]]
  for i in range(len(A)):
   if i!=r and A[i][c]:v=A[i][c];A[i]=[x-v*y for x,y in zip(A[i],A[r])]
  pivs.append(c);r+=1
  if r==len(A):break
 return A,pivs

def kernel(A):
 R,p=rref(A);out=[]
 for c in range(len(A[0])):
  if c in p:continue
  x=[F(0)]*len(A[0]);x[c]=F(1)
  for row,j in zip(R,p):x[j]=-row[c]
  out.append(x)
 return out
basis=[]
for a in range(4):
 B=ms(0,identity);B[a][a]=O;basis.append(B)
for a,b in itertools.combinations(range(4),2):
 for q in (O,ns['I'],ns['J'],ns['K']):
  B=ms(0,identity);B[a][b]=q;B[b][a]=qc(q);basis.append(B)
P=[[[O if a==b==c else Z for b in range(4)] for a in range(4)] for c in range(4)]
Q=ns['projections'](U);R=ns['projections'](W)
A=[[tr(mm(ma(mm(p,q),mm(q,p)),b)) for b in basis] for p,q in itertools.product(P,Q)]
kb=kernel(A);K=[lin(c,basis) for c in kb]
def inK(B):return all(tr(mm(ma(mm(p,q),mm(q,p)),B))==0 for p,q in itertools.product(P,Q))
assert len(K)==13
rows=[]
for signs in ((1,1,-1,-1),(1,-1,1,-1),(1,-1,-1,1)):
 X=lin(signs,R);assert inK(X) and mm(X,X)==identity and tr(X)==0
 dims=[];components=[]
 for sign in (1,-1):
  C=[flat(sub(mm(X,B),ms(sign,B))) for B in K]
  L=[lin(c,K) for c in kernel(list(map(list,zip(*C))))];dims.append(len(L))
  assert L
  for Y in L:
   assert inK(Y) and mm(X,Y)==ms(sign,Y) and mm(Y,X)==ms(sign,Y)
   assert tr(Y)==0
   e=ms(F(1,2),ma(identity,ms(sign,X)))
   assert mm(Y,Y)==ms(tr(mm(Y,Y))/2,e)
  # Choose a physical difference of the two outcomes in this eigenspace.
  inds=[i for i,x in enumerate(signs) if x==sign];Y=sub(R[inds[0]],R[inds[1]])
  assert inK(Y) and mm(X,Y)==ms(sign,Y)
  assert mm(Y,Y)==e
  components.append(Y)
 Y=ma(*components);assert mm(Y,Y)==identity and mm(X,Y)==mm(Y,X) and inK(mm(X,Y))
 rebuilt=[ms(F(1,4),ma(ma(identity,ms(s,X)),ma(ms(t,Y),ms(s*t,mm(X,Y))))) for s,t in itertools.product((1,-1),repeat=2)]
 assert all(tr(T)==1 and mm(T,T)==T for T in rebuilt)
 for T in rebuilt:
  assert inK(sub(T,ms(F(1,4),identity)))
 for a,b in itertools.product(range(4),repeat=2):assert mm(rebuilt[a],rebuilt[b])==(rebuilt[a] if a==b else ms(0,identity))
 assert all(any(T==r for r in R) for T in rebuilt)
 rows.append({'Walsh_signs':signs,'supported_kernel_dimensions':dims,'reconstructed_complete_setting':True})
report={'arithmetic':'exact rational quaternion matrices; standard library only','fixed_fixture':'special i,j mixed strong completion','first_pair_jordan_span_rank':28-len(K),'kernel_dimension':len(K),'balanced_involutions':rows,'scope':'Exact implementation checks of a parameter-free analytic equivalence, not global commuting-pair necessity.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
