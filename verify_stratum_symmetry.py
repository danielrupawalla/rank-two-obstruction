"""Exact admissible chart symmetries and quaternion-embedding convention change."""
import sys,itertools,json
from pathlib import Path
import sympy as S
m=S.symbols('m0:3',real=True);u=S.symbols('u0:3',real=True);v=S.symbols('v0:3',real=True)
zero=(0,0,0,0)
def conjugate(q):return q[0],-q[1],-q[2],-q[3]
def matrix(entries):
 A=[[zero for j in range(4)] for i in range(4)]
 for (i,j),q in entries.items():A[i][j]=q;A[j][i]=conjugate(q)
 return A
def Gamma(m):
 p,q,t=m
 return matrix({(0,1):(p,0,0,0),(0,2):(q,0,0,0),(0,3):(t,0,0,0),(1,2):(0,t,0,0),(1,3):(0,-q,0,0),(2,3):(0,p,0,0)})
def coarse(u,v):
 B=matrix({(0,1):(0,0,u[0],0),(0,2):(0,0,u[1],0),(0,3):(0,0,u[2],0),(2,3):(0,0,0,u[0]),(1,3):(0,0,0,-u[1]),(1,2):(0,0,0,u[2])})
 C=matrix({(0,1):(0,0,0,v[0]),(0,2):(0,0,0,v[1]),(0,3):(0,0,0,v[2]),(2,3):(0,0,-v[0],0),(1,3):(0,0,v[1],0),(1,2):(0,0,-v[2],0)})
 return B,C
def mapped(A,f):return [[f(q) for q in row] for row in A]
def equal(A,B):return all(S.expand(x-y)==0 for row,other in zip(A,B) for q,p in zip(row,other) for x,y in zip(q,p))
X=Gamma(m);B,C=coarse(u,v)
for perm in itertools.permutations(range(3)):
 sign=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3));indices=[0]+[i+1 for i in perm]
 def transform(A):
  A=[[A[i][j] for j in indices] for i in indices]
  return mapped(A,lambda q:(q[0],sign*q[1],q[2],sign*q[3]))
 mp=[m[i] for i in perm];up=[u[i] for i in perm];vp=[sign*v[i] for i in perm];Bp,Cp=coarse(up,vp)
 assert equal(transform(X),Gamma(mp));assert equal(transform(B),Bp);assert equal(transform(C),Cp)
rot=lambda q:(q[0],q[1],q[3],-q[2])
Bp,Cp=coarse(v,[-z for z in u]);assert equal(mapped(C,rot),Bp);assert equal(mapped(B,rot),Cp);assert equal(mapped(X,rot),X)
# Explicit internal-qubit convention change.
I=S.I;paulix=S.Matrix([[0,1],[1,0]]);pauliy=S.Matrix([[0,-I],[I,0]]);pauliz=S.diag(1,-1);H=S.Matrix([[1,1],[1,-1]])/S.sqrt(2)
for old,new in [(I*paulix,I*pauliz),(-I*pauliy,I*pauliy),(I*pauliz,I*paulix)]:assert H*old*H==new
report={'permutations_checked':6,'coarse_support_transformations':'u -> perm(u), v -> sign(perm) perm(v), m -> perm(m)','coarse_exchange_rotation':'(u,v) -> (v,-u); X fixed','quaternion_embedding_hadamard_conjugation':True,'arithmetic':'exact symbolic; no parameter sampling'}
Path(__file__).with_name('stratum_symmetry_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
