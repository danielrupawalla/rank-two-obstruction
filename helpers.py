"""Shared exact quaternion-matrix helpers (standard library only).

U and W are the frame arrays of the fixed mixed strong triple used by
verify_mixed_strong_triple.py and verify_global_coarse_reduction.py.
"""
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
U=[[O,O,O,O],[O,qs(-1,O),J,qs(-1,J)],
   [O,I,tuple(F(t,2) for t in [-1,-1,-1,1]),tuple(F(t,2) for t in [-1,-1,1,-1])],
   [O,qs(-1,I),tuple(F(t,2) for t in [-1,1,-1,-1]),tuple(F(t,2) for t in [-1,1,1,1])]]
W=[[O,O,O,O],[K,K,qs(-1,K),qs(-1,K)],
   [tuple(F(t,2) for t in [-1,-1,1,1]),tuple(F(t,2) for t in [1,1,-1,-1]),tuple(F(t,2) for t in [-1,-1,-1,-1]),tuple(F(t,2) for t in [1,1,1,1])],
   [tuple(F(t,2) for t in [1,-1,1,-1]),tuple(F(t,2) for t in [-1,1,-1,1]),tuple(F(t,2) for t in [1,-1,-1,1]),tuple(F(t,2) for t in [-1,1,1,-1])]]
def mm(A,B):return [[qa(*(qm(A[i][k],B[k][j]) for k in range(4))) for j in range(4)] for i in range(4)]
def ma(A,B):return [[qa(a,b) for a,b in zip(x,y)] for x,y in zip(A,B)]
def ms(c,A):return [[qs(c,z) for z in row] for row in A]
identity=[[O if i==j else Z for j in range(4)] for i in range(4)]
def projections(M):return [[[qs(F(1,4),qm(M[a][b],qc(M[d][b]))) for d in range(4)] for a in range(4)] for b in range(4)]
