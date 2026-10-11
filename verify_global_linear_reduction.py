"""Exact linear reduction of the strong equations for a dephased quaternionic frame, and the DD* identity."""
import sys,itertools,json
from pathlib import Path
import sympy as S

I=S.I
p,q,t=S.symbols('p q t',real=True);A=S.Matrix(S.symbols('A:3',real=True));B=S.Matrix(S.symbols('B:3',real=True));x=S.Matrix(S.symbols('x:3',real=True));y=S.Matrix(S.symbols('y:3',real=True));z=A+I*B;xi=x+I*y
C=S.zeros(4);C[0,0]=S.Rational(1,4)
for i in range(3):C[0,i+1]=S.conjugate(z[i])/4;C[i+1,0]=z[i]/4
for i in range(3):
 for j in range(3):C[i+1,j+1]=z[i]*S.conjugate(z[j])/4+4*xi[i]*S.conjugate(xi[j])
G=S.Matrix([[0,p,q,t],[p,0,I*t,-I*q],[q,-I*t,0,I*p],[t,I*q,-I*p,0]])
actual=[S.expand(S.re(a)) for a in (C*G+G*C).diagonal()]
F=A+A.cross(B)+16*x.cross(y)
for i,a in enumerate(actual):print('Jordan row',i,S.factor(a),flush=True)
# Linear equivalence: first row gives m.A=0; other three reduce m_i F_i=0 up to common first row.
comb=[S.expand(2*actual[i]-sum(actual[1:])+actual[0]) for i in range(1,4)]
for i,a in enumerate(comb):assert S.expand(a-[p,q,t][i]*F[i])==0,(i,a,F[i])
print('frame strong equations equivalence passed',flush=True)
# Decomposable self-dual complex skew matrix.
D=S.zeros(4)
for (i,j),a in zip([(0,1),(0,2),(0,3),(2,3),(1,3),(1,2)],list(xi)+[I*xi[0],-I*xi[1],I*xi[2]]):D[i,j]=a;D[j,i]=-a
N=x.dot(x)+y.dot(y);cross=x.cross(y)
DD=D*D.conjugate().T
expected=N*S.eye(4)-2*S.Matrix([[0,cross[0],cross[1],cross[2]],[cross[0],0,I*cross[2],-I*cross[1]],[cross[1],-I*cross[2],0,I*cross[0]],[cross[2],I*cross[1],-I*cross[0],0]])
assert all(S.expand(z)==0 for z in DD-expected)
print('DD* identity passed without decomposability assumptions',flush=True)
Path(__file__).with_suffix('.json').write_text(json.dumps({'frame_strong_equivalence':True,'skew_DDstar_identity':True},indent=2)+'\n')
