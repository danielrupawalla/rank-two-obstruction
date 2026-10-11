"""Exact reduction of the single-zero (exceptional) stratum: cross, strong and norm equations on both branches."""
import sys,itertools,json
from pathlib import Path
import sympy as S

r,h,c,d,al,be,ga,de,k,p,q,t,l,n=S.symbols('r h c d al be ga de k p q t l n',real=True)
u=S.Matrix([h,-r*c,-r*d]);v=S.Matrix([0,-d,c]);w=S.Matrix([-r,-h*c,-h*d]);x=al*v+be*w;y=ga*u+de*w;A=p*u+q*v;B=t*u+l*v+n*w
# Cross equation from rank-one frame; real x cross A - y cross B = y-k e0, imag x cross B + y cross A = -x.
e0=S.Matrix([1,0,0])
def red(z):
 z=S.expand(z)
 for var,rel in [(d,d*d+c*c-1),(h,h*h+r*r-1)]:z=S.Poly(z,var).rem(S.Poly(rel,var)).as_expr()
 return S.factor(z)
actual=[]
for vec in [x.cross(A)-y.cross(B)-y+k*e0,x.cross(B)+y.cross(A)+x]:
 actual.extend([red(axis.dot(vec)) for axis in [u,v,w]])
expected=[-be*q+de*l-ga+k*h,be*p+ga*n-de*t,-al*p-ga*l-de-k*r,al*n-be*l-de*q,be*t+de*p+al,-al*t+ga*q+be]
assert all(red(z-t)==0 for z,t in zip(actual,expected))
strong=[red(axis.dot(A+A.cross(B)+16*x.cross(y))) for axis in [u,v,w]]
assert strong==[p+q*n+16*al*de,q-p*n+16*be*ga,p*l-q*t-16*al*ga]
print('cross and strong equations independently derived',flush=True)
# beta=0 alpha nonzero branch. delta^2=1/16, gamma=h*k/2 and alpha^2=delta^2+k*r*delta-k^2*h^2/4.
subs={be:0,ga:h*k/2,p:-al/de,q:0,t:0,l:-h*k/(2*de),n:0}
def branch1(z):
 z=S.together(z.subs(subs));num,den=z.as_numer_denom();num=S.Poly(S.expand(num),al).rem(S.Poly(al*al-de*de-k*r*de+k*k*h*h/4,al)).as_expr();num=S.Poly(S.expand(num),de).rem(S.Poly(de*de-S.Rational(1,16),de)).as_expr();return red(num)/den
norms=[A[i]**2+B[i]**2+16*(x[i]**2+y[i]**2) for i in range(3)]
b1=[S.factor(branch1(z)) for z in norms]
assert S.simplify(b1[0].subs(de,S.Rational(1,4))-1)==0 and S.simplify(b1[0].subs(de,-S.Rational(1,4))-1)==0
# Norms 1,2 equal one only if 16*delta*kappa*r=0, impossible since delta^2=1/16, r>0, kappa!=0.
assert all(S.expand(S.cancel(b*16*de**2)-(16*de*k*r+1))==0 for b in b1[1:])  # with 16*delta^2=1 each norm equals 1+16*delta*kappa*r
print('branch1 norms force 16*delta*kappa*r=0: contradiction',flush=True)
# branch 2 uniform linear solution.
D=be*be+de*de;sol={p:-(al*de+be*ga*n)/D,q:(al*de*n-be*(ga-h*k))/D,t:(de*ga*n-al*be)/D,l:(de*(ga-h*k)+al*be*n)/D}
cone=al*al+be*be-ga*ga-de*de+k*(h*ga-r*de)
def branch2(z):
 z=S.together(z.subs(sol).subs(de,-r*k/2));num,den=z.as_numer_denom();num=S.Poly(S.expand(num),al).rem(S.Poly(cone.subs(de,-r*k/2),al)).as_expr();return red(num)/den
assert all(S.factor(branch2(z))==0 for z in expected)
target=-al*(ga*(n*n+16*D+1)-h*k)/D
assert S.factor(branch2(strong[2]-target))==0
print('branch2 solution and forcing factor verified',flush=True)
# alpha=0, n=0, gamma^2=1/16 normalized cone.
sub2={al:0,n:0,de:-r*k/2,p:0,q:-be/ga,t:0,l:-r*k/(2*ga)}
def branch2norm(z):
 z=S.together(z.subs(sub2));num,den=z.as_numer_denom();num=S.Poly(S.expand(num),be).rem(S.Poly(cone.subs({al:0,de:-r*k/2}),be)).as_expr();num=S.Poly(S.expand(num),ga).rem(S.Poly(ga*ga-S.Rational(1,16),ga)).as_expr();return red(num)/den
b2=[S.factor(branch2norm(z)) for z in norms]
assert b2[0]==1
# Norms 1,2 equal one only if 16*gamma*h*kappa=0, impossible since gamma^2=1/16, h>0, kappa!=0.
assert all(S.expand(S.cancel(b*16*ga**2)-(1-16*ga*h*k))==0 for b in b2[1:])  # with 16*gamma^2=1 each norm equals 1-16*gamma*h*kappa
print('branch2 norms force 16*gamma*h*kappa=0: contradiction',flush=True)
Path(__file__).with_suffix('.json').write_text(json.dumps({'cross_equations':True,'strong_equations':True,'branch_two_uniform_solution':True,'branch_two_forcing_factor':True,'branch_contradictions_asserted':True,'scope':'Checks the algebra of the exceptional-stratum reduction; the case split is in the manuscript.'},indent=2)+'\n')
