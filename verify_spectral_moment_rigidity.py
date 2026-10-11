"""Exact fourth/sixth averaged-moment identities; standard library only.
Checks all word types and symbolic Newton identities at arbitrary k,n.
"""
from fractions import Fraction as F
from itertools import product
from collections import Counter
from math import factorial
from pathlib import Path
import json

# Bivariate polynomials in (k,h), with h=1/n.
def add(*polys):
    out={}
    for p in polys:
        for mon,c in p.items(): out[mon]=out.get(mon,F(0))+c
    return {m:c for m,c in out.items() if c}
def scale(a,p): return {m:a*c for m,c in p.items() if a*c}
def mul(p,q):
    out={}
    for (a,b),x in p.items():
        for (c,d),y in q.items(): out[(a+c,b+d)]=out.get((a+c,b+d),F(0))+x*y
    return {m:c for m,c in out.items() if c}
one={(0,0):F(1)};h={(0,1):F(1)}
def choose_k(j):
    p=one
    for i in range(j): p=mul(p,{(1,0):F(1),(0,0):F(-i)})
    return scale(F(1,factorial(j)),p)
e=[one]
for j in range(1,7):
    p=choose_k(j)
    for i in range(j): p=mul(p,{(0,0):F(1),(0,1):F(-i)})
    e.append(scale((-1)**j,p))
power=[{}]
for m in range(1,7):
    power.append(scale(-1,add(scale(m,e[m]),*(mul(e[i],power[m-i]) for i in range(1,m)))))

def canonical(w): return min(tuple(w[i:]+w[:i]) for i in range(len(w))) if w else ()
def reduce_word(word,mum):
    w=list(word); exponent=0
    while w:
        changed=False
        if len(w)>1:
            for i in range(len(w)):
                if w[i]==w[(i+1)%len(w)]:
                    w=w[i:]+w[:i];w.pop(0);changed=True;break
        if changed:continue
        counts=Counter(w)
        for i,x in enumerate(w):
            if counts[x]==1:
                w.pop(i);exponent+=1;changed=True;break
        if changed:continue
        if mum and len(w)>2:
            for i in range(len(w)):
                z=w[i:]+w[:i]
                if z[0]==z[2]:
                    w=[z[0]]+z[3:];exponent+=1;changed=True;break
        if not changed:break
    return canonical(w),exponent

expected={
4:{1:{1:1},2:{2:12},3:{3:36},4:{4:24}},
6:{1:{1:1},2:{2:30,3:30,4:2},3:{3:300,4:234},4:{4:1200,5:360},5:{5:1800},6:{6:720}}}
tables={}; counts={}
for m in (4,6):
    baseline={};rows=[];total=0
    for ell in range(1,m+1):
        const=Counter();residual=Counter()
        for w in product(range(ell),repeat=m):
            if len(set(w))!=ell:continue
            total+=1;z,exponent=reduce_word(w,m==6)
            if not z:const[exponent]+=1
            else:residual[(z,exponent)]+=1
        assert dict(const)==expected[m].get(ell,{})
        local={(0,a):F(c) for a,c in const.items()}
        if residual:
            assert ell==(2 if m==4 else 3)
            for (z,exponent),count in residual.items():
                assert exponent==0 and count>0
                assert len(z)==m and z[:ell]==z[ell:] and len(set(z[:ell]))==ell
            assert sum(residual.values())==(2 if m==4 else 6)
            local=add(local,{(0,3):F(2)} if m==4 else {(0,5):F(12),(0,4):F(-6)})
        baseline=add(baseline,mul(choose_k(ell),local))
        rows.append({'distinct_settings':ell,'constant_trace_terms_by_power_of_1_over_n':dict(const),'irreducible_cycle_words':sum(residual.values())})
    assert baseline==mul(h,power[m]),('symbolic Newton mismatch',m)
    tables[str(m)]=rows;counts[str(m)]=total

# Independent matrix tests, including strict-positive defects.
# Gaussian rational arithmetic uses pairs (real,imag).
Z=(F(0),F(0));O=(F(1),F(0));I=(F(0),F(1))
def ca(*xs):return (sum(x[0] for x in xs),sum(x[1] for x in xs))
def cm(x,y):return (x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def cc(x):return (x[0],-x[1])
def cs(a,x):return (a*x[0],a*x[1])
def mm(A,B):return [[ca(*(cm(A[i][l],B[l][j]) for l in range(4))) for j in range(4)] for i in range(4)]
def ma(*As):return [[ca(*(A[i][j] for A in As)) for j in range(4)] for i in range(4)]
def ms(a,A):return [[cs(a,x) for x in row] for row in A]
def tr(A):
    t=ca(*(A[i][i] for i in range(4)));assert t[1]==0;return t[0]
Id=[[O if i==j else Z for j in range(4)] for i in range(4)]
def mp(A,m):
    B=Id
    for _ in range(m):B=mm(B,A)
    return B
def projections(U):return [[[cm(U[i][a],cc(U[j][a])) for j in range(4)] for i in range(4)] for a in range(4)]
H=[[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]]
U=[[cs(F(x,2),O) for x in row] for row in H]
V=[[cm(a,x) for x in row] for a,row in zip([O,I,I,cs(-1,O)],U)]
P=[[[O if i==j==a else Z for j in range(4)] for i in range(4)] for a in range(4)]
Q=projections(U);R=projections(V)
def avg_moment(settings,m):
    return sum(tr(mp(ma(*selected),m)) for selected in product(*settings))/F(4**len(settings))
def target(k,m):
    return sum(c*F(k**a,4**b) for (a,b),c in power[m].items())
assert target(4,4)==F(257,8) and target(4,6)==F(5441,32)
def pair_defect(A,B):
    return sum(tr(mp(ma(mm(mm(a,b),a),ms(F(-1,4),a)),2)) for a,b in product(A,B))/16
def triple_defect(A,B,C):
    out=F(0)
    for a,b,c in product(A,B,C):
        delta=ma(mm(mm(a,ma(mm(b,c),mm(c,b))),a),ms(F(-1,8),a))
        out+=tr(mp(delta,2))
    return out/64
matrix_tests=[]
for settings in ([P,P],[P,Q,R]):
    k=len(settings);defect=sum(pair_defect(settings[x],settings[y]) for x in range(k) for y in range(x+1,k))
    gap=avg_moment(settings,4)-target(k,4)
    assert gap==2*defect
    matrix_tests.append({'order':4,'settings':k,'gap':str(gap),'pair_defect_sum':str(defect)})
assert matrix_tests[0]['gap']!='0'
for A,B in ((P,Q),(P,R),(Q,R)):assert pair_defect(A,B)==0
delta=triple_defect(P,Q,R);gap=avg_moment([P,Q,R],6)-target(3,6)
assert delta>0 and gap==3*delta
assert triple_defect(Q,P,R)==triple_defect(R,P,Q)==delta
matrix_tests.append({'order':6,'settings':3,'gap':str(gap),'each_anchor_defect_average':str(delta),'strictly_positive':True})
def kron(A,B):return [[cm(A[i//2][j//2],B[i%2][j%2]) for j in range(4)] for i in range(4)]
X=[[Z,O],[O,Z]];Y=[[Z,cs(-1,I)],[I,Z]];PauliZ=[[O,Z],[Z,cs(-1,O)]]
XY=kron(X,Y);YZ=kron(Y,PauliZ);ZX=kron(PauliZ,X)
assert mm(XY,YZ)==ms(-1,ZX) and mm(YZ,XY)==ms(-1,ZX)
E=[ms(F(1,4),ma(Id,ms(s,XY),ms(t,YZ),ms(-s*t,ZX))) for s,t in product((-1,1),repeat=2)]
assert ma(*E)==Id
for a,b in product(range(4),repeat=2):assert mm(E[a],E[b])==(E[a] if a==b else ms(0,Id))
settings=[P,Q,R,E]
for a in range(4):
    for b in range(a+1,4):assert pair_defect(settings[a],settings[b])==0
defect=F(0)
for a in range(4):
    for b in range(a+1,4):
        for c in range(b+1,4):
            A,B,C=[settings[x] for x in (a,b,c)]
            defect+=triple_defect(A,B,C)+triple_defect(B,A,C)+triple_defect(C,A,B)
gap=avg_moment(settings,6)-target(4,6)
assert gap==defect and defect>0
matrix_tests.append({'order':6,'settings':4,'gap':str(gap),'sum_of_all_twelve_anchored_defects':str(defect),'strictly_positive':True})
report={'arithmetic':'Exact rational polynomial and Gaussian-rational matrices, standard library only',
'symbolic_parameters':'k arbitrary; h=1/n; no finite parameter interpolation',
'exhaustive_surjective_word_counts':counts,'word_tables':tables,'matrix_tests':matrix_tests,
'conclusion':'Every spectral k-UM is pairwise MUM; for k>=3 every subtriple is an algebraic 3-UM. At four outcomes/rank two, spectral and algebraic 4-UMs coincide by the common-quaternionic theorem.',
'scope':'Moment identities at arbitrary rank. Four-setting nonexistence at rank two is checked by the global and commuting-pair verifiers; higher-rank equivalence is not claimed.'}
Path(__file__).with_name('spectral_moment_rigidity_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
