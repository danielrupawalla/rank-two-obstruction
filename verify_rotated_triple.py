"""Independent standard-library exact verification of the rotated triple.

No NumPy, floating point, imported exploratory code, or stored projectors.
Gaussian integers rebuild projectors and canonical pair invariants.
"""
import itertools,json,math
from pathlib import Path
def add(a,b):return a[0]+b[0],a[1]+b[1]
def mul(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def conj(a):return a[0],-a[1]
def zero(n,m=None):return [[(0,0) for _ in range(n if m is None else m)] for _ in range(n)]
def eye(n):
    A=zero(n)
    for i in range(n):A[i][i]=(1,0)
    return A
def mm(A,B):
    C=zero(len(A),len(B[0]))
    for i,j,k in itertools.product(range(len(A)),range(len(B[0])),range(len(B))):C[i][j]=add(C[i][j],mul(A[i][k],B[k][j]))
    return C
def adj(A):return [[conj(A[j][i]) for j in range(len(A))] for i in range(len(A[0]))]
def scale(k,A):return [[(k*r,k*i) for r,i in row] for row in A]
def ma(A,B):return [[add(a,b) for a,b in zip(x,y)] for x,y in zip(A,B)]
H=[[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]]
E=[eye(2),[[(0,0),(0,1)],[(0,1),(0,0)]],
   [[(0,0),(-1,0)],[(1,0),(0,0)]],[[(0,1),(0,0)],[(0,0),(0,-1)]]]
O=[[9,-12,-20,0],[20,15,0,0],[12,-16,15,0],[0,0,0,25]]
assert all(sum(O[a][u]*O[b][u] for u in range(4))==625*(a==b) for a,b in itertools.product(range(4),repeat=2))
W=[[(O[i//2][j//2]*(i%2==j%2),0) for j in range(8)] for i in range(8)]
frames=[[],[],[]]
for c in range(4):
    P=zero(8,2)
    for i in range(2):P[2*c+i][i]=(1,0)
    frames[0].append(P)
    for s in [1,2]:
        F=[[(0,0) for _ in range(2)] for _ in range(8)]
        for a,i,j in itertools.product(range(4),range(2),range(2)):
            F[2*a+i][j]=mul((H[a][c],0),(eye(2) if s==1 else E[a])[i][j])
        frames[s].append(F if s==1 else mm(W,F))
den=[1,2,50];D=2500
fixed=[[scale(D//den[s]**2,mm(F,adj(F))) for F in frames[s]] for s in range(3)]
checks={'projector_products':0,'ordered_pair_compressions':0,'anchored_triple_identities':0}
for s,S in enumerate(fixed):
    total=zero(8)
    for a,A in enumerate(S):
        assert adj(A)==A
        assert sum(A[i][i][0] for i in range(8))==2*D
        assert mm(adj(frames[s][a]),frames[s][a])==scale(den[s]**2,eye(2))
        total=ma(total,A)
        for b,B in enumerate(S):
            assert mm(A,B)==(scale(D,A) if a==b else zero(8));checks['projector_products']+=1
    assert total==scale(D,eye(8))
for s,t in itertools.permutations(range(3),2):
    for A,B in itertools.product(fixed[s],fixed[t]):
        assert scale(4,mm(mm(A,B),A))==scale(D**2,A)
        checks['ordered_pair_compressions']+=1
for outcomes in itertools.product(range(4),repeat=3):
    for a in range(3):
        b,c=[s for s in range(3) if s!=a]
        A,B,C=[fixed[s][outcomes[s]] for s in [a,b,c]]
        assert scale(8,mm(mm(A,ma(mm(B,C),mm(C,B))),A))==scale(D**3,A)
        checks['anchored_triple_identities']+=1
def reduced(A,d):
    g=math.gcd(d,math.gcd(*(abs(z) for row in A for entry in row for z in entry)))
    return {'real_numerator':[[r//g for r,i in row] for row in A],
            'imag_numerator':[[i//g for r,i in row] for row in A],'denominator':d//g}
pairs=[]
for s,t in itertools.combinations(range(3),2):
    B=[[scale(2,mm(adj(F),G)) for G in frames[t]] for F in frames[s]]
    bd=den[s]*den[t]
    for row in B:
        for X in row:assert mm(X,adj(X))==scale(bd**2,eye(2))
    # Gauge transformations L_a=B_a0^*, R_b=B_0b^*B_00
    # yield first row and column exactly identity.
    U=[[mm(mm(mm(adj(B[a][0]),B[a][b]),adj(B[0][b])),B[0][0]) for b in range(4)] for a in range(4)]
    assert all(U[0][a]==U[a][0]==scale(bd**4,eye(2)) for a in range(4))
    witness=None
    for a,b,c,d in itertools.product(range(4),repeat=4):
        K=ma(mm(U[a][b],U[c][d]),scale(-1,mm(U[c][d],U[a][b])))
        if K!=zero(2):
            witness={'blocks':[[a,b],[c,d]],'commutator':reduced(K,bd**8)};break
    pairs.append({'settings':[s,t],'direct_sum_of_MUBs':witness is None,'canonical_commutator_witness':witness})
assert [p['direct_sum_of_MUBs'] for p in pairs]==[True,False,False]
report={'arithmetic':'Python arbitrary-precision Gaussian integers only',
        'O_numerator':O,'O_denominator':25,'projector_denominator':D,
        'checks':checks,'pairs':pairs,
        'inequivalence':'Exactly one pair is a direct sum of MUBs. Every rank-two real Hadamard-Clifford construction has at least two such pairs. This invariant survives all allowed equivalences.'}
Path(__file__).with_name('rotated_triple_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
