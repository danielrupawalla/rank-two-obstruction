"""Independent certificate verifier: Python integers only; no NumPy.

Rebuilds the fixed triple from its defining vectors, reconstructs every
selected coefficient row by inserting each Hermitian basis entry into
left/right products, and computes an exact integer determinant.
"""
import itertools
import json
from pathlib import Path

def add(a,b): return (a[0]+b[0],a[1]+b[1])
def mul(a,b): return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def conj(a): return (a[0],-a[1])
def zero(n): return [[(0,0) for _ in range(n)] for _ in range(n)]
def identity(n):
    a=zero(n)
    for i in range(n): a[i][i]=(1,0)
    return a
def mm(a,b):
    n=len(a); out=zero(n)
    for i,j in itertools.product(range(n),repeat=2):
        for k in range(n): out[i][j]=add(out[i][j],mul(a[i][k],b[k][j]))
    return out

E=[identity(2),[[(0,0),(0,1)],[(0,1),(0,0)]],
   [[(0,0),(-1,0)],[(1,0),(0,0)]],
   [[(0,1),(0,0)],[(0,0),(0,-1)]]]
H=[[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]]
fixed=[[],[],[]]
for a in range(4):
    p=zero(8)
    p[2*a][2*a]=p[2*a+1][2*a+1]=(4,0)
    fixed[0].append(p)
for setting,units in [(1,[identity(2)]*4),(2,E)]:
    for outcome in range(4):
        p=zero(8)
        # p = v v*, where v consists of four signed 2x2 blocks.
        for a,b,i,j,k in itertools.product(range(4),range(4),range(2),range(2),range(2)):
            term=mul(units[a][i][k],conj(units[b][j][k]))
            term=mul((H[a][outcome]*H[b][outcome],0),term)
            p[2*a+i][2*b+j]=add(p[2*a+i][2*b+j],term)
        fixed[setting].append(p)

basis=[[(i,i,(1,0))] for i in range(8)]
for i,j in itertools.combinations(range(8),2):
    basis.append([(i,j,(1,0)),(j,i,(1,0))])
    basis.append([(i,j,(0,1)),(j,i,(0,-1))])

certificate=json.loads(Path(__file__).with_name('extension_linear_certificate.json').read_text())
cache={}; rebuilt=[]
for label,stored in zip(certificate['row_labels'],certificate['independent_integer_rows']):
    key=(tuple(label['outcomes']),label['anchor'])
    if key not in cache:
        selected=[fixed[x][a] for x,a in enumerate(key[0])]
        A=selected[key[1]]
        other=[selected[x] for x in range(3) if x!=key[1]]
        pairs=[]
        for order in itertools.permutations(range(3)):
            before=[]; after=[]; passed=False
            for index in order:
                if index==2: passed=True
                elif passed: after.append(other[index])
                else: before.append(other[index])
            L=A
            for factor in before: L=mm(L,factor)
            R=identity(8)
            for factor in after+[A]: R=mm(R,factor)
            pairs.append((L,R))
        # Verify independently that T=I/4 satisfies the affine equations.
        constant=zero(8)
        for L,R in pairs:
            product=mm(L,R)
            for i,j in itertools.product(range(8),repeat=2):
                constant[i][j]=add(constant[i][j],product[i][j])
        assert all(constant[i][j]==mul((24,0),A[i][j])
                   for i,j in itertools.product(range(8),repeat=2))
        cache[key]=pairs
    i,j=divmod(label['entry'],8)
    component=0 if label['component']=='real' else 1
    row=[]
    for entries in basis:
        value=(0,0)
        for L,R in cache[key]:
            for u,v,entry in entries:
                value=add(value,mul(mul(L[i][u],entry),R[v][j]))
        row.append(value[component])
    assert row==stored, label
    rebuilt.append(row)

def determinant_bareiss(matrix):
    a=[row[:] for row in matrix]; n=len(a); sign=1; previous=1
    for k in range(n-1):
        if a[k][k]==0:
            swap=next((i for i in range(k+1,n) if a[i][k]),None)
            if swap is None: return 0
            a[k],a[swap]=a[swap],a[k]; sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=a[i][j]*pivot-a[i][k]*a[k][j]
                assert numerator%previous==0
                a[i][j]=numerator//previous
            a[i][k]=0
        previous=pivot
    return sign*a[-1][-1]

det=determinant_bareiss(rebuilt)
assert len(rebuilt)==64 and det!=0
report={'verified_rows':64,'distinct_anchored_equations':len(cache),
        'arithmetic':'Python arbitrary-precision integers, independent implementation',
        'exact_determinant':str(det),'determinant_mod_1000003':det%1000003,
        'conclusion':'For the fixed published P,Q,R, the anchored algebraic 4-UM equations uniquely force T=I/4. This is not a projector, so no algebraic 4-UM extension exists.'}
Path(__file__).with_name('extension_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
