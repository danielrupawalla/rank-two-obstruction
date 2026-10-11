"""Exact checks supporting the four-outcome Hadamard--Clifford extension.

Standard library only. Verifies the original injectivity certificate,
the simultaneous unitary equivalence of the two irreducible triples,
and the signed-column normal form of every real order-four Hadamard.
The arbitrary-rank representation and projection arguments are in the proof.
"""
import contextlib
import io
import itertools
import json
from pathlib import Path

with contextlib.redirect_stdout(io.StringIO()):
    import verify_extension_certificate as v

def adj(a): return [[v.conj(a[j][i]) for j in range(len(a))] for i in range(len(a))]
def scaled(a,k): return [[v.mul((k,0),x) for x in row] for row in a]
def minus_measurement():
    units=v.E[:]
    units[3]=scaled(units[3],-1)
    out=[]
    for outcome in range(4):
        p=v.zero(8)
        for a,b,i,j,k in itertools.product(range(4),range(4),range(2),range(2),range(2)):
            term=v.mul(units[a][i][k],v.conj(units[b][j][k]))
            term=v.mul((v.H[a][outcome]*v.H[b][outcome],0),term)
            p[2*a+i][2*b+j]=v.add(p[2*a+i][2*b+j],term)
        out.append(p)
    return out

perm=[0,2,1,3]
# Numerator of U=(sigma_x-sigma_y)/sqrt(2).
U0=[[(0,0),(1,1)],[(1,-1),(0,0)]]
W0=v.zero(8)
for a,i,j in itertools.product(range(4),range(2),range(2)):
    W0[2*perm[a]+i][2*a+j]=U0[i][j]
assert v.mm(W0,adj(W0))==scaled(v.identity(8),2)
targets=[v.fixed[0],v.fixed[1],minus_measurement()]
for x,a in itertools.product(range(3),range(4)):
    assert v.mm(v.mm(W0,v.fixed[x][a]),adj(W0))==scaled(targets[x][perm[a]],2)

sign_vectors=list(itertools.product([-1,1],repeat=4))
def orthogonal(a,b): return sum(x*y for x,y in zip(a,b))==0
all_hadamards=set()
for rows in itertools.product(sign_vectors,repeat=4):
    if all(orthogonal(rows[i],rows[j]) for i,j in itertools.combinations(range(4),2)):
        all_hadamards.add(rows)
normal_forms=set()
for row_signs,columns,col_signs in itertools.product(sign_vectors,itertools.permutations(range(4)),sign_vectors):
    normal_forms.add(tuple(tuple(row_signs[a]*v.H[a][columns[b]]*col_signs[b]
                                for b in range(4)) for a in range(4)))
assert all_hadamards==normal_forms and len(all_hadamards)==768
report={'arithmetic':'Python arbitrary-precision Gaussian integers',
        'original_injectivity_certificate_verified':True,
        'opposite_chirality_projector_equivalences_verified':12,
        'all_real_order_four_hadamards_checked':768,
        'signed_column_normal_form_verified':True,
        'scope':'Finite reductions supporting the family theorem; representation decomposition and block-projection argument are analytic.'}
Path(__file__).with_name('family_reduction_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
