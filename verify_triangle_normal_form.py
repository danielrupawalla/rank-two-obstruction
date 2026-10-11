"""Exact unitary-block triangle checks on the published and rotated triples.

The general reduction is analytic; this verifies its constants and the
nonscalar triangle property on both exact constructions.
"""
import contextlib,io,itertools,json
from pathlib import Path
with contextlib.redirect_stdout(io.StringIO()):import verify_rotated_triple as r
def block(F,a):return F[2*a:2*a+2]
def det(A):return r.add(r.mul(A[0][0],A[1][1]),(-r.mul(A[0][1],A[1][0])[0],-r.mul(A[0][1],A[1][0])[1]))
reports=[]
for name,O in [('published',None),('rotated',r.W)]:
    U=[[r.scale(r.H[a][b],r.eye(2)) for b in range(4)] for a in range(4)]
    Wraw=[]
    for a in range(4):
        row=[]
        for c in range(4):
            if O is None:row.append(r.scale(r.H[a][c],r.E[a]))
            else:row.append(block(r.frames[2][c],a))
        Wraw.append(row)
    wd=1 if O is None else 25
    W=[[r.mm(Wraw[a][c],r.adj(Wraw[0][c])) for c in range(4)] for a in range(4)]
    dW=wd**2
    assert all(W[0][c]==r.scale(dW,r.eye(2)) for c in range(4))
    S=[]
    for b in range(4):
        row=[]
        for c in range(4):
            A=r.zero(2)
            for a in range(4):A=r.ma(A,r.mm(r.adj(U[a][b]),W[a][c]))
            row.append(A)
        S.append(row)
    dS=2*dW
    assert all(det(X)==(dW*dW,0) for row in W for X in row)
    assert all(det(X)==(dS*dS,0) for row in S for X in row)
    noncentral=0
    for a,b,c in itertools.product(range(4),repeat=3):
        Z=r.mm(r.mm(U[a][b],S[b][c]),r.adj(W[a][c]));dZ=dS*dW
        assert r.ma(Z,r.adj(Z))==r.scale(dZ,r.eye(2))
        assert r.mm(Z,r.adj(Z))==r.scale(dZ*dZ,r.eye(2))
        assert det(Z)==(dZ*dZ,0)
        assert not (Z[0][1]==Z[1][0]==(0,0) and Z[0][0]==Z[1][1])
        noncentral+=1
    reports.append({'triple':name,'noncentral_triangles':noncentral,'unitarity_hermitian_part_and_determinants_verified':True})
report={'arithmetic':'Python arbitrary-precision Gaussian integers','constructions':reports,
 'scope':'Checks constants and hypotheses at two constructions. The general normal form is analytic; scalar triangles are excluded separately by verify_scalar_triangle_exclusion.py.'}
Path(__file__).with_name('triangle_normal_form_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
