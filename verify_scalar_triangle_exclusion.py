"""Exact finite-field check of the determinant-phase obstruction.

The analytic step deriving cube-root phases and balanced line sums from
unitary frames is proved in the manuscript, not assumed from sampling.
"""
import itertools,json
from pathlib import Path
def rank(A,p=3):
    A=[[x%p for x in row] for row in A];k=0
    for col in range(len(A[0])):
        row=next((i for i in range(k,len(A)) if A[i][col]),None)
        if row is None:continue
        A[k],A[row]=A[row],A[k];q=pow(A[k][col],-1,p);A[k]=[x*q%p for x in A[k]]
        for i in range(len(A)):
            if i!=k and A[i][col]:q=A[i][col];A[i]=[(x-q*y)%p for x,y in zip(A[i],A[k])]
        k+=1
    return k
report=[]
for n in [2,3,4,5]:
    indices=list(itertools.product(range(n),repeat=3));G=[]
    for a,b,c in indices:
        row=[0]*(3*n*n);row[a*n+b]=1;row[n*n+b*n+c]=1;row[2*n*n+a*n+c]=-1;G.append(row)
    sums=[]
    for axis in range(3):
        for x,y in itertools.product(range(n),repeat=2):
            chosen=[row for abc,row in zip(indices,G) if tuple(abc[k] for k in range(3) if k!=axis)==(x,y)]
            sums.append([sum(row[j] for row in chosen) for j in range(3*n*n)])
    rg=rank(G);rs=rank(sums)
    if n%3:assert rg==rs==3*n*n-3*n+1
    else:assert rs<rg
    report.append({'n':n,'phase_array_rank':rg,'line_balance_rank':rs,'line_balance_forces_zero_phase':rg==rs})
assert next(x for x in report if x['n']==4)['phase_array_rank']==37
out={'arithmetic':'Exact Gaussian elimination over F_3','checks':report,
     'four_outcome_conclusion':'Balanced line sums force any phase array z_abc=u_ab+s_bc-w_ac to vanish. Scalar triangles would have phase +/-1, so none can occur.',
     'analytic_requirements':'Unitary rank-two triangles have trace 1 or 1+/-i sqrt(3); completeness gives all line trace sums 4; dephasing gives cube-root block determinants. These steps are proved in the manuscript.'}
Path(__file__).with_name('scalar_triangle_exclusion_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
