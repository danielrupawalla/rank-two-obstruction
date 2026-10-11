"""Exact Pauli-factor checks supporting the analytic fixed-pair normal form.

The exhaustive frame argument is analytic in the manuscript. This checker verifies
the two commuting spin factors, their projector expansions and cross trace.
Standard-library Gaussian integers only; no numerical ranks or optimizers.
"""
import contextlib,io,itertools,json
from pathlib import Path
with contextlib.redirect_stdout(io.StringIO()):
    import verify_extension_certificate as v
    import verify_family_reductions as f
def kron(A,B):return [[v.mul(A[i//len(B)][j//len(B)],B[i%len(B)][j%len(B)]) for j in range(len(A)*len(B))] for i in range(len(A)*len(B))]
def scale(z,A):return [[v.mul(z,a) for a in row] for row in A]
def add(A,B):return [[v.add(a,b) for a,b in zip(x,y)] for x,y in zip(A,B)]
def tr(A):return sum(A[i][i][0] for i in range(len(A))),sum(A[i][i][1] for i in range(len(A)))
p=[v.identity(2),[[(0,0),(1,0)],[(1,0),(0,0)]],
   [[(0,0),(0,-1)],[(0,1),(0,0)]],[[(1,0),(0,0)],[(0,0),(-1,0)]]]
plus=[kron(p[0],p[2]),kron(p[2],p[3]),kron(p[2],p[1])]
minus=[kron(p[3],p[2]),kron(p[2],p[0]),scale((-1,0),kron(p[1],p[2]))]
I=v.identity(4)
for L in [plus,minus]:
    for a in range(3):
        assert v.mm(L[a],L[a])==I
        assert tr(L[a])==(0,0)
        assert all(z[0]==0 for row in L[a] for z in row)
        assert [[z[1] for z in row] for row in L[a]]==[[-L[a][j][i][1] for j in range(4)] for i in range(4)]
        assert v.mm(L[a],L[(a+1)%3])==scale((0,1),L[(a+2)%3])
for A,B in itertools.product(plus,minus):
    assert v.mm(A,B)==v.mm(B,A)
    assert tr(v.mm(A,B))==(0,0)
for L,R0 in [(plus,v.fixed[2][0]),(minus,f.minus_measurement()[0])]:
    rebuilt=v.identity(8)
    for a,sign in enumerate([1,-1,1]):rebuilt=add(rebuilt,scale((sign,0),kron(L[a],p[a+1])))
    assert rebuilt==R0
for L,S in [(plus,v.fixed[2]),(minus,f.minus_measurement())]:
    observed=set()
    for A in S:
        signs=[];rebuilt=v.identity(8)
        for a in range(3):
            r,i=tr(v.mm(kron(L[a],p[a+1]),A));assert i==0 and r in [-8,8]
            signs.append(r//8);rebuilt=add(rebuilt,scale((r//8,0),kron(L[a],p[a+1])))
        assert rebuilt==A and signs[0]*signs[1]*signs[2]==-1
        observed.add(tuple(signs))
    assert len(observed)==4
# Numerators have denominator four, so 8/16=1/2.
assert tr(v.mm(v.fixed[2][0],f.minus_measurement()[0]))==(8,0)
# All mixed spin-component traces vanish, proving the cross-type overlap
# stays 1/2 under arbitrary rotations of the second-factor Pauli axes.
for A,B in itertools.product(plus,minus):
    for a,b in itertools.product(range(1,4),repeat=2):
        assert tr(v.mm(kron(A,p[a]),kron(B,p[b])))==(0,0)
report={'arithmetic':'Python arbitrary-precision Gaussian integers',
 'spin_factor_labels':[['IY','YZ','YX'],['ZY','YI','-XY']],
 'both_spin_factors_obey_Pauli_relations':True,
 'spin_generators_are_pure_imaginary_skew':True,
 'two_factors_commute_and_are_trace_orthogonal':True,
 'both_base_projectors_have_Bell_expansion':True,
 'both_four_outcome_settings_have_all_four_Bell_sign_patterns':True,
 'cross_type_projector_overlap':'1/2',
 'all_rotated_cross_component_traces_checked':81,
 'scope':'Finite structural identities only. Exhaustiveness and maximally-entangled-basis reduction are proved analytically in the manuscript.'}
Path(__file__).with_name('fixedpair_structure_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
