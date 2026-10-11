"""Exact complete strong triple with one noncommuting and two commuting pairs.

Checks the full quaternionic operator equations, hence their C^8 realization.
Illustrates the global commuting-pair theorem: the triple contains commuting canonical pairs although
its first pair does not commute. Standard library only.
"""
from fractions import Fraction as F
from pathlib import Path
import itertools,json
from helpers import *
P=[[[O if i==j==a else Z for j in range(4)] for i in range(4)] for a in range(4)];Q=projections(U);R=projections(W);settings=[P,Q,R]
for setting in settings:
 assert __import__('functools').reduce(ma,setting)==identity
 for a,b in itertools.product(range(4),repeat=2):assert mm(setting[a],setting[b])==setting[a] if a==b else mm(setting[a],setting[b])==ms(0,identity)
for x,y in itertools.permutations(range(3),2):
 for a,b in itertools.product(range(4),repeat=2):assert mm(mm(settings[x][a],settings[y][b]),settings[x][a])==ms(F(1,4),settings[x][a])
for a,b,c in itertools.product(range(4),repeat=3):
 mats=[P[a],Q[b],R[c]]
 for anchor in range(3):
  A=mats[anchor];B,C=[mats[t] for t in range(3) if t!=anchor]
  assert mm(mm(A,ma(mm(B,C),mm(C,B))),A)==ms(F(1,8),A)
S=[[qs(F(1,2),qa(*(qm(qc(U[a][b]),W[a][c]) for a in range(4)))) for c in range(4)] for b in range(4)]
def canonical(M):return [[qm(qm(qm(qc(M[a][0]),M[a][b]),qc(M[0][b])),M[0][0]) for b in range(4)] for a in range(4)]
commuting=[];nonreal=[]
for M in [U,W,S]:
 C=canonical(M);flat=sum(C,[]);commuting.append(all(qm(a,b)==qm(b,a) for a,b in itertools.combinations(flat,2)));witness=next(((a,b,C[a][b]) for a,b in itertools.product(range(4),repeat=2) if any(C[a][b][1:])),None);assert witness;nonreal.append({'row':witness[0],'column':witness[1],'canonical_quaternion':list(map(str,witness[2]))})
assert commuting==[False,True,True]
invariants=[]
for c in range(4):
 u=[W[a][c] for a in range(4)];s=[S[b][c] for b in range(4)];fp=sum(q[0]**2 for q in u[1:]);fq=sum(dot(q,s[0])**2 for q in s[1:]);assert fp==fq==F(1,2)
 pp=[sum(qm(qc(U[a][b]),u[a])[0]**2 for a in range(1,4)) for b in range(4)]
 qq=[sum(qm(qm(U[a][b],s[b]),qc(s[0]))[0]**2 for b in range(1,4)) for a in range(4)]
 assert all(z==F(1,2) for z in pp+qq)
 fullP=[];fullQ=[]
 for a,b in itertools.product(range(4),repeat=2):
  v=[qm(qc(U[i][b]),u[i]) for i in range(4)];tv=[qm(U[a][j],s[j]) for j in range(4)]
  fullP.append(sum(dot(v[i],v[a])**2 for i in range(4) if i!=a));fullQ.append(sum(dot(tv[j],tv[b])**2 for j in range(4) if j!=b))
 assert all(z==F(1,2) for z in fullP+fullQ)
 invariants.append({'outcome':c,'fP':str(fp),'fQ':str(fq),'all_eight_pivot_invariants':list(map(str,pp+qq)),'all_sixteen_P_and_Q_base_pairs_checked':True})
report={'arithmetic':'Exact rational quaternion matrices, standard library only','pvm_settings_checked':3,'ordered_pair_identities_checked':96,'strong_anchor_identities_checked':192,'canonical_pair_commutativity_PQ_PR_QR':commuting,'nonreal_canonical_witnesses':nonreal,'split_invariants':invariants,'conclusion':'A complete strong triple whose first pair is noncommuting. Its other two pairs commute canonically, as the global commuting-pair theorem requires, and the uniform nonextension theorem applies.','scope':'One exact example; the general statements are proved in the manuscript.'}
Path(__file__).with_name('mixed_strong_triple_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
