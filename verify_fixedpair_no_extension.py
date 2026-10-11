"""Exact symbolic bilinear restrictions; all parameters, no sampling."""
import itertools,json
from fractions import Fraction
from pathlib import Path
single={('X','Y'):('Z',1),('Y','Z'):('X',1),('Z','X'):('Y',1),('Y','X'):('Z',3),('Z','Y'):('X',3),('X','Z'):('Y',3)}
labels=[''.join(x) for x in itertools.product('IXYZ',repeat=3)]
table={}
for A,B in itertools.product(labels,repeat=2):
    out='';phase=0
    for a,b in zip(A,B):
        if a=='I':out+=b
        elif b=='I':out+=a
        elif a==b:out+='I'
        else:z,p=single[a,b];out+=z;phase+=p
    table[A,B]=out,[(1,0),(0,1),(-1,0),(0,-1)][phase%4]
def cm(x,y):return x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def add(A,B):
    C=A.copy()
    for k,(r,i) in B.items():a,b=C.get(k,(0,0));C[k]=(a+r,b+i)
    return {k:v for k,v in C.items() if v!=(0,0)}
def mul(A,B):
    C={}
    for a,x in A.items():
        for b,y in B.items():
            z,p=table[a,b];r,i=cm(cm(x,y),p);u,v=C.get(z,(0,0));C[z]=u+r,v+i
    return {k:v for k,v in C.items() if v!=(0,0)}
def chain(*args):
    A={'III':(1,0)}
    for B in args:A=mul(A,B)
    return A
def outcome(u,v,x,y):return {'III':(1,0),x:(u,0),y:(v,0),table[x,y][0]:(u*v,0)}
P=[outcome(u,v,'ZII','IZI') for u,v in itertools.product([1,-1],repeat=2)]
Q=[outcome(u,v,'XII','IXI') for u,v in [(1,1),(1,-1),(-1,1),(-1,-1)]]
def rref(A):
    A=[[Fraction(x) for x in row] for row in A];n=0
    for j in range(len(A[0])):
        k=next((k for k in range(n,len(A)) if A[k][j]),None)
        if k is None:continue
        A[n],A[k]=A[k],A[n];q=A[n][j];A[n]=[x/q for x in A[n]]
        for k in range(len(A)):
            if k!=n and A[k][j]:q=A[k][j];A[k]=[x-q*y for x,y in zip(A[k],A[n])]
        n+=1
    return A[:n]
report={}
for kind,L in [('same',['IY','YZ','YX']),('opposite',['ZY','YI','XY'])]:
    ranks=[]
    for j,sigma in enumerate('XYZ'):
        rows=[]
        for a,b in itertools.product(range(4),repeat=2):
            A=P[a];B=Q[b];columns=[]
            for i,k,l in itertools.product(range(3),repeat=3):
                C={['IY','YZ','YX'][i]+sigma:(1,0)};D={L[k]+'XYZ'[l]:(1,0)};total={}
                for order in itertools.permutations([B,C,D]):total=add(total,chain(A,*order,A))
                assert all(v[1]==0 for v in total.values())
                columns.append(total)
            for out in labels:rows.append([C.get(out,(0,0))[0] for C in columns])
        reduced=rref(rows)
        if kind=='same':
            expected=[]
            row=[0]*27
            for i in range(3):row[i*9+i*3+j]=1
            expected.append(row)
            for i,k in itertools.combinations(range(3),2):
                row=[0]*27;row[i*9+k*3+j]=row[k*9+i*3+j]=1;expected.append(row)
        else:
            expected=[]
            for i,k in itertools.product(range(3),repeat=2):
                row=[0]*27;row[i*9+k*3+j]=1;expected.append(row)
        assert reduced==rref(expected)
        ranks.append(len(reduced))
    report[kind]={'ranks_by_column':ranks,'restrictions':'m.x=0 and m_i x_k+m_k x_i=0' if kind=='same' else 'm_i x_k=0 for all i,k'}
# A row swap exchanges chiralities and permutes the fixed Hadamard outcomes.
H=[[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]]
swap=[0,2,1,3]
assert [[H[swap[a]][b] for a in range(4)] for b in range(4)]==[[H[a][b] for a in range(4)] for b in [0,2,1,3]]
report['conclusion']='All fixed-pair triples are nonextendible: the opposite sector vanishes; a same-sector column can survive only when m has exactly two equal-magnitude nonzero coordinates. An orthogonal 3x3 M cannot have three such mutually orthogonal columns, forcing a zero column in the candidate orthogonal correlation matrix X.'
report['arithmetic']='Exact Gaussian integer Pauli products and rational row reduction'
Path(__file__).with_name('fixedpair_no_extension_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
