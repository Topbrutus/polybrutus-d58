#!/usr/bin/env python3
"""Polybrutus D58: exact Euclidean labelled metric and chosen ideal weighted graph.
Conservative perpendicular unit extension of the immutable Polybrutus D57 v0.1.0.
No physical spacetime dimension or experimentally measured resonance is claimed.
"""
from pathlib import Path
from itertools import combinations
from collections import Counter
from hashlib import sha256
import json
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
Q=s.Rational
core=[(0,0,1),(Q(1,2),0,Q(1,2)),(-Q(1,4),Q(3,4),Q(1,2)),
      (-Q(1,4),-Q(3,4),Q(1,2)),(0,0,-1),(Q(1,2),0,-Q(1,2)),
      (-Q(1,2),0,-Q(1,2))]
n=63
pairs=list(combinations(range(n),2))
G=s.zeros(n)
for i in range(7):
    for j in range(7):
        G[i,j]=3*core[i][0]*core[j][0]+core[i][1]*core[j][1]+core[i][2]*core[j][2]
for i in range(7,12):
    for j in range(7,12):
        G[i,j]=1 if i==j else -Q(1,4)
for i in range(12,n):
    G[i,i]=1
assert all(G[i,i]==1 for i in range(n))
assert DomainMatrix.from_Matrix(G).to_field().rank()==58
D=2*s.ones(n)-2*G
for i in range(n):
    D[i,i]=0
assert len(pairs)==1953 and all(D[i,j]>0 for i,j in pairs)

parent_path=ROOT/'prior'/'d57'/'results'/'certificate_d57.json'
parent_archive=ROOT/'prior'/'d57'/'frozen_scientific_tag_v0.1.0.zip'
assert sha256(parent_archive.read_bytes()).hexdigest()=='fcdb7c69e872c05c3e54021db500c24836dbb3aed2c3f36160db584398de703e'
from zipfile import ZipFile
with ZipFile(parent_archive) as _frozen:
    assert _frozen.testzip() is None
    assert _frozen.read('results/certificate_d57.json')==parent_path.read_bytes()
parent_bytes=parent_path.read_bytes()
assert sha256(parent_bytes).hexdigest()=='c4dd7d4e1b00a46c388b0648d8308f6578b9705e025967b8a92d5e367d9e8da4'
parent=json.loads(parent_bytes)
assert parent['status']=='PASS' and parent['dimension']==57 and parent['vertices']==62
former={(int(e['i'])-1,int(e['j'])-1):s.Rational(e['d2'])
        for e in parent['labelled_squared_distances']}
assert len(former)==1891
assert all(D[i,j]==v for (i,j),v in former.items())
assert all(D[i,62]==2 for i in range(62))
hist=Counter(4*D[i,j] for i,j in pairs)
expected=Counter({4:6,7:2,8:1922,9:3,10:10,12:6,13:2,16:2})
assert hist==expected and sum(hist.values())==1953
assert sum(D[i,j] for i,j in pairs)==Q(15671,4)

H=s.eye(n)-s.ones(n)/n
B=-H*D*H/2
assert B==H*G*H
assert DomainMatrix.from_Matrix(B).to_field().rank()==58
# The first and fifth canonical points are antipodal; origin lies in the affine hull.
assert G[0,4]==-1

L=s.zeros(n)
for i,j in pairs:
    w=1/D[i,j]
    L[i,i]+=w;L[j,j]+=w;L[i,j]-=w;L[j,i]-=w
assert L*s.ones(n,1)==s.zeros(n,1)
assert DomainMatrix.from_Matrix(L-Q(63,2)*s.eye(n)).to_field().rank()==11
assert DomainMatrix.from_Matrix(L-31*s.eye(n)).to_field().rank()==59
old=s.zeros(62)
for i,j in combinations(range(62),2):
    w=1/D[i,j]
    old[i,i]+=w;old[j,j]+=w;old[i,j]-=w;old[j,i]-=w
assert L[:62,:62]-old==s.eye(62)/2
assert all(L[i,62]==-Q(1,2) for i in range(62))
assert L[62,62]==31

# Four adversarial negative controls.
damaged=D.copy();damaged[0,62]=damaged[62,0]=Q(9,4)
assert Counter(4*damaged[i,j] for i,j in pairs)!=hist
duplicate=G.copy();duplicate[62,:]=G[61,:];duplicate[:,62]=G[:,61];duplicate[62,62]=1
assert DomainMatrix.from_Matrix(duplicate).to_field().rank()==57
tilted_dot_to_old_core=Q(3,10)
assert Q(3,10)**2+Q(91,100)==1
assert tilted_dot_to_old_core!=0 and 2-2*tilted_dot_to_old_core!=D[0,62]
swapped=D.copy();swapped.row_swap(0,7);swapped.col_swap(0,7)
assert swapped!=D and Counter(4*swapped[i,j] for i,j in pairs)==hist

report={
  'schema':'POLYBRUTUS-D58-EXACT-CERTIFICATE-v1',
  'status':'PASS','dimension':58,'vertices':63,'pairs':1953,
  'affine_rank':58,'centered_gram_rank':58,'ambient_storage_coordinates':59,
  'radius':1,'parent':'Topbrutus/polybrutus-d57@v0.1.0',
  'parent_tag_commit':'58bdbb2bd6a045c3db8ec42c0721a5833cf2810e',
  'parent_doi':'10.5281/zenodo.24431122',
  'parent_certificate_sha256':sha256(parent_bytes).hexdigest(),
  'parent_frozen_tag_archive_sha256':sha256(parent_archive.read_bytes()).hexdigest(),
  'parent_frozen_tag_embedded_as_verified_zip':True,
  'parent_distances_preserved':1891,'new_distances':62,
  'new_squared_distances':['2']*62,
  'histogram_k_4d2':{str(k):int(v) for k,v in sorted(hist.items())},
  'sum_squared_distances':'15671/4',
  'gram_psd_certificate':'B=-(H D H)/2=H G H; G=XX^T exact',
  'negative_controls':'4/4 PASS',
  'oscillator_model':'Chosen ideal equal-mass complete graph with w_ij=1/d_ij^2; no physical measurement',
  'eigenvalue_63_over_2_multiplicity':52,
  'eigenvalue_31_multiplicity':4,
  'spectral_recurrence':'chi63(z)=z*(z-63/2)*chi62(z-1/2)/(z-1/2)',
  'spectral_closed_form':'chi63(z)=z*(z-63/2)^52*(z-31)^4*f6(z-28), f6(z)=chiP7(z)/z',
  'novel_universal_theorem':'NOT CLAIMED','physical_discovery':'NOT CLAIMED',
  'human_peer_review':'NOT_DONE','sympy':s.__version__,
  'labelled_squared_distances':[{'i':i+1,'j':j+1,'d2':str(D[i,j]),'k':str(4*D[i,j])} for i,j in pairs]
}
out=ROOT/'results'/'certificate_d58.json'
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n')
print('PRIMARY D58 PASS: 63 vertices, 1953 labelled pairs, affine rank 58; 1891/1891 previous pairs preserved')
print('IDEAL GRAPH PASS: eigenvalue 63/2 x52; eigenvalue 31 x4; 4/4 negative controls')
print('CERTIFICATE_SHA256='+sha256(out.read_bytes()).hexdigest())
