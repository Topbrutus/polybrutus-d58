#!/usr/bin/env python3
from pathlib import Path
import json
R=Path(__file__).resolve().parents[1]
c=json.loads((R/'results/certificate_d58.json').read_text('utf8'))
m=json.loads((R/'publication/zenodo_metadata_blueprint.json').read_text('utf8'))
assert c['status']=='PASS' and c['dimension']==m['scientific_dimension']==58
assert c['vertices']==m['norm_one_vertices']==63
assert c['pairs']==m['exact_labelled_squared_pairs']==1953
assert c['affine_rank']==c['centered_gram_rank']==58
assert c['eigenvalue_63_over_2_multiplicity']==m['model_eigenvalue_63_over_2_multiplicity']==52
assert c['eigenvalue_31_multiplicity']==m['model_eigenvalue_31_multiplicity']==4
assert c['parent_tag_commit']=='58bdbb2bd6a045c3db8ec42c0721a5833cf2810e'
assert c['parent_doi']==m['scientific_parent_doi']=='10.5281/zenodo.24431122'
assert m['related_identifiers']==[{'identifier':'10.5281/zenodo.24431122','relation':'isDerivedFrom','resource_type':'publication-preprint','scheme':'doi'}]
assert len(m['keywords'])==9
assert 'affine Euclidean dimension fifty-eight' in m['keywords']
assert 'fifty-twofold ideal eigenvalue' in m['keywords']
assert len(m['references'])==2
assert '10.5281/zenodo.24431122' in m['references'][0]
assert '10.5281/zenodo.24407952' in m['references'][1]
assert not any(x in str(m['keywords']) for x in ('fifty-seven','fifty-onefold','fifty-nine','fifty-threefold'))
note=(R/'publication/POLYBRUTUS_58D_NOTE.md').read_text('utf8')
assert all(x in note for x in ('## References','10.5281/zenodo.24431122','10.5281/zenodo.24407952'))
print('D58_CURRENT_9_KEYWORDS_STRUCTURED_DIRECT_D57_PARENT_AND_SEPARATE_D57_D7_REFERENCES_PASS')
