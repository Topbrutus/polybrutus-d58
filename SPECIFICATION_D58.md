# Polybrutus D58 — 63-point exact Euclidean squared-distance metric

**Sole direct original scientific parent:** Topbrutus/polybrutus-d57@v0.1.0, peeled scientific Git tag \`58bdbb2bd6a045c3db8ec42c0721a5833cf2810e\`, DOI 10.5281/zenodo.24431122, original complete source archive SHA256 \`fcdb7c69e872c05c3e54021db500c24836dbb3aed2c3f36160db584398de703e\`, byte-identical original certificate SHA256 \`c4dd7d4e1b00a46c388b0648d8308f6578b9705e025967b8a92d5e367d9e8da4\`.

## Conservative perpendicular unit-vector extension

For the 62 original norm-one predecessor points p_i, define q_i=(p_i,0), then append perpendicular norm-one q_63. All C(62,2)=1891 predecessor labelled squared distances remain unchanged. All 62 new squared distances equal 2R²; total C(63,2)=1953 unordered labelled pairs.

Construct exact rational Gram G as orthogonal block sum of seven canonical norm-one points rank 3, a five-point centered regular unit simplex with pairwise dot products -1/4 rank 4, and **51** independent mutually orthogonal norm-one axes. Hence G=XXᵀ is PSD, has unit diagonal and rank **3+4+51=58**. Define D_ii=0 and D_ij=2-2G_ij. Let H=I-11ᵀ/63. The Euclidean centered Gram B=-(HDH)/2=HGH is PSD exact rank **58**. An inherited canonical antipodal pair implies the origin is in the affine hull, establishing minimum affine Euclidean realization dimension **58** for this specific metric, using 59 coordinates for storage.

At unit radius the exhaustive histogram k=4d²/R² is 4:6, 7:2, 8:1922, 9:3, 10:10, 12:6, 13:2, 16:2, totaling 1953 pairs. Exact sum of all squared pair distances R=1: **15671/4**.

## Independently chosen ideal weighted complete graph

Choose unit scalar masses and symmetric positive complete graph edges w_ij=1/d_ij². New vertex joins all 62 old vertices with weight 1/2. The dimensionless ideal Laplacian L58=[[L57+I62/2,-1_62/2],[-1_62ᵀ/2,31]].

Characteristic recurrence chi63(z)=z(z-63/2)chi62(z-1/2)/(z-1/2). With the original seven-point f6(z)=chiP7(z)/z, closed form chi63(z)=z(z-63/2)^52(z-31)^4*f6(z-28). Rational Laplacian rank-nullity establishes eigenvalue **63/2 multiplicity 52** and **31 multiplicity 4**. These are model eigenvalues without dimension, not measured physical hertz.

**References:** direct D57 DOI 10.5281/zenodo.24431122, historical D7 DOI 10.5281/zenodo.24407952.
