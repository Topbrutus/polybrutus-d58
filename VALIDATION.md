# Polybrutus D58 exact verification and publication protocol

1. Authenticate direct frozen scientific parent D57@v0.1.0 peeled commit 58bdbb2bd6a045c3db8ec42c0721a5833cf2810e, DOI 10.5281/zenodo.24431122, complete original source ZIP SHA256 fcdb7c69e872c05c3e54021db500c24836dbb3aed2c3f36160db584398de703e, original parent certificate SHA256 c4dd7d4e1b00a46c388b0648d8308f6578b9705e025967b8a92d5e367d9e8da4; CRC, original embedded certificate byte equality and all labels.
2. Independently verify 63 norm-one points, all 1953 unordered labelled squared distances, all 1891 old D57 label values, 62 new 2R² squared distances, Gram PSD and centered Euclidean affine rank 58.
3. Two methods (SymPy Gram vs Python fractions.Fraction independent distance rules), four negative controls, same deterministic JSON certificate. Separately stipulated ideal graph eigenvalues 63/2 x52, 31 x4.
4. GitHub three CI success required (main initial, immutable scientific tag v0.1.0, post-Zenodo documentary main), checking two verifiers, deterministic certificate and metadata QA.
5. Zenodo three original scientific tag ZIP, two-page PDF, exact JSON, public MD5/names/sizes against SHA256-verified originals, nine D58 keywords, structured isDerivedFrom immediate D57 DOI, references separately D57 and historical D7 DOI 10.5281/zenodo.24407952.
6. Preserve unmodified original ancestor ZIPs and disclose D47 historical D46 wrong grandparent textual SHA; do not infer physical dimensions, hertz or new universal theorem.

    python -m pip install sympy==1.14.0
    python src/verify_d58.py
    python tests/independent_d58.py
    python src/check_zenodo_metadata.py
