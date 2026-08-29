# Data integrity hashes

This file did not previously exist. CLAUDE.md and other capstone docs referenced it as
though it already held a hash table for the original 14 processed files; a full
`git log --all --diff-filter=A` search confirmed it was never committed under this or any
other name. Created now, computed directly from the working tree — not from memory.

**Pinned to commit `9381c5e5f0f1aa85463a4e32251e63de0a8ea16b`.** SHA-256, computed by reading
each file directly (`hashlib.sha256`, 1 MiB chunks). All 15 tracked files in
`data/processed/` (the 14 original files plus `analysis_results.json`).

```
12431dad60d07f0bd2bea9a680367007c9e030e9f17c5c20ef0b0694dcb548f9  analysis_results.json         5,428 bytes
b49169bd363a596f400b4faff8c21d354275b70404efe08b9109d38f1bdc0ffd  baseline_results.json         3,415 bytes
165906f07e61819c8fadb2bf3c95a73817e538a22a34f63c415a30222ac49b9f  crispr_effect.labels.json   334,743 bytes
9214efa3ce172079e6ce4ca78853d8bf92fb8f6d4a55d0c6c71e4653b59e8826  crispr_effect.npz        75,224,750 bytes
d18005cc0aec3e4d5f0fd06c748ef66672256cd8c7a6f24ea4c441b0ca785983  expression.labels.json      334,679 bytes
3d5bfa0c3430584f8943fd2365be0eecf8b994b38bfc7d491d59d7b9ff251a2d  expression.npz           65,301,390 bytes
a4b8069cc93af48f01e745bb1a15f4eaf4a7b67c9f92ca44bef3bb9e44c6d0a1  gene_columns.json           932,973 bytes
1e34deeddc539c524cfebe86f68ac419c8c5348a8b35db4832aec1a181d55302  gene_id_map.csv              508,111 bytes
af8ee6d734bea11101d07884f1c72d2b4efaff9875506738a037102a712f1e46  geneformer_embeddings.csv 10,387,948 bytes
1962206fa17646cbd1fec4b642a577cc2586c09c4cabd980541a7e11a8b6f894  head_results.json             3,312 bytes
bf65065b41564156bbc7d9b56957f7f49f7373c46f0b18a3f96bdbcde3894034  join_report.json              3,281 bytes
e6490f012a66b92953d0f2b4c0e0143589c0f2290d8b5a3d7fbf5c73808f0121  join_report.txt               5,841 bytes
1c314197b57c1f8363eb44f8902b3733777e7c304d7f677c76d401e3cabe5180  model_metadata.csv           166,502 bytes
68c8fe39ae8965ce20b04f50870609cc21734386ceeff859f4d0bddd2e5bab35  selective_genes.json          95,170 bytes
f1419abc7cbd31efc173a5857bab9eb318b53f8e535a17048bfcf0ea2f70aeef  splits.json                   33,095 bytes
```

To re-verify from a fresh clone at this commit:

```powershell
py -c "
import hashlib, os
for f in sorted(os.listdir('data/processed')):
    p = os.path.join('data/processed', f)
    if os.path.isfile(p):
        h = hashlib.sha256()
        with open(p, 'rb') as fh:
            for chunk in iter(lambda: fh.read(1<<20), b''):
                h.update(chunk)
        print(f'{h.hexdigest()}  {f}')
"
```

Compare each line against the table above. A mismatch means the checked-out file differs
from what this hash table describes — check `git log -- data/processed/<file>` before
assuming corruption; the file may have legitimately changed in a later commit, in which case
this table needs a new entry pinned to the new commit, not a silent overwrite of this one.

## Phase 2 additions (§9.7)

### Added in commit `d78fbf8d6f8dbf749b8d8e101876d16340a66ad1` (branch `side-v1`)

SHA-256 computed directly from the committed bytes (`hashlib.sha256`), not from memory.

```
652011a1cdb8ecf42812cc5fcd6c55947a77995ebe47893e3b307165467bb711  data/external/sid_osteosarc/BG003082.gene_tpm.gct.gz  778,119 bytes
c5c97ac03dd6bb4e4b9b1becea9c944e5e2742d4bc90af435ec3ddf6d66c1324  data/processed/ensembl_map.csv                       422,382 bytes
```

**`BG003082.gene_tpm.gct.gz`** — Sid Sijbrandij osteosarcoma primary-tumor RNA-seq, CC0 1.0.
Retrieved 2026-08-27 by anonymous read from
`s3://sid-sijbrandij-osteosarc-dataset/rna-seq/reprocessed/BG003082/BG003082.gene_tpm.gct.gz`
(region `us-west-2`). The S3 object's `ETag` `28485ea587fdd33cbde140f1a150a10a` equals the
file's MD5 (single-part upload); a fresh anonymous re-download is byte-identical to the
committed file. `gzip -t` passes; it decompresses to 2,897,595 bytes / 74,631 lines — a GCT
`#1.2` whose dimension line `74628\t1` matches its 74,628 data rows, columns
`Name`/`Description`/`BG003082`, `Name` entirely versioned Ensembl gene IDs, no duplicates,
TPM column non-negative and summing to exactly 1,000,000.

**`ensembl_map.csv`** — Entrez→Ensembl cache for `prepare_geneformer_input.load_ensembl_map`
(schema `entrez,ensembl_id`, both `str`, CRLF). Rebuilt from the static public reference
**NCBI `gene2ensembl.gz`** (`https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz`,
`Last-Modified` 2026-08-26, `Content-Length` 312,693,647,
SHA-256 `0030a333aeecf2151c57069bf3196ad9f89dc9ec0da5ee125b2256f9b068dbbf`, retrieved
2026-08-27), **not** the original `mygene`/Kaggle path — a deliberate provenance difference,
same Entrez join key. Method: filter `tax_id == 9606`; for each of the 18,460 canonical
Entrez IDs in `gene_columns.json`, resolve to its Ensembl gene; when an Entrez had more than
one Ensembl gene (18 cases), keep the one with the most transcript rows in `gene2ensembl`,
breaking further ties by smallest ENSG accession (deterministic, row-order independent);
rows written in `gene_columns.json` order. Coverage **18,459 / 18,460**. The one unmapped
gene is Entrez `79400` (NOX5): NCBI `gene2ensembl` and `Homo_sapiens.gene_info` carry no
Ensembl xref for it, and its true Ensembl ID (`ENSG00000255346`) was **not** patched in, to
keep the file single-provenance — so a local re-run of `prepare_geneformer_input` drops NOX5
(`genes_dropped_no_ensembl == 1`). The frozen Kaggle `geneformer_embeddings.csv` is
unaffected. Verified through the real consumer: `load_ensembl_map` →
`entrez_mapped=18459, entrez_unmapped=1`; `build_frames` →
`genes_mapped_to_ensembl=18459, genes_dropped_no_ensembl=1, ensembl_collisions_dropped=0`.

### Added with the BG003082 Geneformer embedding (Kaggle GPU step)

SHA-256 computed directly from the committed bytes (`hashlib.sha256`, 1 MiB chunks), not
from memory. Both files were produced on Kaggle by `capstone/kaggle_bg003082_embedding.py`
(Tesla T4; `geneformer` pinned to `04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5`) and committed
to `data/processed/` **unmodified** — not regenerated, rounded, or normalised locally.

```
06a4ab9f85e5ac908975268ed502912317503ed277d28eeab1663d8305835080  data/processed/geneformer_bg003082_embedding.csv             18,822 bytes
eb111660731b34e328bc4b07693c405070e24aad3a6866d25c71fb9294e1194c  data/processed/geneformer_bg003082_embedding.provenance.json   6,997 bytes
```

**`geneformer_bg003082_embedding.csv`** — the BG003082 CLS embedding: one data row
`BG003082`, columns `0`–`767`, every value finite. Geneformer-V2-104M_CLcancer, 4,096
context, `special_token=True`, CLS pooling from layer −1 — the same settings as the frozen
1,140-row `geneformer_embeddings.csv`, but a **separate sidecar**. The frozen matrix is
byte-unchanged (`af8ee6d734bea11101d07884f1c72d2b4efaff9875506738a037102a712f1e46`,
verified before/after inside the provenance JSON and re-checked on disk). Token-id list
SHA-256 `f5c33af88f48e9ceba50c0d1b1975dbe20108be6dc5836b9232fda2766b1c2ef` (token length
4,096; starts `<cls>`, ends `<eos>`; every id in vocab; independent norm/rank replication of
the top-50 tokens passed on Kaggle).

**`geneformer_bg003082_embedding.provenance.json`** — environment (Python 3.12.13,
`torch 2.10.0+cu128`, `transformers 4.49.0`, Tesla T4), the pinned and run-time-resolved
Geneformer revision (both `04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5`), the input-build
reconciliation, the tokenisation hard checks, and the frozen-matrix before/after hashes.
Input hashes recorded there reconcile with the committed repo inputs:
`BG003082.gene_tpm.gct.gz` `652011a1…`, `ensembl_map.csv` `c5c97ac0…`, `gene_columns.json`
`a4b80690…`; mapping counts 18,460 canonical → 18,427 mapped + 33 unresolved (NOX5 /
Entrez `79400` among them), 0 canonical-identifier collisions, no symbol fallback.

**Caveats unchanged by this artifact existing.** The embedding is bulk primary-tumour TPM
through a model trained on cultured cell lines (domain shift); the pseudo-count basis is
linear TPM (feasibility doc §4); the Entrez→Ensembl map is NCBI `gene2ensembl`, not the
`mygene` map that produced the training embeddings and is gone; the pinned Geneformer
revision is a fresh pin for this run, **not** the unrecorded Phase 1 training revision. The
embedding is **not** reproducible from the public repository (needs `geneformer` + a GPU)
and is **not** bytewise-commensurable with the Phase 1 embeddings. See
`capstone/geneformer-bg003082-feasibility.md` §3.

The committed DGIdb snapshot will be appended here, pinned to the commit that adds it, once
that file exists. Not yet computed — do not invent placeholder hashes.
