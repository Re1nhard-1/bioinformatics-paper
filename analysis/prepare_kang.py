"""Parse original GSE96583 control counts and deposited annotations, before prediction."""
from __future__ import annotations
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
from pathlib import Path
from datetime import datetime, timezone
from itertools import combinations
import argparse
import gzip
import json
import numpy as np
import pandas as pd
from scipy import sparse, io
from run_deconvolution_pilot import sha, save_json

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('download_manifest', type=Path)
args = parser.parse_args()
download = json.loads(args.download_manifest.read_bytes())
files = {}
for item in download['files']:
    path = ROOT / item['path']
    if sha(path) != item['sha256']:
        raise ValueError('Raw file hash mismatch')
    files[path.name] = path
annotations = pd.read_csv(files['GSE96583_batch2.total.tsne.df.tsv.gz'], sep='\t', index_col=0,
                          dtype={'ind': str})
annotations.index.name = 'cell_id'
control = annotations[annotations.stim.eq('ctrl')].copy()
if control.index.duplicated().any():
    raise ValueError('Duplicate control barcode')
barcodes = pd.read_csv(files['GSM2560248_barcodes.tsv.gz'], header=None, sep='\t')[0].to_numpy()
genes = pd.read_csv(files['GSE96583_batch2.genes.tsv.gz'], header=None, sep='\t', names=['ensembl_id', 'symbol'])
if len(set(barcodes)) != len(barcodes) or set(control.index) != set(barcodes):
    raise ValueError('Barcode mapping mismatch')
if genes.ensembl_id.duplicated().any() or genes.ensembl_id.isna().any():
    raise ValueError('Invalid Ensembl feature identifiers')
with gzip.open(files['GSM2560248_2.1.mtx.gz'], 'rb') as handle:
    original = io.mmread(handle).tocoo()
if original.shape != (len(genes), len(barcodes)):
    raise ValueError('Matrix dimensions disagree with annotations')
if not np.isfinite(original.data).all() or np.any(original.data < 0) or not np.equal(original.data, np.floor(original.data)).all():
    raise ValueError('Not finite, nonnegative, integer counts')
if original.data.max() > np.iinfo(np.int32).max:
    raise ValueError('Counts exceed int32')
matrix = original.astype(np.int32).T.tocsr()
if matrix.nnz != original.nnz:
    raise ValueError('Repeated matrix coordinates')
control = control.loc[barcodes].copy()
control['matrix_row'] = np.arange(len(control))
control['umi_total'] = np.asarray(matrix.sum(axis=1)).ravel()
if control.umi_total.min() <= 0:
    raise ValueError('Zero library in deposited control data')
included = control.multiplets.eq('singlet') & control.cell.notna() & control.ind.notna()
singlets = control[included].copy()
inventory = pd.crosstab(singlets.ind, singlets.cell)
# A deterministic inventory rule, evaluated before any deconvolution outputs.
# Include the largest complete >=50-cell rectangle, requiring >=4 donors and >=3 types.
candidates = []
donors = sorted(inventory.index)
for n in range(4, len(donors) + 1):
    for subset in combinations(donors, n):
        eligible = sorted(inventory.columns[(inventory.loc[list(subset)] >= 50).all(axis=0)])
        if len(eligible) >= 3:
            candidates.append((n * len(eligible), n, tuple(subset), tuple(eligible)))
if not candidates:
    raise ValueError('No eligible four-donor, three-type design')
best = sorted(candidates, key=lambda x: (-x[0], -x[1], x[2], x[3]))[0]
selected_donors, selected_types = list(best[2]), list(best[3])
run = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
processed = ROOT / 'data/processed/kang' / run
out = ROOT / 'results/kang_input' / run
processed.mkdir(parents=True); out.mkdir(parents=True)
genes_path = processed / 'genes.json'
save_json(genes_path, genes.ensembl_id.tolist())
genes.to_csv(out / 'gene_annotations.csv', index=False, lineterminator='\n')
inventory.to_csv(out / 'all_singlet_donor_type_inventory.csv', lineterminator='\n')
control.reset_index().to_csv(out / 'control_cell_audit.csv.gz', index=False,
                           compression={'method': 'gzip', 'mtime': 0}, lineterminator='\n')
records = []
for donor in donors:
    part = singlets[singlets.ind.eq(donor)].copy()
    counts = matrix[part.matrix_row.to_numpy()]
    frame = pd.DataFrame({'donor': donor, 'cell_id': part.index, 'cell_type': part.cell.to_numpy(),
                          'umi_total': part.umi_total.to_numpy(), 'original_matrix_row': part.matrix_row.to_numpy()})
    matpath, obspath = processed / f'{donor}_counts.npz', processed / f'{donor}_cells.csv'
    sparse.save_npz(matpath, counts)
    frame.to_csv(obspath, index=False, lineterminator='\n')
    records.append({'donor': donor, 'matrix_path': matpath.relative_to(ROOT).as_posix(),
                     'matrix_sha256': sha(matpath), 'cells_path': obspath.relative_to(ROOT).as_posix(),
                     'cells_sha256': sha(obspath), 'shape': list(counts.shape),
                     'cell_type_counts': frame.cell_type.value_counts().to_dict(),
                     'selected_for_external_design': donor in selected_donors})
report = {'run_id': run, 'dataset': 'GSE96583 batch2 unstimulated six-hour culture; GSM2560248',
          'script': Path(__file__).relative_to(ROOT).as_posix(), 'script_sha256': sha(__file__),
          'download_manifest': args.download_manifest.as_posix(), 'download_manifest_sha256': sha(args.download_manifest),
          'genes_path': genes_path.relative_to(ROOT).as_posix(), 'genes_sha256': sha(genes_path),
          'counts_shape_genes_by_cells': list(original.shape), 'matrix_nonzeros': int(original.nnz),
          'annotation_rows_all_conditions': len(annotations), 'control_rows': len(control),
          'control_multiplet_categories': control.multiplets.value_counts().to_dict(),
          'control_singlets_retained': len(singlets), 'control_annotation_barcode_match': True,
          'gene_symbols_duplicated': int(genes.symbol.duplicated().sum()),
          'selected_donors': selected_donors, 'selected_cell_types': selected_types,
          'selection_rule': 'Maximize donors*types for complete >=50-cell rectangle, >=4 donors and >=3 types; tie prefer more donors then lexicographic IDs.',
          'selection_uses_prediction_results': False, 'records': records,
          'limitations': ['Eight SLE source donors, not healthy controls; six-hour unstimulated culture.',
            'Deposited GEO singlet/type labels used as operational labels, not newly verified identities.',
            'Donor subset chosen using availability before predictions; conclusions conditional on this subset.',
            'Ensembl IDs preserve original rows; duplicate gene symbols are not silently collapsed.',
            'Not a real-bulk validation or a causal disease comparison.'],
          'artifacts': [{'path': p.relative_to(ROOT).as_posix(), 'sha256': sha(p)} for p in sorted(out.iterdir())]}
save_json(out / 'input_manifest.json', report)
print(json.dumps({'output': str(out), 'control_singlets': len(singlets),
                  'selected_donors': selected_donors, 'selected_types': selected_types,
                  'eligible_inventory': inventory.loc[selected_donors, selected_types].to_dict()}, indent=2))
