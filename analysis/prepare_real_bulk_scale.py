from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/real_bulk"
OUT = ROOT / "data/processed/real_bulk_scale"
OUT.mkdir(parents=True, exist_ok=True)
assert not (OUT / "sources.json").exists(), "Reuse prepared inputs."
bulk = pd.read_csv(ROOT / "data/processed/real_bulk/bulk_counts.tsv.gz", sep="\t", index_col=0)
gtf_path = RAW / "human.gene_sums.G026.gtf.gz"
gtf = pd.read_csv(gtf_path, sep="\t", comment="#", header=None)
assert gtf[2].eq("gene").all()
gene_ids = gtf[8].str.extract(r'gene_id "([^"]+)"')[0].str.replace(r"\.\d+(?=_PAR_Y$|$)", "", regex=True)
lengths = pd.Series(gtf[5].to_numpy(dtype=float), index=gene_ids, name="exonic_bp")
assert lengths.index.is_unique and set(bulk.index) == set(lengths.index)
lengths = lengths.loc[bulk.index]
assert np.isfinite(lengths).all() and lengths.gt(0).all()
lengths.to_csv(OUT / "gene_lengths.csv", index_label="gene_id")
bulk.div(lengths / 1000, axis=0).to_csv(OUT / "length_scaled.tsv.gz", sep="\t", index_label="gene_id")
tpm_path = RAW / "GSE107011_Processed_data_TPM.txt.gz"
tpm = pd.read_csv(tpm_path, sep="\t", index_col=0)
tpm = tpm[[f"{sample}_PBMC" for sample in bulk.columns]]
tpm.columns = bulk.columns
tpm.index = tpm.index.str.replace(r"\.\d+(?=_PAR_Y$|$)", "", regex=True)
assert tpm.index.is_unique and np.isfinite(tpm).all().all() and tpm.ge(0).all().all()
tpm.to_csv(OUT / "author_tpm.tsv.gz", sep="\t", index_label="gene_id")
source = json.loads((ROOT / "data/processed/real_bulk/sources.json").read_text())
manifest = json.loads(Path(source["reference_input_manifest"]).read_text())
refs = pd.read_csv(ROOT / manifest["reference_file"], sep="\t")
refs = refs.loc[refs.block.eq(3) & refs.budget.eq(300) & refs.level.eq("balanced")]
assert len(refs) == 14 and refs.triple_key.is_unique
refs.to_csv(OUT / "references.tsv", sep="\t", index=False)
sources = {
    "prepared_utc": datetime.now(timezone.utc).isoformat(),
    "original_sources": "data/processed/real_bulk/sources.json",
    "reference_input_manifest": source["reference_input_manifest"],
    "annotation": {"url": "https://duffel.rail.bio/recount3/human/annotations/gene_sums/human.gene_sums.G026.gtf.gz", "local_mtime_utc": datetime.fromtimestamp(gtf_path.stat().st_mtime, timezone.utc).isoformat(), "sha256": hashlib.sha256(gtf_path.read_bytes()).hexdigest(), "length_field": "GTF score: gene exonic base-pair length"},
    "author_tpm": {"url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107011/suppl/GSE107011_Processed_data_TPM.txt.gz", "sha256": hashlib.sha256(tpm_path.read_bytes()).hexdigest()},
    "bulk_genes": {"length_scaled": len(bulk), "author_tpm": len(tpm)},
    "reference_ids": refs.reference_id.tolist(),
    "samples": bulk.columns.tolist(),
    "transform": "length_scaled = saved read-equivalent counts / (exonic_bp / 1000); author_tpm = released PBMC TPM; raw sc UMI unchanged",
}
(OUT / "sources.json").write_text(json.dumps(sources, indent=2) + "\n")
print(json.dumps({"references": len(refs), "samples": len(bulk.columns), "genes": sources["bulk_genes"]}))
