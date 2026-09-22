from pathlib import Path
import hashlib
import json
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PREP = ROOT / "data/processed/real_bulk_salmon"
OUT = ROOT / "results/real_bulk_salmon"
assert not (PREP / "bulk_counts.tsv").exists(), "Reuse the prepared counts."
quant = pd.read_csv(OUT / "salmon/quant.sf", sep="\t")
mapping = pd.read_csv(PREP / "tx2gene.tsv", sep="\t")
assert quant.Name.is_unique and mapping.transcript_id.is_unique
assert set(quant.Name).issubset(set(mapping.transcript_id))
missing = set(mapping.transcript_id) - set(quant.Name)
index_log = Path(r"\\wsl.localhost\PBMC-Ubuntu\opt\pbmc\work\real_bulk_salmon\index\pre_indexing.log")
short = set(re.findall(r"Entry with header \[([^]]+)\], had length less", (index_log.read_text() + index_log.with_name("ref_indexing.log").read_text())))
assert missing.issubset(short), missing
assert quant.loc[quant.Name.isin(short), "NumReads"].eq(0).all()

assert np.isfinite(quant[["TPM", "NumReads"]]).all().all()
assert quant.NumReads.ge(0).all()
joined = quant.merge(mapping, left_on="Name", right_on="transcript_id", validate="one_to_one")
counts = joined.groupby("gene_id").NumReads.sum().reindex(sorted(mapping.gene_id.unique()), fill_value=0.0)
assert np.isclose(counts.sum(), quant.NumReads.sum(), rtol=1e-12)
counts.rename("NumReads").to_csv(PREP / "gene_estimated_counts.tsv", sep="\t")
rounded = pd.DataFrame({"453W": np.rint(counts).astype("int64")})
rounded.to_csv(PREP / "bulk_counts.tsv", sep="\t")
meta = json.loads((OUT / "salmon/aux_info/meta_info.json").read_text())
qc = json.loads((OUT / "fastp.json").read_text())
gene_symbols = mapping.drop_duplicates("gene_id").set_index("gene_id").gene_name
is_rpl = gene_symbols.reindex(counts.index).str.match(r"^RP[LS]")
result = {
    "sample_id": "453W", "run": "SRR6298350", "transcripts": len(quant), "genes": len(counts),
    "transcripts_excluded_by_indexer_as_short": sorted(short.intersection(mapping.transcript_id)),
    "annotated_transcripts_absent_from_quant_table": sorted(missing),
    "positive_genes": int(rounded["453W"].gt(0).sum()), "estimated_fragments": float(counts.sum()),
    "rounded_fragments": int(rounded["453W"].sum()),
    "read_retention": qc["summary"]["after_filtering"]["total_reads"] / qc["summary"]["before_filtering"]["total_reads"],
    "mapping_percent": meta["percent_mapped"], "processed_fragments": meta["num_processed"],
    "mapped_fragments": meta["num_mapped"], "library_types": meta["library_types"],
    "rpl_rps_estimated_count_fraction_all_annotated_genes": float(counts[is_rpl].sum() / counts.sum()),
    "bulk_counts_sha256": hashlib.sha256((PREP / "bulk_counts.tsv").read_bytes()).hexdigest(),
}
(OUT / "input_quality.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
