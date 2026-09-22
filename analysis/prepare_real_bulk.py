from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import zipfile
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/real_bulk"
OUT = ROOT / "data/processed/real_bulk"
REF = ROOT / "data/processed/music_mc_extension/20260918T154429497808Z"
OUT.mkdir(parents=True, exist_ok=True)
assert not (OUT / "bulk_counts.tsv.gz").exists(), "Reuse the prepared inputs."

ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
with zipfile.ZipFile(RAW / "mmc7.xlsx") as z:
    strings = ["".join(n.itertext()) for n in ET.fromstring(z.read("xl/sharedStrings.xml"))]
    sheet = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
    rows = []
    for row in sheet.findall("x:sheetData/x:row", ns):
        values = [None] * 61
        for c in row.findall("x:c", ns):
            letters = re.match(r"[A-Z]+", c.attrib["r"])[0]
            col = 0
            for letter in letters:
                col = col * 26 + ord(letter) - 64
            value = c.find("x:v", ns)
            if value is not None:
                values[col - 1] = strings[int(value.text)] if c.attrib.get("t") == "s" else float(value.text)
        rows.append(values)
flow = pd.DataFrame(rows[2:], columns=rows[1]).rename(columns={rows[1][0]: "sample_id"})
flow = flow.dropna(subset=["sample_id"]).set_index("sample_id")
assert len(flow) == 12 and flow.index.is_unique and "DZQV" not in flow.index
flow.to_csv(OUT / "flow_table_S6.csv")
mapping = pd.DataFrame([
    ("B", "B cells", "B Naive+Memory", 44),
    ("T", "CD4 T cells|CD8 T cells", "T cells", 56),
    ("NK", "NK cells", "NK", 22),
    ("Monocytes", "CD14+ Monocytes|FCGR3A+ Monocytes", "Monocytes", 46),
], columns=["lineage", "reference_types", "flow_label", "table_S6_zero_based_column"])
for row in mapping.itertuples():
    assert rows[1][row.table_S6_zero_based_column] == row.flow_label
truth = pd.DataFrame({row.lineage: flow[row.flow_label].astype(float) for row in mapping.itertuples()})
coverage = truth.sum(axis=1)
assert coverage.between(0, 100).all()
truth = truth.div(coverage, axis=0)
truth["represented_PBMC_percent"] = coverage
truth["reported_29_type_sum_percent"] = flow.iloc[:, :29].astype(float).sum(axis=1)
assert np.allclose(truth[mapping.lineage].sum(axis=1), 1)
truth.to_csv(OUT / "flow_truth.csv")
mapping.to_csv(OUT / "cell_type_mapping.csv", index=False)

meta = pd.read_csv(RAW / "sra.sra.SRP125125.MD.gz", sep="\t")
meta = meta.loc[meta.sample_title.str.contains("_PBMC_", regex=False)].copy()
meta["sample_id"] = meta.sample_title.str.split("_").str[0]
assert len(meta) == 13 and meta.sample_id.is_unique
meta["included"] = meta.sample_id.isin(flow.index)
assert set(meta.loc[~meta.included, "sample_id"]) == {"DZQV"}
qc = pd.read_csv(RAW / "sra.recount_qc.SRP125125.MD.gz", sep="\t")
meta = meta.merge(qc[["external_id", "star.average_mapped_length", "star.uniquely_mapped_reads_%"]], on="external_id", validate="one_to_one")
meta[["sample_id", "sample_name", "external_id", "sample_title", "included", "star.average_mapped_length", "star.uniquely_mapped_reads_%"]].to_csv(OUT / "sample_mapping.csv", index=False)
selected = meta.loc[meta.included].set_index("sample_id").loc[truth.index]
lengths = selected.set_index("external_id")["star.average_mapped_length"]
assert (lengths > 0).all()
counts = pd.read_csv(RAW / "sra.gene_sums.SRP125125.G026.gz", sep="\t", comment="#", usecols=["gene_id", *lengths.index]).set_index("gene_id")
counts = np.rint(counts[lengths.index].div(lengths, axis=1)).astype("int64")
counts.index = counts.index.str.replace(r"\.\d+(?=_PAR_Y$|$)", "", regex=True)
assert counts.index.is_unique and counts.ge(0).all().all()
counts.columns = truth.index
genes = (REF / "genes.tsv").read_text().splitlines()
assert len(genes) == len(set(genes))
common = counts.index.intersection(genes)
assert len(common) > 25000
counts.to_csv(OUT / "bulk_counts.tsv.gz", sep="\t", index_label="gene_id")
refs = pd.read_csv(REF / "references.tsv", sep="\t")
assert len(refs) == 1344 and refs.reference_id.is_unique
ref_manifest = json.loads((ROOT / "results/music_mc_extension_input/20260918T154429497808Z/input_manifest.json").read_text())
assert set(mapping.reference_types.str.split("|").explode()) == set(ref_manifest["cell_types"])

base = "https://duffel.rail.bio/recount3/human/data_sources/sra"
sources = {
    "Monaco_2019.xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6367568/fullTextXML",
    "Monaco_supplement.zip": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6367568/supplementaryFiles",
    "GSE107011_Processed_data_TPM.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107011/suppl/GSE107011_Processed_data_TPM.txt.gz",
    "sra.gene_sums.SRP125125.G026.gz": f"{base}/gene_sums/25/SRP125125/sra.gene_sums.SRP125125.G026.gz",
}
for kind in ["sra", "recount_qc", "recount_project"]:
    name = f"sra.{kind}.SRP125125.MD.gz"
    sources[name] = f"{base}/metadata/25/SRP125125/{name}"
manifest = {
    "prepared_utc": datetime.now(timezone.utc).isoformat(),
    "accessions": ["GSE107011", "SRP125125", "GSE107019"],
    "bulk_assembly": "GRCh38 / GENCODE v26",
    "counts_conversion": "round(base_pair_coverage / star.average_mapped_length)",
    "conversion_source": "https://github.com/LieberInstitute/recount3/blob/master/R/transform_counts.R",
    "TPM_used_for_fitting": False,
    "bulk_genes": len(counts), "reference_genes": len(genes), "common_genes": len(common),
    "bulk_samples": len(truth), "reference_configurations": len(refs),
    "reference_input_manifest": str(ROOT / "results/music_mc_extension_input/20260918T154429497808Z/input_manifest.json"),
    "cell_types": ref_manifest["cell_types"],
    "flow_worksheet": "GSE107019 (Our S13 cohort)",
    "represented_PBMC_percent_range": [float(coverage.min()), float(coverage.max())],
    "sources": [{"file": name, "url": url, "downloaded_utc": datetime.fromtimestamp((RAW / name).stat().st_mtime, timezone.utc).isoformat(), "sha256": hashlib.sha256((RAW / name).read_bytes()).hexdigest()} for name, url in sources.items()],
    "extracted_tables": {name: hashlib.sha256((RAW / name).read_bytes()).hexdigest() for name in ["mmc2.xlsx", "mmc7.xlsx"]},
}
(OUT / "sources.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
print(json.dumps({k: manifest[k] for k in ["bulk_genes", "reference_genes", "common_genes", "bulk_samples", "reference_configurations", "represented_PBMC_percent_range"]}, indent=2))
print(truth.to_string())
