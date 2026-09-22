from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import gzip
import hashlib
import json
import re
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/real_bulk_salmon"
PREP = ROOT / "data/processed/real_bulk_salmon"
RAW.mkdir(parents=True, exist_ok=True)
PREP.mkdir(parents=True, exist_ok=True)
BASE = "https://ftp.ensembl.org/pub/grch37/release-87/"
SOURCES = [
    ("genome.fa.gz", BASE + "fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz", 869923173, None),
    ("annotation.gtf.gz", BASE + "gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz", 42400095, None),
    ("SRR6298350_1.fastq.gz", "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR629/000/SRR6298350/SRR6298350_1.fastq.gz", 845999643, "ec687818e372a47e48838563a56b7b1d"),
    ("SRR6298350_2.fastq.gz", "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR629/000/SRR6298350/SRR6298350_2.fastq.gz", 866295103, "b47e08b1a4e8c8f3a4ec9e3bf9d006a6"),
]


def download(item):
    name, url, size, expected_md5 = item
    path = RAW / name
    receipt = RAW / (name + ".json")
    if path.exists() and receipt.exists():
        saved = json.loads(receipt.read_text())
        assert path.stat().st_size == saved["bytes"] == size
        print("Reusing", name, flush=True)
        return saved
    if path.exists():
        raise RuntimeError(f"Unverified complete file already exists: {path}")
    part = path.with_name(path.name + ".part")
    sha, md5 = hashlib.sha256(), hashlib.md5()
    start = part.stat().st_size if part.exists() else 0
    if start:
        with part.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                sha.update(block)
                md5.update(block)
    req = urllib.request.Request(url, headers={"Range": f"bytes={start}-"} if start else {})
    with urllib.request.urlopen(req, timeout=120) as response:
        if start:
            assert response.status == 206 and response.headers["Content-Range"].startswith(f"bytes {start}-")
        with part.open("ab" if start else "wb") as handle:
            for block in iter(lambda: response.read(8 * 1024 * 1024), b""):
                handle.write(block)
                sha.update(block)
                md5.update(block)
    assert part.stat().st_size == size, name
    assert expected_md5 is None or md5.hexdigest() == expected_md5, name
    part.rename(path)
    saved = {"file": str(path.relative_to(ROOT)), "url": url, "bytes": size,
             "md5": md5.hexdigest(), "sha256": sha.hexdigest(), "ena_md5": expected_md5,
             "downloaded_utc": datetime.now(timezone.utc).isoformat()}
    receipt.write_text(json.dumps(saved, indent=2) + "\n", encoding="utf-8")
    print("Downloaded and checked", name, size, flush=True)
    return saved


def annotation():
    target = PREP / "hg19_ensembl87_10x3.filtered.gtf"
    if target.exists():
        print("Reusing prepared annotation", flush=True)
        return
    accepted = {"protein_coding", "lincRNA", "antisense"}
    transcripts = {}
    genes = {}
    exons = set()
    with gzip.open(RAW / "annotation.gtf.gz", "rt") as source, target.open("w", newline="\n") as output:
        for line in source:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            attrs = dict(re.findall(r'(\S+) "([^"]*)";', fields[8]))
            if attrs.get("gene_biotype") not in accepted:
                continue
            output.write(line)
            gene = attrs["gene_id"]
            genes[gene] = attrs.get("gene_name", gene)
            if fields[2] == "exon":
                transcript = attrs["transcript_id"]
                assert transcript not in transcripts or transcripts[transcript] == gene
                transcripts[transcript] = gene
                exons.add((fields[0], int(fields[3]) - 1, int(fields[4])))
    with (PREP / "tx2gene.tsv").open("w", newline="\n") as handle:
        handle.write("transcript_id\tgene_id\tgene_name\n")
        for transcript, gene in sorted(transcripts.items()):
            handle.write(f"{transcript}\t{gene}\t{genes[gene]}\n")
    with (PREP / "exons.bed").open("w", newline="\n") as handle:
        for chrom, start, end in sorted(exons):
            handle.write(f"{chrom}\t{start}\t{end}\n")
    result = {"assembly": "GRCh37", "ensembl_release": 87, "biotypes": sorted(accepted),
              "transcripts": len(transcripts), "genes": len(genes), "exons": len(exons),
              "build_rule_source": "https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-ref-build-steps",
              "build_rule": "Human reference 3.0.0 hg19; rebuild FASTA/GTF contents, not a downloaded Cell Ranger reference package",
              "bed_coordinates": "zero-based half-open; GTF start minus one"}
    (PREP / "annotation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(download, SOURCES))
    annotation()
