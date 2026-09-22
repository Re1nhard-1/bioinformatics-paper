set -euo pipefail
ulimit -c 0
project="$1"
raw="$project/data/raw/real_bulk_salmon"
prep="$project/data/processed/real_bulk_salmon"
out="$project/results/real_bulk_salmon"
work=/opt/pbmc/work/real_bulk_salmon
export PATH=/opt/pbmc/reference-tools/bin:/opt/pbmc/salmon/bin:$PATH
export LC_ALL=C
mkdir -p "$out" "$work"
cd "$work"
/opt/pbmc/bootstrap/bin/micromamba list --prefix /opt/pbmc/reference-tools --explicit > "$project/environment/salmon/reference-tools-explicit.txt"
if [ ! -s genome.fa ]; then
    gzip -dc "$raw/genome.fa.gz" > genome.fa.partial
    mv genome.fa.partial genome.fa
fi
if [ ! -s transcripts.fa ]; then
    gffread "$prep/hg19_ensembl87_10x3.filtered.gtf" -g genome.fa -w transcripts.fa.partial
    mv transcripts.fa.partial transcripts.fa
fi
if [ ! -s gentrome.fa ]; then
    awk '/^>/{print substr($1,2)}' genome.fa > decoys.txt
    cat transcripts.fa genome.fa > gentrome.fa.partial
    mv gentrome.fa.partial gentrome.fa
    sha256sum transcripts.fa gentrome.fa decoys.txt > "$prep/reference_sha256.txt"
fi
if [ ! -s index/info.json ]; then
    salmon index -t gentrome.fa -d decoys.txt -i index -k 23 -p 6 --keepDuplicates --sparse
fi
if [ "${2:-all}" = reference ]; then
    exit 0
fi
if [ ! -s "$out/fastp.json" ]; then
    fastp --in1 "$raw/SRR6298350_1.fastq.gz" --in2 "$raw/SRR6298350_2.fastq.gz" --out1 reads_1.fastq.gz --out2 reads_2.fastq.gz --detect_adapter_for_pe --length_required 25 --thread 6 --json "$out/fastp.json" --html "$out/fastp.html"
fi
if [ ! -s "$out/salmon/quant.sf" ]; then
    salmon quant -i index -l A -1 reads_1.fastq.gz -2 reads_2.fastq.gz --validateMappings --seqBias --gcBias -p 6 -o "$out/salmon"
fi
test -s "$out/salmon/aux_info/meta_info.json"
