"""Fetch only identified processed GSE96583 batch-2 control inputs."""
from pathlib import Path
from datetime import datetime, timezone
from urllib.request import Request, urlopen
import argparse
import json
import shutil
from run_deconvolution_pilot import sha, save_json

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('series_audit', type=Path)
parser.add_argument('sample_audit', type=Path)
args = parser.parse_args()
series = json.loads(args.series_audit.read_bytes())
samples = json.loads(args.sample_audit.read_bytes())
series_record = next(r for r in series['records'] if r['accession'] == 'GSE96583')
control = next(r for r in samples['records'] if r['accession'] == 'GSM2560248')
assert 'stimulated with: none (control)' in control['parsed']['characteristics_ch1']
urls = [u for u in series_record['parsed']['supplementary_files'] if '_batch2.' in u]
urls += control['parsed']['supplementary_files']
assert len(urls) == 4
run = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
raw = ROOT / 'data/raw/GSE96583' / run
out = ROOT / 'results/kang_download' / run
raw.mkdir(parents=True); out.mkdir(parents=True)
if shutil.disk_usage(ROOT).free < 5 * 1024**3:
    raise ValueError('Less than 5 GiB free')
files = []
for source in urls:
    url = source.replace('ftp://ftp.ncbi.nlm.nih.gov/', 'https://ftp.ncbi.nlm.nih.gov/', 1)
    assert url.startswith('https://ftp.ncbi.nlm.nih.gov/geo/')
    with urlopen(Request(url, method='HEAD'), timeout=30) as response:
        expected = int(response.headers.get('Content-Length', 0))
    if expected > 512 * 1024**2:
        raise ValueError('Unexpectedly large processed input')
    path = raw / url.rsplit('/', 1)[1]
    received = 0
    with urlopen(url, timeout=60) as response, path.open('xb') as handle:
        final_url = response.geturl()
        while chunk := response.read(1024**2):
            received += len(chunk)
            if received > 512 * 1024**2:
                raise ValueError('Exceeded download ceiling')
            handle.write(chunk)
    if expected and expected != received:
        raise ValueError('Download length mismatch')
    files.append({'deposited_url': source, 'url': url, 'final_url': final_url,
                  'retrieved_utc': datetime.now(timezone.utc).isoformat(),
                  'path': path.relative_to(ROOT).as_posix(), 'bytes': received, 'sha256': sha(path)})
    print(f'{path.name}: {received} bytes', flush=True)
save_json(out / 'download.json', {'dataset': 'GSE96583 batch2 control GSM2560248',
          'script': Path(__file__).relative_to(ROOT).as_posix(), 'script_sha256': sha(__file__),
          'series_audit': args.series_audit.as_posix(), 'series_audit_sha256': sha(args.series_audit),
          'sample_audit': args.sample_audit.as_posix(), 'sample_audit_sha256': sha(args.sample_audit),
          'files': files})
print('DOWNLOAD_MANIFEST=' + str(out / 'download.json'), flush=True)
