"""Small file helpers; package and input paths are explicit, never workstation-specific."""
from pathlib import Path
import hashlib,json
PACKAGE=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*2**20),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(Path(p).read_bytes())
def save(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
def verify(p,h):
 if sha(p)!=h:raise ValueError('SHA256 mismatch: '+Path(p).name)
def verify_manifest():
 path=PACKAGE/'MANIFEST.json'
 if not path.exists():raise FileNotFoundError('Package must be sealed before execution')
 for r in read(path)['files']:verify(PACKAGE/r['path'],r['sha256'])
