"""Collect/hash frozen processed inputs from a project copy, or verify a relocated bundle."""
import argparse,os,shutil
from pathlib import Path
from common import PACKAGE,read,save,verify,sha,verify_manifest
def main():
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
 q=sub.add_parser('collect');q.add_argument('--project-root',type=Path,required=True);q.add_argument('--output',type=Path,required=True);q.add_argument('--hardlink',action='store_true',help='Same-volume read-only validation only; shares file contents')
 q=sub.add_parser('verify');q.add_argument('--inputs',type=Path,required=True)
 a=p.parse_args();verify_manifest();plan=read(PACKAGE/'STUDY.json');records=[]
 if a.command=='collect':
  dest=a.output.resolve();dest.mkdir(parents=True,exist_ok=False)
  for cohort,v in plan['cohorts'].items():
   (dest/cohort).mkdir()
   for item in v['files']:
    src=a.project_root/item['project_source'];target=dest/cohort/item['filename'];verify(src,item['sha256'])
    if a.hardlink and item['filename'].endswith('.mtx.gz'):os.link(src,target)
    else:shutil.copy2(src,target)
    verify(target,item['sha256']);records.append(dict(path=target.relative_to(dest).as_posix(),bytes=item['bytes'],sha256=item['sha256'],mode='hardlink' if a.hardlink and item['filename'].endswith('.mtx.gz') else 'copy'))
  save(dest/'collected.json',dict(status='passed',files=records,source_package_manifest_sha256=sha(PACKAGE/'MANIFEST.json')))
 else:
  for cohort,v in plan['cohorts'].items():
   for item in v['files']:verify(a.inputs/cohort/item['filename'],item['sha256']);records.append(item)
 print('Verified',len(records),'input files;',sum(x['bytes'] for x in records),'bytes')
if __name__=='__main__':main()
