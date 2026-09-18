"""Verify a sealed package, without running analyses or accessing any network."""
from common import PACKAGE,read,verify_manifest
if __name__=='__main__':
 verify_manifest();print('Passed',len(read(PACKAGE/'MANIFEST.json')['files']),'file hashes')
