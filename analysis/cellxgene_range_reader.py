"""Read-only, identity-checked HTTP ranges with an auditable disk block cache.

Uses h5py's documented binary file-like driver. Close the HDF5 handle before
closing this object. No service threads; no whole-file download or full-file
hash is claimed. Each cached block has its own SHA256 and source byte range.
"""
import io,json,hashlib,time
from pathlib import Path
from urllib.request import Request,urlopen
from datetime import datetime,timezone
from collections import OrderedDict

class RangeReader(io.RawIOBase):
    def __init__(self,url,cache,*,expected_size,expected_etag,block_size=1048576,max_network_bytes=805306368,deadline_seconds=900):
        self.url=url;self.cache=Path(cache);self.cache.mkdir(parents=True,exist_ok=True)
        self.size=int(expected_size);self.etag=expected_etag;self.block_size=block_size
        self.max_network_bytes=max_network_bytes;self.started=time.monotonic();self.deadline_seconds=deadline_seconds
        self.position=0;self.network_bytes=0;self.memory=OrderedDict();self.used={};self.requests=0
        identity={'url':url,'size':self.size,'etag':self.etag,'block_size':self.block_size}
        ip=self.cache/'identity.json'
        if ip.exists():assert json.loads(ip.read_text())==identity,'Cache identity mismatch'
        else:ip.write_text(json.dumps(identity,indent=2)+'\n',encoding='utf-8')
        super().__init__()
    def readable(self):return True
    def seekable(self):return True
    def writable(self):return False
    def tell(self):return self.position
    def seek(self,offset,whence=0):
        base={0:0,1:self.position,2:self.size}[whence];pos=base+int(offset)
        if pos<0:raise ValueError('negative seek')
        self.position=pos;return pos
    def readinto(self,buffer):
        data=self.read(len(buffer));buffer[:len(data)]=data;return len(data)
    def _block(self,index):
        if index in self.memory:
            self.memory.move_to_end(index);return self.memory[index]
        start=index*self.block_size;end=min(self.size,start+self.block_size)-1
        bp=self.cache/f'block_{index:06d}.bin';mp=self.cache/f'block_{index:06d}.json'
        if bp.exists() and mp.exists():
            rec=json.loads(mp.read_text());data=bp.read_bytes()
            assert rec['start']==start and rec['end']==end and rec['etag']==self.etag
            assert hashlib.sha256(data).hexdigest()==rec['sha256'] and len(data)==end-start+1
        else:
            if self.network_bytes+end-start+1>self.max_network_bytes:raise RuntimeError('Bounded range-transfer byte budget exhausted')
            if time.monotonic()-self.started>self.deadline_seconds:raise TimeoutError('Bounded range-transfer deadline exceeded')
            request=Request(self.url,headers={'Range':f'bytes={start}-{end}','If-Match':self.etag,'Accept-Encoding':'identity','User-Agent':'bioinformatics-paper-range-audit/1.0'})
            with urlopen(request,timeout=30) as response:
                assert response.status==206,('Server did not honor range',response.status)
                assert response.headers.get('Content-Range')==f'bytes {start}-{end}/{self.size}'
                assert response.headers.get('ETag')==self.etag
                data=response.read(end-start+2)
                assert len(data)==end-start+1
            self.network_bytes+=len(data);self.requests+=1
            rec={'start':start,'end':end,'bytes':len(data),'etag':self.etag,'sha256':hashlib.sha256(data).hexdigest(),'retrieved_utc':datetime.now(timezone.utc).isoformat()}
            part=bp.with_suffix('.part');part.write_bytes(data);part.replace(bp)
            mp.write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8')
            if self.requests%25==0:print(f'Fetched {self.network_bytes/2**20:.1f} MiB in {self.requests} range requests',flush=True)
        self.used[index]=rec;self.memory[index]=data
        if len(self.memory)>32:self.memory.popitem(last=False)
        return data
    def read(self,size=-1):
        size=self.size-self.position if size<0 else min(size,self.size-self.position)
        if size<=0:return b''
        if size>256*2**20:raise RuntimeError('Refusing one oversized whole-file-like read')
        chunks=[];remaining=size
        while remaining:
            index=self.position//self.block_size;offset=self.position%self.block_size
            block=self._block(index);take=min(remaining,len(block)-offset)
            if take<=0:raise EOFError('Unexpected block boundary')
            chunks.append(block[offset:offset+take]);self.position+=take;remaining-=take
        return b''.join(chunks)
    def report(self):
        return {'url':self.url,'file_bytes':self.size,'etag':self.etag,'block_size':self.block_size,'network_bytes_this_run':self.network_bytes,'range_requests_this_run':self.requests,'used_unique_blocks':len(self.used),'used_bytes':sum(r['bytes'] for r in self.used.values()),'elapsed_seconds':time.monotonic()-self.started,'blocks':[dict(index=k,**v) for k,v in sorted(self.used.items())],'whole_file_downloaded':False,'whole_file_sha256':None}
