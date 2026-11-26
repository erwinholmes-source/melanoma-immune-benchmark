#!/usr/bin/env python3
"""Download original public processed sources to a private directory and verify hashes."""
import argparse,hashlib,json,urllib.request
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--destination',type=Path,required=True);a=p.parse_args();a.destination.mkdir(parents=True,exist_ok=True)
manifest=Path(__file__).resolve().parents[1]/'resources/source_files.json'
for r in json.loads(manifest.read_text()):
 f=a.destination/r['file']
 if not f.exists():
  temp=f.with_name(f.name+'.partial')
  with urllib.request.urlopen(r['url'],timeout=180) as src,temp.open('wb') as dst:
   while block:=src.read(1024*1024):dst.write(block)
  temp.replace(f)
 h=hashlib.sha256(f.read_bytes()).hexdigest()
 if h!=r['sha256']:raise RuntimeError('Checksum mismatch for '+r['file']+'; retain the file separately and obtain the documented source version. No substitute was silently accepted.')
 print(r['file'],'verified')
