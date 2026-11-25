"""Build the anonymous source archive, then compile clean copies with both engines."""
from pathlib import Path
import argparse
import json
import os
import re
import shutil
import subprocess
import tempfile
from zipfile import ZipFile, ZIP_DEFLATED

ROOT = Path(__file__).resolve().parents[1]
WARNING = re.compile(r'(?:LaTeX|Package \S+|Font) Warning|Overfull |Underfull |undefined (?:references|citations)|^warning:', re.I | re.M)

def sources():
    names = {Path('iopjournal.cls')}
    def visit(name):
        if name in names:
            return
        names.add(name)
        source = (ROOT/name).read_text(encoding='ascii')
        for dep in re.findall(r'\\input\{([^}]+)\}', source):
            visit(Path(dep).with_suffix('.tex'))
        for fig in re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}',source):
            names.add(Path('figures')/Path(fig).with_suffix('.pdf'))
    visit(Path('manuscript_ip.tex'))
    return sorted(names)

def run(args,where):
    proc=subprocess.run(args, cwd=where,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    if proc.returncode:
        raise RuntimeError(proc.stdout[-7000:])
    return proc.stdout

def main():
    pr=argparse.ArgumentParser()
    pr.add_argument('--pdflatex',default=shutil.which('pdflatex'))
    pr.add_argument('--tectonic',default=shutil.which('tectonic'))
    args=pr.parse_args()
    if not args.pdflatex or not args.tectonic:
        pr.error('provide paths to pdflatex and tectonic')
    names=sources()
    out=ROOT/'IP_Anonymous_LaTeX.zip'
    with ZipFile(out,'w',ZIP_DEFLATED,compresslevel=5) as z:
        for name in names:
            z.write(ROOT/name,name)
    report=dict(source_files=[n.as_posix() for n in names],engines={})
    with tempfile.TemporaryDirectory(prefix='material-sets-submit-') as base:
        for engine,path in [('pdflatex',args.pdflatex),('tectonic',args.tectonic)]:
            d=Path(base)/engine
            d.mkdir()
            with ZipFile(out) as z:
                z.extractall(d)
            if engine=='pdflatex':
                cmd=[path,'-interaction=nonstopmode','-halt-on-error','manuscript_ip.tex']
                for _ in range(3): result=run(cmd,d)
            else:
                result=run([path,'--keep-logs','manuscript_ip.tex'],d)
            log=(d/'manuscript_ip.log').read_text(errors='replace')
            warnings=WARNING.findall(log+'\n'+result)
            if warnings:
                raise RuntimeError(engine+' warnings: '+str(warnings[:12])+'\n'+'\n'.join(l for l in log.splitlines() if WARNING.search(l))[:3000])
            pdf=d/'manuscript_ip.pdf'
            meta=run(['pdfinfo',str(pdf)],d)
            pages=int(re.search(r'^Pages:\s*(\d+)',meta,re.M).group(1))
            if not 20<=pages<=24:
                raise ValueError('unexpected page count '+str(pages))
            author=re.search(r'^Author:[ \t]*([^\n]*)',meta,re.M)
            if author and author.group(1).strip():raise ValueError('nonempty author metadata')
            rendered=run(['pdftotext','-layout',str(pdf),'-'],d)
            if '??' in rendered:raise ValueError('unresolved PDF reference')
            if engine=='pdflatex':
                shutil.copy2(pdf,ROOT/'manuscript_ip.pdf')
            report['engines'][engine]=dict(pages=pages,references_resolved=True,layout_warnings=warnings,source='empty directory, source archive only')
    report['ascii_tex_and_class']=all((ROOT/n).read_bytes().isascii() for n in names if n.suffix in {'.tex','.cls'})
    if not report['ascii_tex_and_class']:raise ValueError('non-ASCII source')
    (ROOT/'results/typesetting.json').write_text(json.dumps(report,indent=2)+'\n')
    print(report['engines'])

if __name__=='__main__':main()
