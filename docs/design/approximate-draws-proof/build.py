"""Concatenate texparts/p1..p4.tex into proof.tex, filling the placeholders from scripts/fill.txt
(lines 'NAME<TAB>text'; FLOORTABLE and WGRIDROWS from files) and run pdflatex twice.

Usage: python3 build.py [OUTDIR]. proof.tex, proof.pdf and the LaTeX auxiliaries go to OUTDIR (default: the
current directory), never beside this script."""
import subprocess, os, sys
root = os.path.dirname(os.path.abspath(__file__))
out = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.getcwd())
os.makedirs(out, exist_ok=True)
s = ''.join(open(f'{root}/texparts/p{i}.tex').read() for i in (1, 2, 3, 4))
fill = {}
for line in open(f'{root}/scripts/fill.txt'):
    if '\t' in line:
        k, v = line.rstrip('\n').split('\t', 1)
        fill[k] = v
for k in sorted(fill, key=len, reverse=True):
    v = fill[k]
    if v.startswith('@file:'):
        v = open(f'{root}/scripts/' + v[6:]).read()
    s = s.replace('\\' + k + '{}', v).replace('\\' + k, v).replace(k, v)
open(f'{out}/proof.tex', 'w').write(s)
env = dict(os.environ)
for _ in range(2):
    r = subprocess.run(['pdflatex', '-interaction=nonstopmode', 'proof.tex'], cwd=out, env=env,
                       capture_output=True, text=True)
log = open(f'{out}/proof.log', errors='replace').read()
import re
for m in re.finditer(r'^(! .*|LaTeX Warning: (Reference|Citation|Label).*|.*undefined.*|Overfull.*)$', log, re.M):
    print(m.group(0))
print('pages:', re.findall(r'Output written on proof.pdf \((\d+) pages', log))
