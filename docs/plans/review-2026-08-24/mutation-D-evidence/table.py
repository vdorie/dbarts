import re,sys,os
rows=[]
tot={}
for line in open("static-table.txt"):
    line=line.rstrip()
    if not line or line.startswith("TOTAL"): continue
    mm=re.match(r"^(\S+)\s+(\d+)\s*(.*)$", line)
    if not mm: continue
    name,n,rest=mm.group(1),int(mm.group(2)),mm.group(3)
    d={}
    for k,v in re.findall(r"([A-Za-z_:]+)=(\d+)", rest): d[k]=int(v)
    g=lambda *ks: sum(d.get(k,0) for k in ks)
    rows.append((name.replace("test-","").replace(".R",""), n,
        g("content_pin"), g("shape_pin","shape"), g("expr_vs_expr"),
        g("all_true","all_finite"),
        g("bound","finite","bool_other","smoke_silent","other:expect_inherits"),
        g("refusal_pat","refusal_bare")))
for r in rows:
    for i,k in enumerate(["n","pin","shp","cmp","all","smk","ref"]):
        tot[k]=tot.get(k,0)+r[i+1]
import json
json.dump(rows, open("static-rows.json","w"))
print("%-52s %4s %3s %3s %3s %3s %3s %3s"%("file","n","pin","shp","cmp","all","smk","ref"))
for r in rows: print("%-52s %4d %3d %3d %3d %3d %3d %3d"%r)
print("TOTAL", tot)
