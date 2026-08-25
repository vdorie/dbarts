import re,glob,json,os
aim=json.load(open("aim.json"))
rows=json.load(open("static-rows.json"))
byfile={}; bymut={}
for r in sorted([x for x in glob.glob("logs/run*.out") if "snapshot" not in x and "base" not in x]):
    for line in open(r):
        if "\tKILLERS=" not in line: continue
        parts=line.rstrip("\n").split("\t")
        mid=parts[0]; det=parts[-1]
        n=int(re.search(r"KILLERS=(\d+)",line).group(1))
        files=set(x[len("test-"):-2] for x in re.findall(r"(test-[A-Za-z0-9._-]+\.R)\(", det))
        if "SESSION-CRASH" in det or "HANG" in det: files=set(["<all>"])
        bymut[mid]=(n,sorted(files),parts[1])
        for f in files: byfile.setdefault(f,set()).add(mid)
P={}
for mid,fs in aim.items():
    if mid in bymut:
        for f in fs: P[f]=P.get(f,0)+1
out=[]
print("%-52s %4s %3s %3s %3s %3s %3s %3s %3s %3s %3s  %s"%("file","n","pin","shp","cmp","all","smk","ref","P","M","Mn","verdict"))
ok=weak=gap=0
for r in rows:
    name=r[0]
    mset=byfile.get(name,set())
    narrow=set(k for k in mset if len(bymut[k][1])<=10)
    M=len(mset); Mn=len(narrow)
    p=P.get(name,0)
    v = "OK" if Mn>=2 else ("WEAK" if Mn==1 else "GAP")
    if v=="OK": ok+=1
    elif v=="WEAK": weak+=1
    else: gap+=1
    print("%-52s %4d %3d %3d %3d %3d %3d %3d %3d %3d %3d  %s"%(name,r[1],r[2],r[3],r[4],r[5],r[6],r[7],p,M,Mn,v))
print()
print("OK %d, WEAK %d, GAP %d"%(ok,weak,gap))
print("mutations scored:",len(bymut))
print()
print("ZERO-KILLER:")
for k,(n,files,note) in sorted(bymut.items()):
    if n==0: print("  %-5s %s"%(k,note))
