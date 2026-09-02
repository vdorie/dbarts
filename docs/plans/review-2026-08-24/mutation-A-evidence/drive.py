import json,os,subprocess,sys,time,shutil
W="<private-scratch-prefix>/mutA-37158"
SRC="/Users/vdorie/Repositories/dbarts"
BUILD=W+"/build"; PRIST=W+"/pristine"; LIB=W+"/lib"; LOGS=W+"/logs"

def touch_sources():
    for f in os.listdir(BUILD+"/src"):
        if f.endswith(".cpp") or f.endswith(".c"):
            os.utime(BUILD+"/src/"+f, None)

def restore(path):
    shutil.copy2(os.path.join(PRIST,path), os.path.join(BUILD,path))
    os.utime(os.path.join(BUILD,path), None)
    if path.startswith("src/"):
        touch_sources()

def apply_mut(m):
    p=os.path.join(BUILD,m["file"])
    t=open(p).read()
    n=t.count(m["find"])
    if n!=1:
        return "ANCHOR_%d"%n
    open(p,"w").write(t.replace(m["find"],m["repl"]))
    if m["file"].startswith("src/"):
        touch_sources()
    return None

def install():
    env=dict(os.environ); env["MAKEFLAGS"]="-j10"
    r=subprocess.run(["R","CMD","INSTALL","-l",LIB,"."],cwd=BUILD,capture_output=True,text=True,env=env)
    return r.returncode, r.stdout[-3000:]+r.stderr[-3000:]

def runtests(out):
    env=dict(os.environ); env["SRCROOT"]=SRC; env["R_LIBS"]=LIB
    try:
        r=subprocess.run(["Rscript",W+"/runfails.R",os.environ.get("LISTFILE",W+"/allscope.txt"),out],capture_output=True,text=True,env=env,timeout=int(os.environ.get("TESTTIMEOUT","420")))
    except subprocess.TimeoutExpired:
        return 99,"TIMEOUT"
    return r.returncode, r.stdout+r.stderr

def parse(path):
    d={}
    for line in open(path):
        p=line.rstrip("\n").split("\t")
        if len(p)>=3:
            d[p[0]]=(p[1],p[2],p[3] if len(p)>3 else "")
    return d

def main():
    muts=[json.loads(l) for l in open(sys.argv[1]) if l.strip() and not l.startswith("//")]
    only=sys.argv[2].split(",") if len(sys.argv)>2 else None
    base=parse(os.environ.get("BASEFILE",LOGS+"/base.txt"))
    res=[]
    for m in muts:
        if only and m["id"] not in only: continue
        t0=time.time()
        restore(m["file"])
        err=apply_mut(m)
        if err:
            print(m["id"],"ANCHOR-FAIL",err,m["file"]); restore(m["file"]); continue
        rc,log=install()
        if rc!=0:
            open(LOGS+"/%s-install.log"%m["id"],"w").write(log)
            print(m["id"],"BUILD-FAIL (see logs)"); restore(m["file"]); continue
        outp=LOGS+"/%s.txt"%m["id"]
        rct,_=runtests(outp)
        if rct==99:
            restore(m["file"]); print("%s\t%s\tHANG\tKILLERS=hang"%(m["id"],m.get("note","")[:60]),flush=True); continue
        cur=parse(outp)
        killers=[]
        if len(cur) < len(base):
            killers.append("SESSION-CRASH(ran %d of %d files)"%(len(cur),len(base)))
        for k,v in cur.items():
            b=base.get(k)
            if b is None: continue
            if v[2]!=b[2] or v[1]!=b[1]:
                killers.append("%s(%s%s)"%(k,v[1],(" "+v[2][:120]) if v[2] else ""))
        restore(m["file"])
        el=time.time()-t0
        line="%s\t%s\t%ds\tKILLERS=%d\t%s"%(m["id"],m.get("note","")[:60],int(el),len(killers),"; ".join(killers))
        print(line,flush=True)
        res.append(line)
    open(LOGS+"/summary-%s.txt"%time.strftime("%H%M%S"),"w").write("\n".join(res)+"\n")

main()
