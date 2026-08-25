import json,os,subprocess,sys,time,shutil
W=os.environ["W"]
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

def runtests(listfile,out,timeout):
    env=dict(os.environ); env["SRCROOT"]=SRC; env["R_LIBS"]=LIB
    try:
        r=subprocess.run(["Rscript",W+"/runfails.R",listfile,out],capture_output=True,text=True,env=env,timeout=timeout)
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
    only=sys.argv[2].split(",") if len(sys.argv)>2 and sys.argv[2] else None
    wide = os.environ.get("WIDE","0")=="1"
    base=parse(LOGS+"/base.txt")
    for m in muts:
        if only and m["id"] not in only: continue
        t0=time.time()
        restore(m["file"])
        err=apply_mut(m)
        if err:
            print(m["id"],"ANCHOR-FAIL",err,m["file"],flush=True); restore(m["file"]); continue
        rc,log=install()
        if rc!=0:
            open(LOGS+"/%s-install.log"%m["id"],"w").write(log)
            print(m["id"],"BUILD-FAIL (see logs)",flush=True); restore(m["file"]); continue
        if wide:
            listfile=W+"/half-paths.txt"; nexp=51; to=420
        else:
            listfile=LOGS+"/%s-list.txt"%m["id"]
            open(listfile,"w").write("\n".join("inst/tinytest/"+f for f in m["files"])+"\n")
            nexp=len(m["files"]); to=300
        outp=LOGS+"/%s%s.txt"%(m["id"],"-wide" if wide else "")
        rct,_=runtests(listfile,outp,to)
        if rct==99:
            restore(m["file"]); print("%s\t%s\tHANG"%(m["id"],m.get("note","")[:60]),flush=True); continue
        cur=parse(outp) if os.path.exists(outp) else {}
        killers=[]
        if len(cur) < nexp:
            killers.append("SESSION-CRASH(ran %d of %d)"%(len(cur),nexp))
        for k,v in cur.items():
            b=base.get(k)
            if b is None: continue
            if v[2]!=b[2] or v[1]!=b[1]:
                killers.append("%s(%s%s)"%(k,v[1],(" "+v[2][:150]) if v[2] else ""))
        restore(m["file"])
        el=time.time()-t0
        print("%s\t%s\t%ds\tK=%d\t%s"%(m["id"],m.get("note","")[:60],int(el),len(killers),"; ".join(killers)),flush=True)

main()
