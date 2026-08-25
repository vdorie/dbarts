import re,sys,os
root="/Users/vdorie/Repositories/dbarts"

def split_args(s):
    out=[];d=0;q=None;cur=""
    for ch in s:
        if q:
            cur+=ch
            if ch==q and not cur.endswith("\\"+q): q=None
            continue
        if ch in "\"'": q=ch;cur+=ch;continue
        if ch in "([{": d+=1
        if ch in ")]}": d-=1
        if ch=="," and d==0: out.append(cur);cur="";continue
        cur+=ch
    if cur.strip(): out.append(cur)
    return [a.strip() for a in out]

def calls(text):
    res=[]
    for m in re.finditer(r'\bexpect_[A-Za-z_]+\s*\(', text):
        name=text[m.start():m.end()-1].strip()
        i=m.end();d=1;q=None
        while i<len(text) and d>0:
            ch=text[i]
            if q:
                if ch=="\\": i+=2; continue
                if ch==q: q=None
            elif ch in "\"'": q=ch
            elif ch in "([{": d+=1
            elif ch in ")]}": d-=1
            i+=1
        line=text.count("\n",0,m.start())+1
        res.append((name,line,text[m.end():i-1]))
    return res

LIT=re.compile(r'^\s*(-?[0-9][0-9eE\.\-+_]*L?|TRUE|FALSE|NA|NA_real_|NA_integer_|NA_character_|NULL|Inf|-Inf|NaN|"[^"]*"|c\(\s*(-?[0-9][^()]*)\)|list\(\s*\)|character\(0\)|integer\(0\)|numeric\(0\))\s*$')
SHAPE=re.compile(r'\b(dim|length|nrow|ncol|names|dimnames|class|is\.numeric|is\.character|is\.logical|is\.matrix|is\.array|is\.list|is\.function|is\.factor|inherits|typeof|levels|colnames|rownames|nlevels)\s*\(')
FIN=re.compile(r'\b(is\.finite|is\.nan|is\.infinite|is\.na)\s*\(')

def classify(name,args,src):
    a=split_args(args)
    a0=a[0] if a else ""
    a1=a[1] if len(a)>1 else ""
    if name in ("expect_error","expect_warning","expect_message","expect_condition"):
        return "refusal_pat" if len(a)>1 and ('"' in a1 or "pattern" in args) else "refusal_bare"
    if name=="expect_silent": return "smoke_silent"
    if name in ("expect_null",): return "shape"
    if name=="expect_true" or name=="expect_false":
        s=a0
        if re.match(r'^\s*all\s*\(',s):
            inner=s.strip()[4:-1] if s.strip().endswith(")") else s
            if FIN.search(s): return "all_finite"
            return "all_true"
        if FIN.search(s): return "finite"
        if SHAPE.search(s): return "shape"
        if re.search(r'[<>]=?|abs\s*\(|max\s*\(|mean\s*\(|sqrt\s*\(',s): return "bound"
        return "bool_other"
    if name in ("expect_equal","expect_identical","expect_equivalent"):
        lit0=bool(LIT.match(a0)); lit1=bool(LIT.match(a1))
        if lit0 or lit1:
            if SHAPE.search(a0) or SHAPE.search(a1): return "shape_pin"
            return "content_pin"
        if SHAPE.search(a0) and SHAPE.search(a1): return "shape_pin"
        # both sides expressions
        return "expr_vs_expr"
    return "other:"+name

files=[l.strip() for l in open(sys.argv[1]) if l.strip()]
order=["content_pin","shape_pin","shape","expr_vs_expr","bound","all_true","all_finite","finite","bool_other","refusal_pat","refusal_bare","smoke_silent"]
tot={}
for f in files:
    p=os.path.join(root,f)
    t=open(p).read()
    cs=calls(t)
    cnt={}
    for name,line,args in cs:
        k=classify(name,args,t)
        cnt[k]=cnt.get(k,0)+1
        tot[k]=tot.get(k,0)+1
    n=len(cs)
    print("%-52s %4d  %s"%(os.path.basename(f),n," ".join("%s=%d"%(k,cnt[k]) for k in order if k in cnt)+(" | "+" ".join("%s=%d"%(k,v) for k,v in cnt.items() if k not in order) if any(k not in order for k in cnt) else "")))
print()
print("TOTAL", sum(tot.values()), sorted(tot.items(), key=lambda x:-x[1]))
