#!/usr/bin/env python3
"""预测 (i) 决定性检验（修正解析）。

XML 真实结构：<mg name="GE_ART"><m name="Wave">BASE64</m><m name="Points">180</m>...</mg>
（无 </waveform> 标签——首版正则因此全空。这是第 3 次踩同类坑：先看清结构再写解析。）
"""
import sys, os, csv, json, base64, re, numpy as np
from concurrent.futures import ProcessPoolExecutor

RAW=os.environ.get("MOVER_WAVE_ROOT", "/mnt/share/ljfang/exp/raw/wave_pilot_20260909")
QC =os.environ.get("MOVER_ECG_QC", "/mnt/share/ljfang/exp/results/run-20260912-audit/20_wp2_ecg/ecg_file_qc_full.csv")
MG = re.compile(r'<mg name="(GE_ART|GE_ECG)">(.*?)</mg>', re.S)
WAVE = re.compile(r'<m name="Wave">([^<]+)</m>')

def lag1(x):
    x=np.asarray(x,float)
    if x.size<50: return None
    x=x-x.mean(); den=np.dot(x,x)
    if den<=0: return None
    return float(np.dot(x[:-1],x[1:])/den)

def classify(path, want="GE_ART", npmax=8):
    try:
        with open(path,encoding="utf-8",errors="replace") as f: s=f.read()
    except Exception: return None
    acs=[]
    for name,body in MG.findall(s):
        if name!=want: continue
        m=WAVE.search(body)
        if not m: continue
        try:
            raw=base64.b64decode(m.group(1))
            v=np.frombuffer(raw,dtype="<i2").astype(float)
            d=lag1(v)
            if d is not None: acs.append(d)
        except Exception: pass
        if len(acs)>=npmax: break
    return float(np.median(acs)) if acs else None

def work(a):
    rel,alias=a
    p=os.path.join(RAW,rel)
    return (rel,alias,classify(p,"GE_ART"),classify(p,"GE_ECG"))

if __name__=="__main__":
    rows=[]
    with open(QC) as f:
        for r in csv.DictReader(f):
            try: a=float(r["alias_pct"])
            except Exception: continue
            rows.append((r["relpath"],a))
    rng=np.random.default_rng(20260911)
    sel=[rows[i] for i in rng.choice(len(rows),size=200,replace=False)]
    print(f"抽样 {len(sel)} 例")
    out=[]
    with ProcessPoolExecutor(max_workers=16) as ex:
        for k,r in enumerate(ex.map(work, sel, chunksize=4)):
            if r[2] is not None and r[3] is not None: out.append(r)
            if (k+1)%50==0: print(f"  ...{k+1} 有效 {len(out)}",flush=True)
    print(f"两通道均有效 {len(out)}")
    if out:
        aac=np.array([o[2] for o in out]); eac=np.array([o[3] for o in out]); al=np.array([o[1] for o in out])
        artA=aac<-0.5; ecgA=al>50
        both=int((artA&ecgA).sum()); neither=int(((~artA)&(~ecgA)).sum())
        oA=int((artA&~ecgA).sum()); oE=int((~artA&ecgA).sum()); n=len(out)
        po=(both+neither)/n
        pe=(artA.mean()*ecgA.mean())+((1-artA.mean())*(1-ecgA.mean()))
        kappa=(po-pe)/(1-pe) if pe<1 else float('nan')
        print(f"\n  GE_ART 受损 {int(artA.sum())}/{n} ({100*artA.mean():.1f}%)   ac1 中位 {np.median(aac):+.3f}")
        print(f"  GE_ECG 受损 {int(ecgA.sum())}/{n} ({100*ecgA.mean():.1f}%)   ac1 中位 {np.median(eac):+.3f}")
        print(f"  都坏 {both} | 都好 {neither} | 只 ART 坏 {oA} | 只 ECG 坏 {oE}")
        print(f"  Po={po:.3f}  Pe={pe:.3f}  **κ={kappa:+.3f}**")
        print(f"\n  解释：κ>0.8 强一致（支持'同一打包路径'）；κ≈0 独立（反驳）")
        json.dump(dict(n=n,art_bad=int(artA.sum()),ecg_bad=int(ecgA.sum()),both=both,
                       neither=neither,only_art=oA,only_ecg=oE,po=float(po),pe=float(pe),
                       kappa=float(kappa),art_ac1=float(np.median(aac)),
                       ecg_ac1=float(np.median(eac))),
                  open("/tmp/chanpred_result.json","w"),indent=2)
        print("wrote /tmp/chanpred_result.json")
