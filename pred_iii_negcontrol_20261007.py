#!/usr/bin/env python3
"""预测 (iii)：另一条通道集（CB- 文件）应无此痕迹 —— 同语料内的负对照。"""
import sys, os, json, base64, re, glob, collections, numpy as np
RAW=os.environ.get("MOVER_WAVE_ROOT", "/mnt/share/ljfang/exp/raw/wave_pilot_20260909")
MG=re.compile(r'<mg name="([A-Z0-9_]+)">(.*?)</mg>',re.S)
WAVE=re.compile(r'<m name="Wave">([^<]+)</m>')
def lag1(x):
    x=np.asarray(x,float)
    if x.size<50: return None
    x=x-x.mean(); d=np.dot(x,x)
    return float(np.dot(x[:-1],x[1:])/d) if d>0 else None
def cls(path,want,n=4):
    try: s=open(path,encoding="utf-8",errors="replace").read()
    except Exception: return None
    acs=[]
    for name,body in MG.findall(s):
        if name!=want: continue
        m=WAVE.search(body)
        if not m: continue
        try:
            v=np.frombuffer(base64.b64decode(m.group(1)),dtype="<i2").astype(float)
            d=lag1(v)
            if d is not None: acs.append(d)
        except Exception: pass
        if len(acs)>=n: break
    return float(np.median(acs)) if acs else None

cbs=[]
for part in ("epic_wave_1_v2","epic_wave_2_v2","epic_wave_3_v2"):
    cbs+=glob.glob(f"{RAW}/{part}/**/*CB-*.xml",recursive=True)
    if len(cbs)>=200: break
print(f"CB- 文件候选 {len(cbs)}")
rng=np.random.default_rng(7)
pick=[cbs[i] for i in rng.choice(len(cbs),size=min(60,len(cbs)),replace=False)]
byc=collections.defaultdict(list)
for i,p in enumerate(pick):
    for ch in ("ECG1","INVP1","PLETH","GE_ART","GE_ECG"):
        v=cls(p,ch,4)
        if v is not None: byc[ch].append(v)
    if (i+1)%20==0: print(f"  ...{i+1}/{len(pick)}",flush=True)
print("\nCB- 文件内各通道分类（负对照）：")
for ch,vs in sorted(byc.items()):
    a=np.array(vs); print(f"  {ch:>6}: n={len(a):>3}  ac1 中位={np.median(a):+.3f}  受损={int((a<-0.5).sum())}")
json.dump({k:[float(x) for x in v] for k,v in byc.items()}, open("/tmp/chanpred_cb.json","w"),indent=1)
print("wrote /tmp/chanpred_cb.json")
