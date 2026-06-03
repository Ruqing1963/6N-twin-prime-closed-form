#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Closed-form absolute right-centre survival  (resolves the Part VI open problem).

  P(N+d twin | N twin, omega) = K(d) * prod_q f_q(d, N)

where the product runs over small primes q>3, f_q is the PURE CRT conditional
safety of the right centre M=N+d given the left centre's divisibility by q:
    q | N : f_q = 1 if (d mod q) is q-safe (d not in dead(q)) else 0
    q ! N : f_q = (# admissible nonzero residues r with (r+d) safe) / (# admissible nonzero r)
        with dead(q) = {+-6^{-1} mod q},
and K(d) is an omega-INDEPENDENT tail constant (primes beyond the POOL). No
fitted parameters: f_q is closed-form CRT; K is fixed once from the omega-merged
ratio. Combined with the bridge constant C(d) of Part VI, this gives a fully
closed-form r(d|omega) = (K/C) * <prod_q f_q>_omega.

NOTE: f_qnot uses a uniform-over-admissible-residues approximation; any residual
(a few % in S9) is expected to be the Part-I enrichment correction to that
uniformity. Check whether S10 tightens the error and whether K(d) is the same for
42 and 210 (i.e. d-independent).

USAGE: python abs_surv_closed_S10.py   (default S10); MAXK=9 for S9.
Requires: numpy.
"""
# Closed-form absolute survival: P(N+d twin|omega) = K * prod_q f_q(d,N),
# where f_q is the PURE CRT conditional safety:
#   q|N: f = 1 if d mod q not in dead(q) else 0
#   q∤N: f = (#admissible nonzero r with (r+d) safe)/(#admissible nonzero r)
# K is an omega-independent tail constant (primes > max(POOL)). Test vs measured.
import numpy as np, math, os
def primes_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(math.isqrt(n))+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0].astype(np.int64)
MAXK=int(os.environ.get("MAXK",10))
LO=10**(MAXK-1)//6+1; HI=10**MAXK//6; SEG=4_000_000
PB=int(math.isqrt(6*HI+250))+1; BP=primes_upto(PB)
POOL=[5,7,11,13,17,19,23,29,31,37,41,43,47]
poolbit={q:1<<i for i,q in enumerate(POOL)}
def dead(q):
    inv=pow(6,-1,q); return {inv%q,(-inv)%q}
def f_qN(d,q): return 1.0 if (d%q) not in dead(q) else 0.0
def f_qnot(d,q):
    dq=dead(q); adm=[r for r in range(q) if r not in dq and r!=0]
    return sum(1 for r in adm if (r+d)%q not in dq)/len(adm) if adm else 0.0
N_list=[]; om_list=[]; mask_list=[]
n=LO
while n<=HI:
    nh=min(n+SEG,HI+1); sz=nh-n
    rem=np.arange(n,nh,dtype=np.int64); ob=np.zeros(sz,np.int16); mk=np.zeros(sz,np.int32)
    for p in BP:
        if p*p>nh-1: break
        f=((n+p-1)//p)*p
        if f>=nh: continue
        idx=np.arange(f-n,sz,p)
        if idx.size==0: continue
        sub=rem[idx]; m=(sub%p)==0
        while m.any(): sub[m]//=p; m=(sub%p)==0
        rem[idx]=sub
        if p>3:
            ob[idx]+=1
            if p in poolbit: mk[idx]|=poolbit[p]
    ob[rem>1]+=1
    Narr=np.arange(n,nh,dtype=np.int64)
    vlo=6*n-1; vhi=6*(nh-1)+1; span=vhi-vlo+1
    comp=np.zeros(span,bool); sq=int(math.isqrt(vhi))+1
    for p in BP:
        if p>sq: break
        st=max(p*p,((vlo+p-1)//p)*p)
        if st>vhi: continue
        comp[st-vlo:span:p]=True
    tw=(~comp[(6*Narr-1)-vlo])&(~comp[(6*Narr+1)-vlo])
    pos=np.nonzero(tw)[0]
    N_list.append(Narr[pos]); om_list.append(ob[pos]); mask_list.append(mk[pos])
    n=nh
N_arr=np.concatenate(N_list); om_arr=np.concatenate(om_list).astype(np.int16); mask_arr=np.concatenate(mask_list)
print(f"S{MAXK} twins {len(N_arr):,}")
def is_twin(vals):
    idx=np.searchsorted(N_arr,vals); idx=np.clip(idx,0,len(N_arr)-1)
    return N_arr[idx]==vals
for d in [7,35]:
    sixd=6*d
    FqN={q:f_qN(d,q) for q in POOL}; Fqnot={q:f_qnot(d,q) for q in POOL}
    # determine K from omega-merged: K = mean over all twins of [measured M-twin] / [pure product]
    prod_all=np.ones(len(N_arr))
    for q in POOL:
        qb=poolbit[q]; has=(mask_arr&qb)>0
        prod_all*=np.where(has,FqN[q],Fqnot[q])
    measM_all=is_twin(N_arr+d)
    K=measM_all.mean()/prod_all.mean()
    print(f"\n=== 6dN={sixd}: closed-form P = K * prod(CRT f_q),  K={K:.5f} ===")
    print(f"  {'omega':>5}{'measured P':>12}{'closed-form':>13}{'err%':>8}")
    maxe=0
    for om in range(1,7):
        band=(om_arr==om)
        if band.sum()<20000: continue
        Nb=N_arr[band]; mkb=mask_arr[band]
        measP=is_twin(Nb+d).mean()
        prod=np.ones(len(Nb))
        for q in POOL:
            qb=poolbit[q]; has=(mkb&qb)>0
            prod*=np.where(has,FqN[q],Fqnot[q])
        cf=K*prod.mean()
        e=100*(cf-measP)/measP
        if abs(e)>abs(maxe): maxe=e
        print(f"  {om:>5}{measP:>12.5f}{cf:>13.5f}{e:>8.1f}")
    print(f"  max |err| = {abs(maxe):.1f}%")

# ---- emit CSV (closed_S{K}_data.csv) ----
import csv as _csv
with open(f'closed_S{MAXK}_data.csv','w',newline='') as _f:
    _w=_csv.writer(_f); _w.writerow(['gap','omega','measured_P','closed_form','err_pct'])
    for d in [7,35]:
        FqN={q:f_qN(d,q) for q in POOL}; Fqnot={q:f_qnot(d,q) for q in POOL}
        prod_all=np.ones(len(N_arr))
        for q in POOL:
            qb=poolbit[q]; has=(mask_arr&qb)>0
            prod_all*=np.where(has,FqN[q],Fqnot[q])
        K=is_twin(N_arr+d).mean()/prod_all.mean()
        for om in range(1,8):
            band=(om_arr==om)
            if band.sum()<20000: continue
            Nb=N_arr[band]; mkb=mask_arr[band]
            measP=is_twin(Nb+d).mean()
            prod=np.ones(len(Nb))
            for q in POOL:
                qb=poolbit[q]; has=(mkb&qb)>0
                prod*=np.where(has,FqN[q],Fqnot[q])
            cf=K*prod.mean(); e=100*(cf-measP)/measP if measP>0 else 0
            _w.writerow([6*d,om,f'{measP:.5f}',f'{cf:.5f}',f'{e:.1f}'])
print(f"\n[ok] wrote closed_S{MAXK}_data.csv")
