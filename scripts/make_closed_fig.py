#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build the 3-panel closed-form figure from ../data/closed_S10_data.csv (produced
by closed_form.py with default MAXK=10). Left/centre: measured survival vs the
closed form K*prod f_q for 42 and 210. Right: residuals (<=0.5% for omega<=5,
edge ~2.5% at omega=6).
"""
import csv, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
rows=list(csv.DictReader(open('../data/closed_S10_data.csv')))
def series(gap):
    om=[];m=[];c=[];e=[]
    for x in rows:
        if int(x['gap'])==gap:
            om.append(int(x['omega'])); m.append(float(x['measured_P']))
            c.append(float(x['closed_form'])); e.append(float(x['err_pct']))
    return map(np.array,(om,m,c,e))
fig,axes=plt.subplots(1,3,figsize=(17,4.8))
for ax,gap,col in [(axes[0],42,'#185FA5'),(axes[1],210,'#c0392b')]:
    om,m,c,e=series(gap)
    ax.plot(om,m,'o-',color=col,lw=2.2,ms=8,label='measured $P(N{+}d\\,\\mathrm{twin}\\mid\\omega)$')
    ax.plot(om,c,'s--',color='#e8845b',lw=1.7,ms=7,label='closed form  $K\\cdot\\prod_q f_q$')
    me=np.max(np.abs(e))
    ax.set_title(f'$6\\Delta N={gap}$  (max err {me:.1f}%)',fontsize=11)
    ax.set_xlabel(r'$\omega_{>3}(N)$',fontsize=11); ax.set_ylabel(r'absolute right-centre survival $P$',fontsize=10.5)
    ax.legend(fontsize=9,loc='best'); ax.grid(alpha=.25); ax.set_xticks(range(1,7))
om1,m1,c1,e1=series(42); om2,m2,c2,e2=series(210)
axes[2].axhline(0,color='gray',lw=1)
axes[2].axhspan(-1,1,color='#2ca25f',alpha=.10,label='$\\pm1\\%$ band')
axes[2].plot(om1,e1,'o-',color='#185FA5',lw=2,ms=7,label='42')
axes[2].plot(om2,e2,'D-',color='#c0392b',lw=2,ms=7,label='210')
axes[2].set_title('closed-form residual: $\\leq0.5\\%$ for $\\omega\\leq5$, edge at $\\omega{=}6$',fontsize=11)
axes[2].set_xlabel(r'$\omega_{>3}(N)$',fontsize=11); axes[2].set_ylabel('error  (%)',fontsize=11)
axes[2].legend(fontsize=9); axes[2].grid(alpha=.25); axes[2].set_xticks(range(1,7)); axes[2].set_ylim(-4,4)
plt.suptitle('Closed-form absolute survival in $S_{10}$:  $P(N{+}d\\,\\mathrm{twin}\\mid\\omega)=K\\cdot\\prod_q f_q(d,N)$,  pure CRT factors, no fitted parameters',fontsize=12,y=1.02)
plt.tight_layout()
plt.savefig('fig_paper7_closed.pdf',bbox_inches='tight')
plt.savefig('fig_paper7_closed.png',dpi=160,bbox_inches='tight')
print("figure saved")
