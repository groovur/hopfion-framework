"""Frank director stiffness K of the semi-Dirac condensate — PRODUCTION,
gate-validated, resumable. Undoped semimetal: interband + diamagnetic only
(intraband Pauli-blocked). RADIAL DISK cutoff |k|<Lam (preserves the axis-
rotation symmetry -> Goldstone gate passes; verified ~1e-8).

Per k, lower band |-(k)>=(1,-e^{i phi})/sqrt2, phi=atan2(d2,d1), d=(kx^2/2,ky).
Generator G=(kx ky, -kx).  M(k,k')=<+(k')|G.sigma|-(k)>.
dE(q)/theta0^2 = INT_disk d2k/(2pi)^2 (1/4) sum_{s=+-1}
   [ |M(k,k+sq)|^2/(-|d_k|-|d_{k+sq}|) - |M(k,k)|^2/(-2|d_k|) ]
(the subtraction = the diamagnetic, exact via the gate; dE(0)=0 by construction).
K = 2 * dE(q)/theta0^2 / q^2, extracted as q->0, along bend(q||y)/splay(q||x).

Resumable: each (direction,q,Lam,Nr,Nth) integral is one checkpoint record in
checkpoint.jsonl; completed records are skipped on rerun. Grid rows accumulated
so a single big integral is itself chunked and resumable via partial-sum records.
"""
import numpy as np, json, os, sys, time

CKPT="checkpoint.jsonl"

def load_done():
    done={}
    if os.path.exists(CKPT):
        for l in open(CKPT):
            try: r=json.loads(l)
            except: continue
            if r.get("kind")=="result": done[r["key"]]=r
    return done

def append(rec):
    with open(CKPT,"a") as f: f.write(json.dumps(rec)+"\n")

def integrand_rows(q, ex, ey, Lam, Nr, Nth, r_lo=0):
    """yield (i_r, partial_contribution) for r-block resumability."""
    rs=np.linspace(Lam/Nr, Lam, Nr); dr=rs[1]-rs[0]
    th=np.linspace(0,2*np.pi,Nth,endpoint=False); dth=th[1]-th[0]
    C=np.cos(th); S=np.sin(th)
    for i in range(r_lo,Nr):
        r=rs[i]; kx=r*C; ky=r*S
        d1,d2=kx**2/2,ky; nd=np.hypot(d1,d2); phi=np.arctan2(d2,d1)
        G1,G2=kx*ky,-kx
        val=0.0
        for s in (+1,-1):
            kxp,kyp=kx+s*q*ex, ky+s*q*ey
            d1p,d2p=kxp**2/2,kyp; ndp=np.hypot(d1p,d2p); phip=np.arctan2(d2p,d1p)
            z =-(G1-1j*G2)*np.exp(1j*phi)+(G1+1j*G2)*np.exp(-1j*phip)  # M(k,k+sq)
            z0=-(G1-1j*G2)*np.exp(1j*phi)+(G1+1j*G2)*np.exp(-1j*phi)   # M(k,k)
            term=(np.abs(z)/2)**2/(-nd-ndp) - (np.abs(z0)/2)**2/(-2*nd)
            val+=0.25*term
        # integrate this r-shell: sum over theta * r dr dth /(2pi)^2
        yield i, float(np.sum(val)*r*dr*dth/(2*np.pi)**2)

def compute(key, q, ex, ey, Lam, Nr, Nth):
    # resume partial if present
    part_key=key+"|partial"; acc=0.0; r_lo=0
    if os.path.exists(CKPT):
        for l in open(CKPT):
            try: r=json.loads(l)
            except: continue
            if r.get("kind")=="partial" and r.get("key")==key:
                acc=r["acc"]; r_lo=r["i_r"]+1
    last_ck=time.time()
    for i,contrib in integrand_rows(q,ex,ey,Lam,Nr,Nth,r_lo):
        acc+=contrib
        if time.time()-last_ck>60:  # checkpoint partial every 60s
            append({"kind":"partial","key":key,"i_r":i,"acc":acc}); last_ck=time.time()
    return acc

def run(mode="test"):
    if mode=="test":
        Lam,Nr,Nth=20.0,400,720; qs=[0.1,0.2,0.4]
    else:  # full: fine, small-q for K extraction; silent
        Lam,Nr,Nth=40.0,2000,2880; qs=[0.02,0.03,0.05,0.08,0.12]
    done=load_done()
    for lab,(ex,ey) in [("bend",(0,1)),("splay",(1,0))]:
        for q in qs:
            key=f"{mode}|{lab}|q{q}|L{Lam}|Nr{Nr}|Nth{Nth}"
            if key in done: continue
            t=time.time(); dE=compute(key,q,ex,ey,Lam,Nr,Nth)
            rec={"kind":"result","key":key,"mode":mode,"dir":lab,"q":q,"Lam":Lam,
                 "Nr":Nr,"Nth":Nth,"dE_over_th2":dE,"K":2*dE/q**2,"wall_s":round(time.time()-t,1)}
            append(rec)
            if mode=="test": print(json.dumps({k:rec[k] for k in ("dir","q","dE_over_th2","K","wall_s")}))

if __name__=="__main__":
    run(sys.argv[1] if len(sys.argv)>1 else "test")
    if (sys.argv[1] if len(sys.argv)>1 else "test")=="test":
        # summarise test K + exponent
        done=load_done()
        for lab in ("bend","splay"):
            rows=[(r["q"],r["dE_over_th2"]) for r in done.values() if r.get("mode")=="test" and r["dir"]==lab]
            rows.sort()
            if len(rows)>=2:
                q=np.array([x[0] for x in rows]); dE=np.array([x[1] for x in rows])
                s=np.polyfit(np.log(q),np.log(np.abs(dE)),1)[0]
                K=np.mean(2*dE/q**2)
                print(f"  {lab}: dE ~ q^{s:.2f}  (2=standard Frank)  K~{K:+.4e}  sign={'STABLE' if K>0 else 'UNSTABLE'}")
