"""Static (omega=0), undoped polarization Pi(q) for the 2D semi-Dirac two-band
Hamiltonian H(k) = (k_x^2/2m) sigma_x + (v k_y) sigma_y  [hbar=1, m=v=1].
E(k)=|d|, d=(k_x^2/2, k_y). Interband Lindhard at the node:
  Pi(q) = int d^2k/(2pi)^2 * (1 - dhat_k . dhat_{k+q}) / (E_k + E_{k+q})
Pi(0)=0 (vanishing DOS). We extract the small-q behaviour along the two axes:
prediction from anisotropic scaling (k_x~L^{1/2}, k_y~L, E~L, d^2k~L^{3/2}, Pi~L^{1/2}):
  along quadratic axis (q=q_x): Pi ~ q_x^1 (linear, like a Dirac direction)
  along linear axis   (q=q_y): Pi ~ q_y^{1/2} (square-root -> stronger IR)
The square-root is the key: more singular than Dirac's linear-q screening.
"""
import numpy as np

def dvec(kx,ky):
    return np.stack([kx**2/2.0, ky], axis=-1)

def Epair(kx,ky):
    d=dvec(kx,ky); return np.sqrt((d**2).sum(-1)), d

def Pi(qx,qy, Kmax=40.0, N=1400):
    # integrate over k-grid; use non-uniform (denser near 0) to resolve the node
    u=np.linspace(-1,1,N)
    kx=np.sign(u)*(Kmax*u**2); ky=np.sign(u)*(Kmax*u**2)
    KX,KY=np.meshgrid(kx,ky,indexing='ij')
    dkx=np.gradient(kx); dky=np.gradient(ky)
    Wx,Wy=np.meshgrid(dkx,dky,indexing='ij'); W=Wx*Wy
    E1,d1=Epair(KX,KY); E2,d2=Epair(KX+qx,KY+qy)
    n1=d1/(E1[...,None]+1e-30); n2=d2/(E2[...,None]+1e-30)
    cosdt=(n1*n2).sum(-1)
    integ=(1-cosdt)/(E1+E2+1e-30)
    return float(np.sum(integ*W)/(2*np.pi)**2)

if __name__=="__main__":
    print("q            Pi(q_x-axis)     Pi(q_y-axis)")
    qs=[0.02,0.04,0.08,0.16,0.32]
    px=[]; py=[]
    for q in qs:
        a=Pi(q,0.0); b=Pi(0.0,q); px.append(a); py.append(b)
        print(f"{q:<8}{a:>14.5e}{b:>17.5e}")
    lq=np.log(qs)
    sx=np.polyfit(lq,np.log(np.abs(px)),1)[0]
    sy=np.polyfit(lq,np.log(np.abs(py)),1)[0]
    print(f"\nfitted exponents:  quadratic-axis Pi~q^{sx:.2f} (expect ~1)   linear-axis Pi~q^{sy:.2f} (expect ~0.5)")
    print("interpretation: linear-axis sqrt(q) => static dielectric eps(q)=1+V0(q)Pi(q)")
    print("  with V0~1/q^2 (3D grav analog): eps ~ 1 + q^{-1.5} -> strong IR modification of V(r)")
