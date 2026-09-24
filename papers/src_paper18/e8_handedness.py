#!/usr/bin/env python3
r"""
Handedness of the twist assignment: is the gen-3 sign asymmetry a fundamental chirality, or the
quark mass hierarchy in disguise? And what integer framing does it imply? (2026-09-03).
See notes/quark_generation_e8_ribbon_twist.md.

Twist Tw(m)=m-15 (integer, from e8_twist_framing.py). Used exponents {1,7,11,13,17,19}, unused {23,29}.
"""
exps_all = [1,7,11,13,17,19,23,29]
assign = {'u':(1,19),'c':(2,11),'t':(3,1),'d':(1,17),'s':(2,13),'b':(3,7)}
pdg    = {'u':2.16,'d':4.67,'s':93.4,'c':1275,'b':4180,'t':172570}   # MeV
Tw = lambda m: m-15

print("="*72); print("HANDEDNESS CHECK: twist sign, parity, and the mass hierarchy"); print("="*72)

print("\n[1] PARITY of the twist (depth for the integer-framing step):")
print("    E_8 exponents are all ODD -> Tw=m-15 is always EVEN.")
for m in exps_all:
    print(f"      m={m:>2}  Tw={Tw(m):>+3}  {'even' if Tw(m)%2==0 else 'ODD'}")
print("    => every framing is an EVEN integer. A ribbon with even framing closes with a consistent")
print("       orientation/spin structure -- so the twist coordinate is a genuine (even) framing, not a")
print("       half-integer artefact. This is the constraint the integer-framing step needs.")

print("\n[2] TWIST vs MASS (is the sign chirality, or heavy/light?):")
print("     q  gen  Tw    PDG mass(MeV)")
for q,(g,m) in sorted(assign.items(), key=lambda kv: Tw(kv[1][1])):
    print(f"    {q:>2}   {g}  {Tw(m):>+3}   {pdg[q]:>9.1f}")
print("    => monotone: MORE-NEGATIVE twist = HEAVIER quark (t at -14 heaviest, u at +4 lightest).")
print("       The twist SIGN is just light(+)/heavy(-), i.e. the mass ordering -- NOT a separate chirality.")

print("\n[3] USED vs UNUSED as the mass floor:")
print("    used twists   {-14,-8,-4,-2,+2,+4} = the 6 quarks (heavy end down to u,d).")
print("    unused twists {+8,+14} (exps 23,29) = the POSITIVE (light) extremes -> quarks LIGHTER than u,d.")
print("    => the 'gen-3 handedness' is the mass FLOOR: there are very heavy quarks (t,b at -8,-14) but")
print("       no ultra-light mirror partners, so the positive extremes are empty. It is the flavour")
print("       hierarchy asymmetry in twist language, NOT a fundamental left/right asymmetry. [deflates]")

print("\n[4] The REAL structural feature: gen-1 <-> gen-2 are exact twist-mirrors, gen-3 is unpaired:")
mirror = {}
for q,(g,m) in assign.items():
    for q2,(g2,m2) in assign.items():
        if Tw(m2)==-Tw(m) and q2!=q: mirror[q]=q2
print(f"    twist-mirror partners (Tw <-> -Tw): {mirror}")
print("    (u,c)=(+4,-4) and (d,s)=(+2,-2) are exact mirror pairs = gen-1 <-> gen-2, SAME isospin.")
print("    gen-3 (t,b)=(-14,-8) have NO partner (mirrors 23,29 empty).")
print("    PARALLEL [speculative]: the (1,2) generations mix most strongly (Cabibbo ~13deg) while the")
print("    3rd generation mixes weakly (small CKM to t,b). Twist-mirror-paired <-> strongly-mixed;")
print("    twist-unpaired <-> weakly-mixed. Suggestive that the twist pairing seeds the mixing hierarchy.")

print("\n" + "="*72)
print("VERDICT: the gen-3 'handedness' is the mass hierarchy restated (heavy quarks = negative twist,")
print("no ultra-light mirrors) -- NOT a new chirality. BUT two real things fall out:")
print("  (a) all framings are EVEN integers  -> the integer-framing step should use even Calugareanu")
print("      framings (consistent spin structure);")
print("  (b) gen-1<->gen-2 twist-mirror pairing with gen-3 unpaired -> a possible seed of the CKM")
print("      mixing hierarchy (paired=strongly mixed, unpaired=weakly mixed). [lead, not established]")
print("="*72)
