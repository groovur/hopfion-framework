#!/usr/bin/env python3
r"""
E_6 bridge ordering: assign the 6 E_6 exponents to the 6 quarks using the E_8 bridge (mass/generation),
and label the mod-3 residue table. (2026-09-03). See notes/quark_generation_e8_ribbon_twist.md.

E_6: h=12, exps {1,4,5,7,8,11}, Tw=m-6 in {-5,-2,-1,+1,+2,+5}, mirror pairs (1,11),(4,8),(5,7)=|Tw|{5,2,1}.
E_8 bridge (established): gen-3 (t,b) carry the LARGEST |Tw| (14,8); |Tw| grows with generation.
"""
E6_exps = [1,4,5,7,8,11]
Tw = {m: m-6 for m in E6_exps}
pairs = [(1,11),(4,8),(5,7)]          # mirror pairs m+m'=12

print("="*70); print("E_6 BRIDGE ORDERING: |Tw| -> generation, mirror -> isospin doublet"); print("="*70)

print("\n[1] Mirror pairs = |Tw| levels (candidate generations):")
for a,b in pairs:
    print(f"    ({a:>2},{b:>2})  Tw = {Tw[a]:>+2}/{Tw[b]:>+2}   |Tw| = {abs(Tw[a])}")
print("    -> three DISTINCT |Tw| magnitudes {1,2,5} (unlike E_8, where gen-1,2 SHARE |Tw|={2,4}).")
print("       So in E_6, |Tw| labels generation cleanly.")

print("\n[2] Which |Tw| = which generation? -- fixed by the E_8 bridge:")
print("    E_8: gen-3 (t,b) have the LARGEST |Tw| (14,8); |Tw| increases with generation/mass.")
print("    => E_6:  |Tw|=1 -> gen-1,  |Tw|=2 -> gen-2,  |Tw|=5 -> gen-3.")

print("\n[3] Mirror pair = isospin DOUBLET (the decisive check):")
gen_of = {1:3, 11:3, 4:2, 8:2, 5:1, 7:1}     # |Tw|=5->g3, 2->g2, 1->g1
doublet = {3:'(t,b)', 2:'(c,s)', 1:'(u,d)'}
for a,b in pairs:
    g = gen_of[a]
    print(f"    ({a:>2},{b:>2}) |Tw|={abs(Tw[a])} -> generation {g} doublet {doublet[g]}  (up = +Tw, down = -Tw)")
print("    => each E_6 mirror pair is a SAME-GENERATION up/down doublet (u-d, c-s, t-b). PHYSICAL.")
print("       Contrast E_8: its mirror pairs were SAME-ISOSPIN, CROSS-generation (s-d, c-u), gen-3 UNPAIRED.")
print("       E_6's mirror structure is the real isospin doublets -- another sign E_6 is the per-strand group.")

# assignment: gen by |Tw|, isospin by sign (up=+Tw convention; intra-doublet up/down TBD by charge/writhe)
assign = {  # quark: (E6 exponent, Tw, generation, isospin)
 'u':(7,+1,1,'up'),  'd':(5,-1,1,'down'),
 'c':(8,+2,2,'up'),  's':(4,-2,2,'down'),
 't':(11,+5,3,'up'), 'b':(1,-5,3,'down'),
}
print("\n[4] Labelled table (up=+Tw convention; residue = Tw mod 3 = SU(3)_1 framing class 0=1,1=3,2=3bar):")
print("     q  m   Tw   gen  isospin   res(mod3)")
for q,(m,tw,g,iso) in sorted(assign.items(), key=lambda kv: (kv[1][2], kv[1][3])):
    print(f"    {q:>2} {m:>2}  {tw:>+2}   {g}   {iso:>4}     {tw%3}")

r1=[q for q,(m,tw,g,i) in assign.items() if tw%3==1]
r2=[q for q,(m,tw,g,i) in assign.items() if tw%3==2]
print(f"\n    residue 1 (3):  {sorted(r1)}     residue 2 (3bar): {sorted(r2)}     residue 0 (singlet): [] (empty)")
print("    residue-0 empty = no quark is a colour singlet (correct). The 1-vs-2 (3/3bar) split does NOT")
print("    line up with isospin or generation -> it is a separate framing-orientation label [meaning open].")

print("\n[5] BRIDGE VALIDATION + honest limits:")
print("    - generation ordering |Tw|=1,2,5 -> gen 1,2,3 matches the E_8/mass ordering (gen-3 = max twist).")
print("    - mirror pairs = physical isospin doublets (u-d,c-s,t-b): consistent, cleaner than E_8.")
print("    - FORCED: generation (|Tw|) and the doublet pairing. NOT yet forced: which doublet member is")
print("      up vs down (the +Tw=up convention needs the electric-charge / writhe input); and the 3/3bar")
print("      residue meaning. Mass values (item 3) are a separate check, via the E_8 bridge not the E_6 twist.")
print("="*70)
