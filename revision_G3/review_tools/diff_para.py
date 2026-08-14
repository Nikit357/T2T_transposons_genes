"""Sentence-level diff between two paragraph-per-line text extracts."""
import sys, re, difflib

def sentences(path):
    out = []
    for i, line in enumerate(open(path)):
        line = line.strip()
        if not line:
            continue
        # split on sentence enders followed by space+capital / digit
        parts = re.split(r'(?<=[.:;])\s+(?=[A-Z0-9(])', line)
        for p in parts:
            p = p.strip()
            if p:
                out.append(p)
    return out

a = sentences(sys.argv[1])
b = sentences(sys.argv[2])
sm = difflib.SequenceMatcher(None, a, b, autojunk=False)
for tag, i1, i2, j1, j2 in sm.get_opcodes():
    if tag == 'equal':
        continue
    print(f"\n@@@ {tag.upper()}  claude[{i1}:{i2}]  daniil[{j1}:{j2}]")
    for s in a[i1:i2]:
        print("  - " + s)
    for s in b[j1:j2]:
        print("  + " + s)
