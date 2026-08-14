"""Word-level diff with context, for two paragraph-per-line extracts."""
import sys, re, difflib

def words(path):
    t = open(path).read()
    return re.findall(r"\S+|\n", t)

a = words(sys.argv[1])
b = words(sys.argv[2])
sm = difflib.SequenceMatcher(None, a, b, autojunk=False)
CTX = 8
for tag, i1, i2, j1, j2 in sm.get_opcodes():
    if tag == 'equal':
        continue
    if (i2 - i1) > 400 or (j2 - j1) > 400:
        print(f"\n@@@ {tag.upper()} [large block: {i2-i1} -> {j2-j1} words, skipped]")
        continue
    pre = " ".join(w for w in a[max(0, i1 - CTX):i1] if w != "\n")
    post = " ".join(w for w in a[i2:i2 + CTX] if w != "\n")
    old = " ".join(w for w in a[i1:i2] if w != "\n")
    new = " ".join(w for w in b[j1:j2] if w != "\n")
    print(f"\n@@@ {tag.upper()}")
    print(f"  ctx< ...{pre}")
    print(f"  -    {old}")
    print(f"  +    {new}")
    print(f"  >ctx {post}...")
