"""Extract the order of first citation of every Figure / Figure S / File S / Table token."""
import re, sys

pat = re.compile(r'(Figure S\d+[A-Z]?(?:[-,–] ?\d?[A-Z])?|Figure \d+[A-Z]?|File S\d+|Supplementary File \d+|Supplementary Figure \d+|Table \d+|Table S\d+)')

def scan(path):
    text = open(path).read()
    hits = pat.findall(text)
    return hits

for path in sys.argv[1:]:
    hits = scan(path)
    print("=" * 70)
    print(path, len(hits), "tokens")
    # normalise to the container (drop panel letter)
    norm = []
    for h in hits:
        m = re.match(r'((?:Figure S|Figure |File S|Table S|Table |Supplementary Figure |Supplementary File )\d+)', h)
        norm.append(m.group(1) if m else h)
    seen = []
    for n in norm:
        if n not in seen:
            seen.append(n)
    print("first-appearance order:")
    print("  " + " | ".join(seen))
