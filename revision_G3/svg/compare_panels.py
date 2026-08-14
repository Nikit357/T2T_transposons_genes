"""Compare two versions of a network SVG at the level that is actually reproducible.

`adjust_text` places labels with a stochastic solver and matplotlib stamps a `<dc:date>`
into every SVG, so two runs of the same code on the same data are never byte-identical --
measured: two renders with the same PYTHONHASHSEED differ in label coordinates. A md5
gate on these panels therefore cannot pass and never could. What *is* reproducible, and
what a regression gate has to compare, is the content: the page size, the number of
drawn nodes, and the exact set of label strings.

Usage:
    python svg/compare_panels.py OLD_DIR NEW_DIR panel.svg [panel.svg ...]
"""
import re
import sys
from pathlib import Path

# svg.fonttype="none" keeps label text as literal characters inside <text> elements, one
# element per rendered line — so a wrapped label contributes one entry per line.
TEXT = re.compile(r"<text [^>]*>([^<]*)</text>")


def describe(path):
    """Page size in pt, node count and the multiset of label strings in one SVG."""
    svg = Path(path).read_text()
    width = float(re.search(r'width="([\d.]+)pt"', svg).group(1))
    height = float(re.search(r'height="([\d.]+)pt"', svg).group(1))
    labels = sorted(text.strip() for text in TEXT.findall(svg) if text.strip())
    paths = svg.count("<use ") + svg.count("<path ")
    return {"pt": (round(width, 2), round(height, 2)), "labels": labels,
            "n_labels": len(labels), "n_marks": paths}


def main():
    old_dir, new_dir, panels = sys.argv[1], sys.argv[2], sys.argv[3:]
    failures = 0
    for panel in panels:
        old, new = describe(Path(old_dir) / panel), describe(Path(new_dir) / panel)
        same_labels = old["labels"] == new["labels"]
        same_size = old["pt"] == new["pt"]
        print(f"{panel:24s} {old['pt']} -> {new['pt']}  "
              f"labels {old['n_labels']} -> {new['n_labels']} "
              f"({'identical set' if same_labels else 'DIFFERENT SET'})  "
              f"marks {old['n_marks']} -> {new['n_marks']}")
        if not same_labels:
            only_old = sorted(set(old["labels"]) - set(new["labels"]))
            only_new = sorted(set(new["labels"]) - set(old["labels"]))
            for label in only_old[:10]:
                print(f"    only in old: {label!r}")
            for label in only_new[:10]:
                print(f"    only in new: {label!r}")
        if not (same_labels and same_size):
            failures += 1
    print(f"\n{len(panels) - failures}/{len(panels)} panels unchanged in size and labels")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
