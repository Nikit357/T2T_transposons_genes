"""Extract the FINAL text (all tracked changes accepted) of a .docx, one paragraph per line.

Usage: python extract_docx.py <in.docx> <out.txt>
"""
import sys, zipfile, re
from lxml import etree

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def para_final_text(p):
    """Text of a paragraph as it would read with all tracked changes accepted."""
    parts = []
    for node in p.iter():
        tag = node.tag
        if tag == W + "t":
            # skip text inside a <w:del> (deleted -> gone when accepted)
            anc = node.getparent()
            deleted = False
            while anc is not None:
                if anc.tag == W + "del":
                    deleted = True
                    break
                if anc.tag == W + "p":
                    break
                anc = anc.getparent()
            if not deleted:
                parts.append(node.text or "")
        elif tag == W + "delText":
            pass
        elif tag == W + "tab":
            parts.append("\t")
        elif tag == W + "br":
            parts.append(" ")
    return "".join(parts)


def para_original_text(p):
    """Text as it would read with all tracked changes REJECTED."""
    parts = []
    for node in p.iter():
        tag = node.tag
        if tag in (W + "t", W + "delText"):
            anc = node.getparent()
            inserted = False
            while anc is not None:
                if anc.tag == W + "ins":
                    inserted = True
                    break
                if anc.tag == W + "p":
                    break
                anc = anc.getparent()
            if not inserted:
                parts.append(node.text or "")
        elif tag == W + "tab":
            parts.append("\t")
    return "".join(parts)


def main():
    src, dst = sys.argv[1], sys.argv[2]
    mode = sys.argv[3] if len(sys.argv) > 3 else "final"
    fn = para_final_text if mode == "final" else para_original_text
    with zipfile.ZipFile(src) as z:
        xml = z.read("word/document.xml")
    root = etree.fromstring(xml)
    body = root.find(W + "body")
    out = []
    for el in body.iter(W + "p"):
        t = fn(el).strip()
        if t:
            out.append(t)
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"{src}: {len(out)} non-empty paragraphs -> {dst}")


main()
