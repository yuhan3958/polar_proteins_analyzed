"""Print selected facts from locally cached public records, without network access."""

import html
import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

CACHE = Path(__file__).resolve().parents[1] / "data" / "raw" / "research"
sys.stdout.reconfigure(encoding="utf-8")

for name in ["cold_gh2_catalytic", "meso_reviewed"]:
    path = CACHE / (name + ".json")
    if not path.exists():
        continue
    entries = json.loads(path.read_text(encoding="utf-8"))["results"]
    for e in entries:
        organism = e["organism"]["scientificName"]
        if name == "cold_gh2_catalytic" and not ("sp." in organism and e["sequence"]["length"] >= 950):
            continue
        if name == "meso_reviewed" and not any(term in organism for term in ["agglomerans", "cloacae", "ATCC BAA-895", "ATCC 39541", "ATCC 7966", "IL1403", "xylosus", "adolescentis", "megaterium", "Klebsiella pneumoniae"]):
            continue
        pfam = {x["id"]: x.get("properties", [{}])[0].get("value", "") for x in e.get("uniProtKBCrossReferences", []) if x["database"] == "Pfam"}
        print(e["primaryAccession"], organism, e["sequence"]["length"], pfam, "bad", re.sub("[ACDEFGHIKLMNPQRSTVWY]", "", e["sequence"]["value"]))

body = (CACHE / "cazy_gh2.html").read_text(encoding="utf-8")
for row in re.findall(r"<tr\b.*?</tr>", body, flags=re.S | re.I):
    if any(term in row for term in ["Arthrobacter", "P81650", "Q6RGS1", "A1SWB8", "Q6LL68"]):
        print("CAZY:", re.sub(r"\s+", " ", html.unescape(re.sub("<[^>]+>", " ", row))).strip())

for pmc in ["PMC193751", "PMC91242"]:
    path = CACHE / (pmc + ".xml")
    if path.exists():
        root = ET.parse(path).getroot()
        for p in root.iter("p"):
            text = "".join(p.itertext())
            if re.search(r"grew at|optimal growth|optimal temperature|accession", text, re.I):
                print(pmc, text)
