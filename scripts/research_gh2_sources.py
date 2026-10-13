"""Cache public source records for selection of a GH2 beta-galactosidase dataset."""

import concurrent.futures
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "raw" / "research"


def fetch(item):
    name, url = item
    path = CACHE / name
    if not path.exists():
        for attempt in range(3):
            try:
                request = urllib.request.Request(url, headers={"User-Agent": "GH2Research/1.0 (educational protein sequence collection)"})
                with urllib.request.urlopen(request, timeout=45) as response:
                    body = response.read()
                path.write_bytes(body)
                break
            except Exception as error:
                if attempt == 2:
                    return {"file": name, "url": url, "error": str(error)}
                time.sleep(2 * (attempt + 1))
    return {"file": name, "url": url, "bytes": path.stat().st_size}


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    CACHE.mkdir(parents=True, exist_ok=True)
    sources = [("cazy_gh2.html", "https://www.cazy.org/GH2_characterized.html")]
    sources.append(("ACI41243.1.fasta", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=ACI41243.1&rettype=fasta&retmode=text"))
    sources.append(("ACI41243.1.gp", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=ACI41243.1&rettype=gp&retmode=text"))
    for pmc in ["PMC7832830", "PMC11677235", "PMC3268748", "PMC6315854", "PMC193751", "PMC91242"]:
        sources.append((pmc + ".xml", f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML"))
    for accession in ["P81650", "P00722", "Q59140", "Q8KRF6", "Q08KN3", "D9IU58", "A0A023UGN9", "B6EAX7", "A1SWB8", "Q6LL68"]:
        sources.append((accession + ".json", f"https://rest.uniprot.org/uniprotkb/{accession}.json"))
    queries = {
        "cold_gh2_catalytic": '(protein_name:"beta-galactosidase") AND (organism_name:"Arthrobacter" OR organism_name:"Pseudoalteromonas" OR organism_name:"Paracoccus") AND (xref:pfam-PF02836)',
        "meso_reviewed": '(protein_name:"beta-galactosidase") AND (xref:pfam-PF00703) AND (reviewed:true)',
    }
    for name, query in queries.items():
        url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode({"query": query, "format": "json", "size": 500})
        sources.append((name + ".json", url))
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        results = list(pool.map(fetch, sources))
    (CACHE / "fetch_log.json").write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    for result in results:
        print(json.dumps(result, ensure_ascii=False))
    for name in queries:
        path = CACHE / (name + ".json")
        if path.exists():
            entries = json.loads(path.read_text(encoding="utf-8")).get("results", [])
            print(name, "records:", len(entries))
            if name.startswith("cold"):
                for entry in entries:
                    if "sp." in entry["organism"]["scientificName"] or "psychrolactophilus" in entry["organism"]["scientificName"]:
                        print(name, entry["primaryAccession"], entry["organism"]["scientificName"], entry["sequence"]["length"])


if __name__ == "__main__":
    main()
