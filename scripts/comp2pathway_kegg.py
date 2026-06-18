import requests

def get_kegg_compound_info(kegg_compound_id):
    """Fetches KEGG compound info, including pathway links."""
    url = f"https://rest.kegg.jp/get/{kegg_compound_id}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    return None

def parse_kegg_pathways(kegg_text):
    """Parses KEGG text output to extract pathway information."""
    pathways = []
    in_pathway_section = False
    for line in kegg_text.splitlines():
        if line.startswith("PATHWAY"):
            in_pathway_section = True
            parts = line.split()
            if len(parts) > 1:
                pathway_id = parts[1]
                pathway_name = " ".join(parts[2:])
                pathways.append({"id": pathway_id, "name": pathway_name})
        elif in_pathway_section and not line.startswith("  "): # Stop if not indented
            in_pathway_section = False
    return pathways

# --- Main logic ---
# Assuming you have a list of KEGG Compound IDs from HMDB
# Example: HMDB00001 maps to C00001, C00002 etc. in KEGG
# For demonstration, let's assume we have one KEGG ID
kegg_id = "C00001" # Example: Acetate
kegg_data = get_kegg_compound_info(kegg_id)

if kegg_data:
    pathway_info = parse_kegg_pathways(kegg_data)
    print(f"KEGG ID: {kegg_id}")
    print("Associated Pathways:")
    for pw in pathway_info:
        print(f"  - ID: {pw['id']}, Name: {pw['name']}")
else:
    print(f"Could not retrieve data for KEGG ID: {kegg_id}")

# You would repeat this for all KEGG Compound IDs obtained from HMDB
# and then combine this with your HMDB IDs.