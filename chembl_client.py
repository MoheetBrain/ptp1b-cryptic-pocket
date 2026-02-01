import os
import time
import pandas as pd
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors

# --- CONFIGURATION ---
TARGET_ID = "CHEMBL335"  # PTP1B
MAX_RESULTS = 1000       # Cap for safety, increase later
OUTPUT_DIR = "data/processed"
OUTPUT_FILE = "ptp1b_high_confidence.csv"

def fetch_and_clean_data():
    print(f"🔧 [MINER] Connecting to ChEMBL for Target {TARGET_ID}...")
    
    activity = new_client.activity
    molecule = new_client.molecule
    
    # 1. THE FILTER (The Perplexity Guardrails)
    # We filter server-side as much as possible to save time
    print(f"🔍 [MINER] Querying for IC50 < 1000nM, High Confidence, Binding/Functional assays...")
    activities = activity.filter(
        target_chembl_id=TARGET_ID,
        standard_type="IC50",
        standard_units="nM",           # STRICT: Nanomolar only
        confidence_score__gte=7,       # STRICT: High confidence data only
        assay_type__in=['B', 'F']      # Binding or Functional assays only
    ).filter(
        standard_value__lte=1000       # Potency limit (1 uM)
    ).only(
        'molecule_chembl_id', 'standard_value', 'standard_units', 
        'standard_relation', 'type'
    )

    # 2. LOCAL PROCESSING & DEDUPLICATION
    unique_hits = {}
    print("⏳ [MINER] Downloading and deduplicating (keeping best IC50 per molecule)...")
    
    count = 0
    for act in activities:
        if len(unique_hits) >= MAX_RESULTS:
            break
            
        chembl_id = act.get('molecule_chembl_id')
        val = act.get('standard_value')
        
        if not chembl_id or not val:
            continue
            
        try:
            ic50 = float(val)
        except (ValueError, TypeError):
            continue
            
        # Logic: If we haven't seen this ID, OR if this new IC50 is lower (better) than what we have, take it.
        if chembl_id not in unique_hits or ic50 < unique_hits[chembl_id]['ic50_nm']:
            unique_hits[chembl_id] = {
                'chembl_id': chembl_id,
                'ic50_nm': ic50,
                'units': act.get('standard_units'),
                'relation': act.get('standard_relation'),
                'smiles': None # Placeholder
            }
            
        count += 1
        if count % 100 == 0:
            print(f"   ...scanned {count} records. Unique hits so far: {len(unique_hits)}")

    print(f"💎 [MINER] Retrieved {len(unique_hits)} unique high-quality hits. Fetching SMILES...")

    # 3. ENRICH WITH STRUCTURES (SMILES)
    final_data = []
    # Convert dict to list for processing
    hits_list = list(unique_hits.values())
    
    for i, item in enumerate(hits_list):
        try:
            # Respect API rate limits gently
            if i % 10 == 0: time.sleep(0.1)
            
            mol_data = molecule.get(item['chembl_id'])
            if mol_data and mol_data.get('molecule_structures'):
                smiles = mol_data['molecule_structures'].get('canonical_smiles')
                if smiles:
                    item['smiles'] = smiles
                    
                    # 4. RDKit FIRST PASS (Add MW/LogP immediately)
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        item['mw'] = Descriptors.MolWt(mol)
                        item['logp'] = Descriptors.MolLogP(mol)
                        final_data.append(item)
                        
        except Exception as e:
            print(f"⚠️ Error processing {item['chembl_id']}: {e}")
            continue

    return pd.DataFrame(final_data)

def save_data(df):
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        print(f"📁 Created directory: {OUTPUT_DIR}")
        
    full_path = os.path.join(OUTPUT_DIR, OUTPUT_FILE)
    df.to_csv(full_path, index=False)
    print(f"✅ [SUCCESS] Data saved to: {full_path}")
    print("\n--- PREVIEW ---")
    print(df.head())
    print(f"Total Molecules: {len(df)}")

if __name__ == "__main__":
    df = fetch_and_clean_data()
    if not df.empty:
        save_data(df)
    else:
        print("❌ [FAILURE] No data found.")