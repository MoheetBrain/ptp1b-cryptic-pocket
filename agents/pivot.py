import pandas as pd
from chembl_webresource_client.new_client import new_client

class PivotAgent:
    def __init__(self):
        print("🔄 [PIVOT] Initializing Allosteric Hunter...")
        self.mechanism = new_client.mechanism
        self.assay = new_client.assay
        self.activity = new_client.activity

    def hunt_allosteric(self):
        print("   ...Scanning ChEMBL for keyword 'Allosteric' on PTP1B...")
        
        # 1. Search for Mechanisms explicitly labeled "Allosteric"
        # Target CHEMBL335 is PTP1B
        mechs = self.mechanism.filter(
            target_chembl_id="CHEMBL335", 
            mechanism_of_action__icontains="allosteric"
        )
        
        candidates = []
        seen_ids = set()

        # Process Mechanism Hits
        print(f"   ...Found {len(mechs)} confirmed allosteric records via Mechanism.")
        for m in mechs:
            cid = m['molecule_chembl_id']
            if cid not in seen_ids:
                candidates.append({
                    'chembl_id': cid,
                    'reason': 'Confirmed Mechanism: ' + m['mechanism_of_action'],
                    'source': 'Mechanism API'
                })
                seen_ids.add(cid)

        # 2. Search for Assays containing "Allosteric" (Deep Dig)
        # This finds compounds tested in assays designed for allosteric binding
        assays = self.assay.filter(
            target_chembl_id="CHEMBL335", 
            description__icontains="allosteric"
        )
        
        print(f"   ...Found {len(assays)} assays mentioning 'Allosteric'. Digging for compounds...")
        
        for i, assay in enumerate(assays[:5]): # Check first 5 relevant assays to save time
            assay_id = assay['assay_chembl_id']
            acts = self.activity.filter(assay_chembl_id=assay_id, standard_type="IC50", standard_value__lte=10000)
            
            for act in acts:
                cid = act['molecule_chembl_id']
                if cid not in seen_ids:
                    candidates.append({
                        'chembl_id': cid,
                        'reason': f"Found in Allosteric Assay: {assay['description'][:50]}...",
                        'source': 'Assay Description'
                    })
                    seen_ids.add(cid)

        # 3. Report
        if candidates:
            df = pd.DataFrame(candidates)
            print(f"\n✅ [PIVOT] Found {len(df)} TRUE Allosteric Candidates.")
            print(df.head(10).to_string(index=False))
            
            # Save for the Architect
            df.to_csv("data/processed/ptp1b_true_allosteric.csv", index=False)
        else:
            print("❌ No allosteric candidates found.")

if __name__ == "__main__":
    bot = PivotAgent()
    bot.hunt_allosteric()
    