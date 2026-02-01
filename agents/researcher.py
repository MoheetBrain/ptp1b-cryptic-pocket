import pandas as pd
from chembl_webresource_client.new_client import new_client

# --- CONFIG ---
INPUT_FILE = "data/processed/ptp1b_architect_designs.csv"

class ResearcherAgent:
    def __init__(self):
        print("🕵️‍♂️ [RESEARCHER] Initializing Intel Unit...")
        self.document = new_client.document
        self.activity = new_client.activity

    def get_paper_info(self, chembl_id):
        """Finds the scientific paper associated with a compound."""
        try:
            # 1. Get activities for this molecule to find the Document ID
            activities = self.activity.filter(molecule_chembl_id=chembl_id, target_chembl_id="CHEMBL335")
            if not activities:
                return None
            
            # Take the first activity record (usually the primary paper)
            doc_id = activities[0]['document_chembl_id']
            
            # 2. Fetch Document Metadata
            doc = self.document.get(doc_id)
            return {
                'title': doc.get('title', 'Unknown Title'),
                'journal': doc.get('journal', 'Unknown Journal'),
                'year': doc.get('year', 'N/A'),
                'abstract': doc.get('abstract', 'No abstract available.'),
                'doi': doc.get('doi', 'N/A')
            }
        except Exception as e:
            print(f"   ⚠️ Error fetching intel for {chembl_id}: {e}")
            return None

    def execute_research_mission(self):
        print(f"📂 [RESEARCHER] Reading Top Candidates from {INPUT_FILE}...")
        
        try:
            df = pd.read_csv(INPUT_FILE)
        except FileNotFoundError:
            print("❌ File not found.")
            return

        # Take the Top 3 "Gold" Candidates
        top_hits = df.head(3)
        
        print("\n🔎 --- INTELLIGENCE REPORT ---")
        for index, row in top_hits.iterrows():
            cid = row['chembl_id']
            ic50 = row['ic50_nm']
            
            print(f"\n🧪 Target: {cid} (IC50: {ic50} nM)")
            print("   ...Tracing source paper...")
            
            paper = self.get_paper_info(cid)
            
            if paper:
                print(f"   📄 Title:   {paper['title']}")
                print(f"   📅 Journal: {paper['journal']} ({paper['year']})")
                print(f"   🔗 DOI:     https://doi.org/{paper['doi']}")
                # print(f"   📝 Abstract: {paper['abstract'][:200]}...") # Uncomment for full text
            else:
                print("   ❌ Source Classified/Unknown.")

if __name__ == "__main__":
    agent = ResearcherAgent()
    agent.execute_research_mission()