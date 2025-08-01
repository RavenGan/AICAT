import os
import argparse
# os.chdir("/Users/david/Desktop/CelltypeAgenticAI")
os.chdir("/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/CelltypeAgenticAI")
import pandas as pd
import sys
import json
# sys.path.append(os.path.abspath("/Users/david/Desktop/CelltypeAgenticAI/code/Build_Agent"))
sys.path.append(os.path.abspath("/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/CelltypeAgenticAI/code/Build_Agent"))
from agent import CellTypeAgent, AnnotationEvaluation
from utils import format_query_prompt, chat_to_markdown

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--tissue", type=str, help="A chosen tissue from the data MCA")
parser.add_argument("--save_path_root", type=str, required=True)
args = parser.parse_args()

# set up the OpenAI API key CONFIDENTIAL
openai_api_key = "OpenAI_API_Key_Confidential"

tissue = args.tissue
save_path_root = args.save_path_root

res_path = f'./res/2025_0712_MCA_AICAT/CollectAnnoSingle_{tissue.lower()}_res_df.csv'

res_tab = pd.read_csv(res_path)

# Format prompts
query_ls = []
for idx, row in res_tab.iterrows():
    query = format_query_prompt(row)
    query_ls.append(query)

# Set up the agent
agent = CellTypeAgent(api_key=openai_api_key, 
                        tools = [], #[wiki_tool],
                        ResponseFormat=AnnotationEvaluation,
                        verbose=False,
                        mode="assess")

res_dict_res = {}
final_score_ls = []

for idx, query in enumerate(query_ls):
    try:
        print(f"Processing query {idx + 1}/{len(query_ls)}")
        response = agent.run(query)
        response_parsed = agent.parse_output(response)
        
        # Convert Annotate object to dict for serialization
        res_dict = response_parsed.model_dump()
        res_dict_res[idx] = res_dict
        
        final_score = res_dict.get('final_score', None)
        final_score_ls.append(final_score)

    except Exception as e:
        print(f"Error processing query {idx + 1}: {e}")
        res_dict_res[idx] = {"error": str(e)}
        final_score_ls.append(None)

# organize and save the chat history
chat_hist = agent.chat_history
chat_hist = chat_to_markdown(chat_hist)

save_folder = f"{save_path_root}/{tissue}"
os.makedirs(save_folder, exist_ok=True)

res_tab['final_score'] = final_score_ls
res_tab.to_csv(f"{save_folder}/CollectAnnoSingle_{tissue}_WithScore.csv", index=False)

with open(f"{save_folder}/{tissue}_chat_history.md", "w") as f:
     f.write(chat_hist)

with open(f"{save_folder}/{tissue}_res_dict.json", "w") as f:
        json.dump(res_dict_res, f, indent=4)


