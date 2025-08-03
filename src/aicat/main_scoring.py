from .agent import CellTypeAgent, AnnotationEvaluation
from .utils import format_query_prompt, chat_to_markdown
import os
import pandas as pd
import json

def scoring_annotation(openai_api_key,
                       tissue,
                       save_path_root,
                       AnnoSingle_res_path):
    
    res_tab = pd.read_csv(AnnoSingle_res_path)

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