from langchain_openai import ChatOpenAI
from langchain.agents import create_tool_calling_agent, AgentExecutor
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.output_parsers import PydanticOutputParser
from pydantic import BaseModel
from typing import Optional, List, Dict
import re

class Annotate(BaseModel):
    yes_or_no: str
    cell_lineage: List[str]
    cell_type: List[str]
    summary: Optional[str] = ""
    cell_lineage_marker_map: Dict[str, List[str]]
    cell_type_marker_map: Dict[str, List[str]]
    sources: Optional[List[str]] = []
    tools_used: Optional[List[str]] = []

class single_pair_res(BaseModel):
    group_ids: List[str]
    group: List[str]
    merged_name: str
    reason: str

class group_clusters(BaseModel):
    results: List[single_pair_res]

class cross_Annotate_single(BaseModel):
    yes_or_no: str
    input_annotation: str
    summary: Optional[str] = ""
    renewed_annotation: Optional[str] = ""
    renewed_annotation_marker_map: Optional[Dict[str, List[str]]] = {}
    sources: Optional[List[str]] = []
    tools_used: Optional[List[str]] = []

class cross_Annotate_pair(BaseModel):
    merge_yes_or_no: str
    target_cell_type: Optional[str] = ""
    query_cell_type: Optional[str] = ""
    summary: Optional[str] = ""
    merged_cell_type_name: Optional[str] = ""
    merged_cell_type_marker_map: Optional[Dict[str, List[str]]] = {}
    sources: Optional[List[str]] = []
    tools_used: Optional[List[str]] = []

class AnnotationEvaluation(BaseModel):
    lineage_score: Optional[float] = None
    cell_type_score: Optional[float] = None
    functional_score: Optional[float] = None

    final_score: float

    lineage_justification: Optional[str] = None
    cell_type_justification: Optional[str] = None
    functional_justification: Optional[str] = None
    final_justification: str

    sources: Optional[List[str]] = []
    tools_used: Optional[List[str]] = []


class CellTypeAgent:
    def __init__(self, api_key, tools, 
                 ResponseFormat=Annotate,
                 model_name="gpt-4o-mini", 
                 verbose=True, 
                 mode="single"):
        self.api_key = api_key
        self.tools = tools
        self.ResponseFormat = ResponseFormat
        self.model_name = model_name
        self.verbose = verbose
        self.mode = mode
        self.chat_history = [] # store (role, content) tuples or Message objects
        self._initialize_agent()

    def _initialize_agent(self):
        self.llm = ChatOpenAI(model=self.model_name, 
                              api_key=self.api_key,
                              temperature=1.0 # default 1.0; lowering the temperture will make the model more deterministic
                              )
        self.parser = PydanticOutputParser(pydantic_object=self.ResponseFormat)

        if self.mode == "single":
            system_prompt = (
                "You are an expert in cell type annotation. "
                "**First, try to answer the user's question using your own knowledge. "
                "Provide source link/paper.** "
                "**Use tools if more information is needed.** "
                "Focus on decisive marker genes for cell lineage and type identification. "
                "If there is no definite marker gene, use majority vote to decide cell lineages and types. "
                "Wrap the output in this *strict* format \n{format_instructions}\n "
                "Each cell type/lineage must have its own list of marker genes in the "
                "'cell_type_marker_map' and 'cell_lineage_marker_map' field. "
                "You must return a **strict JSON object** that includes *all required fields*, even if values are null."
            )
        elif self.mode == "group_clusters":
            system_prompt = (
                "You are an expert in cell type annotation.\n"
                "Given a list of cell type names, your task is to identify which cell types in the list refer to the same biological populations.\n"
                "Group highly similar cell types together if they are alternate names.\n"
                "Return the result as a JSON list. Each group should include:\n"
                "- \"group_ids\": the list of IDs from the input that belong to the same group\n"
                "- \"group\": the list of input names that refer to similar cell types\n"
                "- \"merged_name\": a suggested unified name\n"
                "- \"reason\": a brief explanation for the grouping\n\n"
                "Only consider relationships **within the provided list**. Ignore synonyms not present in the input.\n"
                "Return strictly formatted JSON like \n{format_instructions}\n, no extra commentary.\n\n"
                "Example input:\n- cluster_1: Basal keratinocyte\n- cluster_2: Keratinocyte (basal)\n- cluster_3: T cell\n- cluster_4: CD8+ T cell\n\n"
                "Example output:\n"
                "[\n"
                "  {{\n"
                "    \"group_ids\": [\"cluster_1\", \"cluster_2\"],\n"
                "    \"group\": [\"Basal keratinocyte\", \"Keratinocyte (basal)\"],\n"
                "    \"merged_name\": \"Basal keratinocyte\",\n"
                "    \"reason\": \"Both refer to basal-layer keratinocytes.\"\n"
                "  }}\n"
                "]"
            )
        elif self.mode == "cross":
            system_prompt = (
                "You are an expert in cell type annotation. "
                "**Begin by answering the user's question using your own knowledge whenever possible. "
                "Cite relevant sources such as papers or databases.** "
                "**If additional information is needed, use the available tools.** "
                "Be aware that some cell types may share similar marker genes—distinguish them carefully. "
                "Base your judgment on the specific differences between marker gene profiles. "
                "Wrap the output in this *strict* format \n{format_instructions}\n. "
                "For the 'renewed_annotation' field, provide a concise but informative name. "
                "You must return a **strict JSON object** that includes *all required fields*, even if values are null."
            )
        elif self.mode == "assess":
            system_prompt = ("""
                You are a scientific annotation expert specializing in immunology and single-cell biology. 
                Your task is to assess how well a set of annotations (lineage, cell type, and functional annotation) 
                match a known ground truth cell type.

                Annotations may include one, two, or all three levels.

                Use the following scoring system for each available level:
                - 1 = correct (matches or is a subtype of the ground truth)
                - 0.5 = partially correct (biologically related but too broad or imprecise)
                - 0 = incorrect (biologically unrelated)
                - -1 = seriously wrong (contradicts known lineage or function)

                Rules for computing the final score:
                - If any provided level scores 1 → final score = 1
                - Else → final score = highest of the available scores
                - If no annotations are provided, final score = -1

                For each available level, provide:
                - A numeric score
                - A brief biological justification

                Then provide:
                - final_score
                - final_justification
                - sources: optional list of data sources or ontologies used
                - tools_used: optional list of tools, APIs, or packages used
                Wrap the output in this *strict* format \n{format_instructions}\n.
                """
            )

        self.prompt = ChatPromptTemplate.from_messages(
            [
                ("system", system_prompt),
                ("placeholder", "{chat_history}"),
                ("human", "{query}"),
                ("placeholder", "{agent_scratchpad}")
            ]
        ).partial(format_instructions=self.parser.get_format_instructions())

        agent = create_tool_calling_agent(
            llm=self.llm,
            prompt=self.prompt,
            tools=self.tools,
        )

        self.agent_executor = AgentExecutor(agent=agent, 
                                            tools=self.tools, 
                                            verbose=self.verbose)
    
    def update_ResponseFormat(self, ResponseFormat):
        self.ResponseFormat = ResponseFormat
        self.parser = PydanticOutputParser(pydantic_object=self.ResponseFormat)

        # Update the prompt with new format instructions
        self.prompt = self.prompt.partial(format_instructions=self.parser.get_format_instructions())

        # Rebuild agent with updated prompt
        agent = create_tool_calling_agent(
            llm=self.llm,
            prompt=self.prompt,
            tools=self.tools,
        )
        self.agent_executor = AgentExecutor(agent=agent,
                                        tools=self.tools,
                                        verbose=self.verbose)

    def run(self, query):
        input_data = {
            "query": query,
            "chat_history": self.chat_history,
        }
        # Update chat history with the current query
        response = self.agent_executor.invoke(input_data)
        # Update chat history with the response
        self.chat_history.append(("user", query))
        self.chat_history.append(("ai", response["output"]))
        return response

    def parse_output(self, response):
        output = response["output"]
        json_str = self._extract_json_from_response(output)
        parsed_res = self.ResponseFormat.model_validate_json(json_str)
        return parsed_res
    
    def _extract_json_from_response(self, response_text):
        # Match content between ```json ... ```
        match = re.search(r"```json\s*(\{.*?\})\s*```", response_text, re.DOTALL)
        
        if match:
            json_str = match.group(1)
        else:
            # Try fallback: if no code block, maybe it's already plain JSON
            return response_text.strip()
        
        # Clean control characters (ASCII 0–31 except \t, \n, \r)
        json_str_cleaned = re.sub(r'[\x00-\x08\x0B\x0C\x0E-\x1F]', '', json_str)
        
        return json_str_cleaned
