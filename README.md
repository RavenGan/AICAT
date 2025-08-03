# AICAT: Agentic-AI Cell-type Annotation Tool

**AICAT** (Agentic-AI Cell-type Annotation Tool) is a computational method for automatic cell type annotation of single-cell RNA-seq data. It leverages the model via OpenAI's API to generate deep and context-aware cell type predictions based on cluster-level expression profiles.

---

## Installation

To install the package in development mode:

```bash
git clone https://github.com/RavenGan/AICAT.git
cd AICAT
pip install -e .
```

## API Key Requirement
AICAT requires access to the OpenAI GPT model. Users need to provide their own [OpenAI API](https://platform.openai.com/api-keys) key. To avoid the risk of exposing the API key or committing the key to browsers, users need to set up the API key as a system environment variable before running AICAT. 

Users can generate the API key in the [OpenAI](https://openai.com/) account webpage: log in to OpenAI. In the pop-up windows, click on “->” next to “API”; next, click on the left-hand-side icon of “API key”; then click on “Create new secret key” to create your key which directs you to the API key page. Copy the key and paste it on a note for further use. Avoid sharing your API key with others or uploading it to public spaces. Make sure it’s not visible in browsers or any client-side scripts. Finally, on the left bar, click “Settings”; on the break-down list click on “Billing”, and make sure you have non-zero credit balance.

With the API key, there are two options to setup the key:
#### Option 1: Using `.env` file
```bash
OPENAI_API_KEY=your-api-key-here
```
#### Option 2: Using terminal (temporary)
```bash
export OPENAI_API_KEY=your-api-key-here
```

## Supported Input
The main input is a clustered `.h5ad` file containing single-cell gene expression data.

#### Required Structure
- The `.obs` DataFrame should contain a cluster label column.
- The expression matrix should be preprocessed (e.g., log-normalized) to allow differential gene calculation later in the program.
- The `.h5ad` file can be generated via the python package `Scanpy`.

## Usage
### CLI Usage
#### Indepth annotation CLI
Run the following from your terminal. Here a small test data is used as an example. This data can also be found under the folder `./tests/data`.
```bash
aicat-indepth \
  --openai_api_key "OPENAI_API_KEY" \
  --adata_path "tests/data/CRC_SMC05-T_processed.h5ad" \
  --species "human" \
  --tissue "primary colorectal cancer" \
  --cluster_col_name "Cell_type" \
  --data_name "CRC_SMC05-T" \
  --save_path "tests/res/CRC_SMC05-T" 
```

Display CLI options with 
```bash
aicat-indepth --help
```

#### Subcluster CLI (optional)
To perform subclustering annotation, run the following from the terminal. The argument `--AnnoSingle_res_path` requires the output from the previous command `aicat-indepth`.
```bash
aicat-subcluster \
  --openai_api_key "OPENAI_API_KEY" \
  --adata_path "tests/data/CRC_SMC05-T_processed.h5ad" \
  --tissue "primary colorectal cancer" \
  --cluster_col_name 'Cell_type' \
  --chosen_cluster "B cells" \
  --AnnoSingle_res_path "tests/res/CRC_SMC05-T/AnnoSingle_primary colorectal cancer_res_dict.json" \
  --save_path "tests/res/CRC_SMC05-T_subcluster"
```
Display CLI options with 
```bash
aicat-subcluster --help
```

### Programmatic Usage
Users can also call `aicat` from python:
```bash
from aicat.main_indepth import indepth_annotation
# Indepth annotation
indepth_annotation(api_key, 
                    adata_path,
                    species,
                    tissue,
                    cluster_col_name)

# subclustering annotation
subcluster_annotation(api_key, 
                        adata_path,
                        tissue,
                        cluster_col_name,
                        chosen_cluster,
                        AnnoSingle_res_path)
```

## Output
AICAT will output the following:
- A `.json` file with predicted cell types
- The full conversation history including the prompts and GPT response content
- A log of the full running progress

## Testing
Users can also run unit tests using 
```bash
pytest tests/
```
To test a specific function manually:
```bash
python tests/test_indepth.py
```
Make sure `OPENAI_API_KEY` is set in your environment.

## Dependencies
Major dependencies used in this package:
- `langchain`, `langchain-openai`, `langchain-core`
- `pydantic`, `pydantic-settings`
- `scanpy`, `pandas`, `dotenv`, `requests`
- `argparse`, `json`, `os`, `re`, `subprocess`
See `pyproject.toml` for full details about the package versions.

## License
This project is licensed under the MIT License.

## Contributions

Contributions are welcome! If you’d like to add features or fix bugs, feel free to fork the repository and open a pull request.