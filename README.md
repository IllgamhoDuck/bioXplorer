# bioXplorer
exploratory data analysis of cell data using cellxgene. Everything works based on anndata

## Installation
```python
git clone https://github.com/IllgamhoDuck/bioXplorer.git
pip install ./bioXplorer
```

## How to use?

### Pseudo-bulk Creation
```python
from bioxplorer.pseudo_bulk import convert_to_pseudo_bulk

pdata = convert_to_pseudo_bulk(adata)
```

### Embedding Extraction
> Environment setting (downloading, installation etc) will automatically processed. only support 3 model currently

```python
from bioxplorer.embeddings import Embeddings

emb = Embeddings()
uce_data = emb.extract('uce', adata)
scgpt_adata = emb.extract('scgpt', adata)
geneformer_adata = emb.extract('geneformer', adata)
```

