

import os
import json
import subprocess
import warnings
from scipy.sparse import csr_matrix

from bioxplorer.environment import TEMP_DIR
from bioxplorer.environment import install_uce
from bioxplorer.environment import install_geneformer
from bioxplorer.environment import install_scgpt


warnings.filterwarnings("ignore")


class Embeddings:
    _initialized = None

    def __new__(cls):
        if not hasattr(cls, '_instance'):
            cls._instance = super(Embeddings, cls).__new__(cls)
        return cls._instance

    def __init__(self) -> None:
        if self._initialized:
            return
        self._initialized = True

        self.install_fns = {
            'uce': install_uce,
            'scgpt': install_scgpt,
            'geneformer': install_geneformer
        }
        self.setup_fns = {
            'uce': setup_uce,
            'scgpt': setup_scgpt,
            'geneformer': setup_geneformer,
        }
        self.is_initialized = {fn: False for fn in self.setup_fns}
        self.extract_fns = {
            'uce': extract_uce,
            'scgpt': extract_scgpt,
            'geneformer': extract_geneformer
        }

        # where model store necessary information
        self.cache = {fn: None for fn in self.setup_fns}

    def __repr__(self):
        return "Embeddings()"

    def install(self, model):
        """
        install the model only when it is used
        """
        if model in self.install_fns:
            self.install_fns[model]()
        else:
            raise ValueError(f"Model {model} not supported.")

    def setup(self, model):
        if model in self.setup_fns:
            return self.setup_fns[model]()
        else:
            raise ValueError(f"Model {model} not supported.")


    def extract(self, model, adata, verbose=False, *args, **kwargs):
        if not self.is_initialized[model]:
            self.install(model)
            self.cache[model] = self.setup(model)
            self.is_initialized[model] = True

        if model in self.extract_fns:
            return self.extract_fns[model](
                adata,
                cache=self.cache[model],
                verbose=verbose,
                *args,
                **kwargs,
            )
        else:
            raise ValueError(f"Model {model} not supported.")


def setup_uce():
    return None

def setup_scgpt():
    model_dir = os.path.join(TEMP_DIR, "scgpt", "model")
    return {
        'model_dir': model_dir,
    }

def setup_geneformer():
    from geneformer import (
        EmbExtractor,
        TranscriptomeTokenizer,
    )

    tokenizer = TranscriptomeTokenizer(custom_attr_name_dict={"joinid": "joinid"})
    model_dir = os.path.join(TEMP_DIR, "geneformer", "models", "fine_tuned_geneformer")
    label_mapping_dict_file = os.path.join(model_dir, "label_to_cell_subclass.json")

    with open(label_mapping_dict_file) as fp:
        label_mapping_dict = json.load(fp)

    n_classes = len(label_mapping_dict)
    emb_extractor = EmbExtractor(
        model_type="CellClassifier",
        num_classes=n_classes,
        max_ncells=None,
        emb_label=["joinid"],
        emb_layer=0,
        forward_batch_size=30,
        nproc=8,
    )
    return {
        'tokenizer': tokenizer,
        'model_dir': model_dir,
        'emb_extractor': emb_extractor,
    }


def extract_uce(adata, cache=None, verbose=False, *args, **kwargs):
    adata = adata.copy()
    return adata


# check if scGPT is already installed
IS_SCGPT_IMPORTED = False

def extract_scgpt(adata, gene_col="feature_name", cache=None, verbose=False, *args, **kwargs):
    global IS_SCGPT_IMPORTED
    if not IS_SCGPT_IMPORTED:
        import scgpt
        IS_SCGPT_IMPORTED = True
    
    adata = adata.copy()
    adata = scgpt.tasks.embed_data(
        adata,
        cache['model_dir'],
        gene_col=gene_col,
        batch_size=64,
    )

    return adata

def extract_geneformer(
    adata,
    ensembl_id_col='feature_id',
    cache=None,
    verbose=False,
    *args,
    **kwargs,
):
    adata = adata.copy()

    # preprocess for Geneformer
    adata.var['ensembl_id'] = adata.var[ensembl_id_col]
    adata.obs['n_counts'] = adata.X.sum(axis=1)
    adata.obs["joinid"] = list(range(adata.n_obs))
    adata.var.index = adata.var.ensembl_id

    # prepare huggingface dataset
    h5ad_dir = os.path.join(TEMP_DIR, "geneformer", "h5ad")
    token_dir = os.path.join(TEMP_DIR, "geneformer", "tokenized_data")
    embs_dir = os.path.join(TEMP_DIR, "geneformer", "embeddings")

    if os.path.exists(h5ad_dir):
        subprocess.run(["rm", "-r", h5ad_dir])
    if os.path.exists(token_dir):
        subprocess.run(["rm", "-r", token_dir])
    if os.path.exists(embs_dir):
        subprocess.run(["rm", "-r", embs_dir])

    os.makedirs(h5ad_dir)
    os.makedirs(token_dir)
    os.makedirs(embs_dir)

    # dense matrix in h5ad format will cause error for fancy indexing
    if isinstance(adata.X, list):
        adata.X = csr_matrix(adata.X)

    adata.write(os.path.join(h5ad_dir, "geneformer.h5ad"))
    cache['tokenizer'].tokenize_data(
        data_directory=h5ad_dir,
        output_directory=token_dir,
        output_prefix="geneformer",
        file_format="h5ad",
    )

    # extract embeddings
    embs = cache['emb_extractor'].extract_embs(
        model_directory=cache['model_dir'],
        input_data_file=os.path.join(token_dir, "geneformer.dataset"),
        output_directory=embs_dir,
        output_prefix="emb",
    )
    embs = embs.sort_values("joinid")
    adata.obsm["geneformer"] = embs.drop(columns="joinid").to_numpy()

    return adata