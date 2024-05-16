

import os
import json
import subprocess
import warnings
import anndata
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
    return {}

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


def extract_uce(
    adata,
    ensembl_name_col='feature_name',
    batch_size=24,
    layers=33,
    filter=False,
    cache=None,
    verbose=False,
    *args,
    **kwargs
):
    adata = adata.copy()

    current_dir = os.getcwd()

    # change folder to UCE
    package_path = os.path.join(TEMP_DIR, 'uce', 'UCE')
    os.chdir(package_path)

    # preprocess for UCE
    adata.var.reset_index(inplace=True)
    adata.var.set_index(ensembl_name_col, inplace=True)
    adata.var.rename(columns={'index': 'ensembl_id'}, inplace=True)

    # save h5ad
    data_dir = os.path.join(TEMP_DIR, "uce", "data")
    if os.path.exists(data_dir):
        subprocess.run(["rm", "-r", data_dir])
    os.makedirs(data_dir)

    data_path = os.path.join(data_dir, f'uce.h5ad')
    adata.write_h5ad(data_path)

    # setup UCE args
    assert layers in [4, 33], "Only support 4 or 33 layers"
    if layers == 4:
        model_dir = os.path.join(TEMP_DIR, "uce", "model", "model4l")
    elif layers == 33:
        model_dir = os.path.join(TEMP_DIR, "uce", "model", "model33l")

    # process with UCE
    cmds = [
        'python', os.path.join(package_path, 'eval_single_anndata.py'),
        '--adata_path', data_path,
        '--dir', data_dir + '/', # becaues uce adds two path with just `+`
        '--model_loc', model_dir,
        '--species', 'human',
        '--batch_size', str(batch_size),
        '--nlayers', str(layers),
    ]
    if filter:
        cmds.append('--filter')
    subprocess.run(cmds)

    save_path = os.path.join(data_dir, 'uce_uce_adata.h5ad')
    adata = anndata.read_h5ad(save_path)

    os.chdir(current_dir)

    return adata


# check if scGPT is already installed
IS_SCGPT_IMPORTED = False

def extract_scgpt(
    adata,
    ensembl_name_col="feature_name",
    batch_size=32,
    cache=None,
    verbose=False,
    *args,
    **kwargs
):
    global IS_SCGPT_IMPORTED
    if not IS_SCGPT_IMPORTED:
        import scgpt
        IS_SCGPT_IMPORTED = True

    adata = adata.copy()
    adata = scgpt.tasks.embed_data(
        adata,
        cache['model_dir'],
        gene_col=ensembl_name_col,
        batch_size=batch_size,
    )

    return adata

def extract_geneformer(
    adata,
    ensembl_id_col='feature_id',
    batch_size=None,
    cache=None,
    verbose=False,
    *args,
    **kwargs,
):
    adata = adata.copy()

    # if batch_size is provided, use it
    if batch_size is not None:
        cache['emb_extractor'].forward_batch_size = batch_size

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