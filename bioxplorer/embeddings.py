
from bioxplorer.environment import TEMP_DIR
from bioxplorer.environment import install_uce
from bioxplorer.environment import install_geneformer
from bioxplorer.environment import install_scgpt



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

        self.extract_fns = {
            'uce': extract_uce,
            'scgpt': extract_scgpt,
            'geneformer': extract_geneformer
        }
        self.install_fns = {
            'uce': install_uce,
            'scgpt': install_scgpt,
            'geneformer': install_geneformer
        }
        self.init_fns = {fn: False for fn in self.install_fns}

    def __repr__(self):
        return "Embeddings()"

    def install(self, model):
        """
        install the model only when it is used
        """
        if model in self.install_fns:
            if self.init_fns[model]:
                return

            self.install_fns[model]()
            self.init_fns[model] = True
        else:
            raise ValueError(f"Model {model} not supported.")

    def extract(self, model, adata):
        self.install(model)
        if model in self.extract_fns:
            return self.extract_fns[model](adata)
        else:
            raise ValueError(f"Model {model} not supported.")


def extract_uce(adata):
    return adata

def extract_scgpt(adata):
    return adata

def extract_geneformer(adata):
    return adata