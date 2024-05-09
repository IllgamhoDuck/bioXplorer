
import os
import subprocess
import sys


TEMP_DIR = os.path.join(os.path.expanduser('~'), '.bioxplorer')

def install_scgpt(
    data_dir: str = os.path.join(TEMP_DIR, 'scGPT', 'data'),
    model_dir: str = os.path.join(TEMP_DIR, 'scGPT', 'model'),
):
    """post-install for scgpt"""

    # check if scGPT is already installed
    try:
        import scgpt
        print("scGPT is already installed.")
        return
    except ImportError:
        pass

    subprocess.check_call([
        sys.executable, '-m', 'pip', 'install', '-U',
        'scgpt', 'wandb', 'louvain', 'torchtext', 'gdown'
    ])

    os.makedirs(data_dir, exist_ok=True)
    data_file = os.path.join(data_dir, "human_pancreas_norm_complexBatch.h5ad")
    if not os.path.exists(data_file):
        subprocess.check_call([
            'wget', '--content-disposition',
            'https://figshare.com/ndownloader/files/24539828',
            '-O', data_file
        ])

    if not os.path.exists(model_dir):
        os.makedirs(model_dir, exist_ok=True)
        import gdown
        gdown.download_folder(
            "https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y",
            output=model_dir,
        )

    print("scGPT has been installed successfully.")
