
import os
import subprocess
import sys
import warnings
from tqdm.notebook import tqdm
from bioxplorer.utils import update_code
from bioxplorer.utils import figshare_download

warnings.filterwarnings("ignore")


TEMP_DIR = os.path.join(os.path.expanduser('~'), '.bioxplorer')

def install_scgpt(
    package_dir: str = os.path.join(TEMP_DIR, 'scgpt'),
):
    """post-install for scgpt"""

    # check if scGPT is already installed
    try:
        import scgpt
        print("scGPT is already installed.")
    except ImportError:
        subprocess.check_call([
            sys.executable, '-m', 'pip', 'install', '-U',
            'scgpt', 'wandb', 'louvain', 'torch', 'torchtext', 'gdown'
        ])

    data_dir: str = os.path.join(package_dir, 'data')
    model_dir: str = os.path.join(package_dir, 'model')

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

    print("scGPT Environment is set up successfully.")


def install_geneformer(
    package_dir: str = os.path.join(TEMP_DIR, 'geneformer'),
):
    """Post-installation for geneformer"""

    # Check if geneformer is already installed
    try:
        import geneformer
        print("Geneformer is already installed.")
    except ImportError:
        try:
            subprocess.check_call(['git', 'lfs', 'install'])
        except subprocess.CalledProcessError:
            lfs_msg = (
                "Git LFS is not installed or recognized. Please install Git LFS manually.\n"
                "For Debian/Ubuntu-based systems, you can install Git LFS using the following commands:\n"
                "  sudo apt-get update\n"
                "  sudo apt-get install git-lfs\n"
                "After installation, run 'git lfs install' to initialize Git LFS or run this code again\n"
            )

            # package install error
            raise ImportError(lfs_msg)

        # Clone the geneformer repository
        os.makedirs(package_dir, exist_ok=True)
        package_path = os.path.join(package_dir, 'Geneformer')
        if not os.path.exists(package_path):
            os.makedirs(package_path, exist_ok=True)
            subprocess.check_call([
                'git', 'clone', 'https://huggingface.co/ctheodoris/Geneformer',
                package_path,
            ])

        # Update code for tokenizer
        file_path = os.path.join(package_path, 'geneformer', 'tokenizer.py')
        old_line = 'tokenized_dataset.save_to_disk(output_path)'
        new_line = 'tokenized_dataset.save_to_disk(str(output_path))'
        update_code(file_path, old_line, new_line)

        # Install the local geneformer package
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', package_path])

    # Sync model from AWS
    model_dir: str = os.path.join(package_dir, 'models')
    os.makedirs(model_dir, exist_ok=True)
    subprocess.check_call([
        'aws', 's3', 'sync', '--no-sign-request', '--no-progress', '--only-show-errors',
        's3://cellxgene-contrib-public/models/geneformer/2023-12-15/homo_sapiens/fined-tuned-model/',
        os.path.join(model_dir, 'fine_tuned_geneformer'),
    ])

    print("Geneformer Environment is set up successfully")


def install_uce(
    package_dir: str = os.path.join(TEMP_DIR, 'uce'),
):
    """Post-installation for UCE"""

    # Setup directory for UCE
    os.makedirs(package_dir, exist_ok=True)
    package_path = os.path.join(package_dir, 'UCE')
    if not os.path.exists(package_path):
        os.makedirs(package_path, exist_ok=True)
        # Assuming UCE is to be cloned from a repository
        subprocess.check_call([
            'git', 'clone', 'https://github.com/snap-stanford/UCE.git',
            package_path,
        ])

    model4l_url = "https://figshare.com/ndownloader/files/42706576"
    model33l_url = "https://figshare.com/ndownloader/files/43423236"

    model_dir = os.path.join(package_dir, 'model')
    model4l_path = os.path.join(model_dir, "model4l")
    model33l_path = os.path.join(model_dir, "model33l")
    if not os.path.exists(model4l_path):
        figshare_download(model4l_url, model4l_path)
    if not os.path.exists(model33l_path):
        figshare_download(model33l_url, model33l_path)

    print("UCE setup is completed successfully.")

