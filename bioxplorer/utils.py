
import os
import tarfile
import requests
from tqdm import tqdm


def update_code(file_path, old_line, new_line):
    """Updates a specific line in a file with a new line."""
    with open(file_path, 'r') as file:
        lines = file.readlines()
    lines = [line.replace(old_line, new_line) if line.strip() == old_line else line for line in lines]
    with open(file_path, 'w') as file:
        file.writelines(lines)


def figshare_download(url, save_path):
    """
    Figshare download helper with progress bar
    
    This function is a copy of UCE utils.py figshare_download function.
    Reference - https://github.com/snap-stanford/UCE/blob/main/utils.py#L72

    Args:
        url (str): the url of the dataset
        path (str): the path to save the dataset
    """

    if os.path.exists(save_path):
        return
    else:
        # Check if directory exists
        if not os.path.exists(os.path.dirname(save_path)):
            os.makedirs(os.path.dirname(save_path))
        print("Downloading " + save_path + " from " + url + " ..." + "\n")
        response = requests.get(url, stream=True)
        total_size_in_bytes = int(response.headers.get('content-length', 0))
        block_size = 1024
        progress_bar = tqdm(
            total=total_size_in_bytes,
            unit='iB',
            unit_scale=True
        )
        with open(save_path, 'wb') as file:
            for data in response.iter_content(block_size):
                progress_bar.update(len(data))
                file.write(data)
        progress_bar.close()

    # If the downloaded filename ends in tar.gz then extraact it
    if save_path.endswith(".tar.gz"):
       with tarfile.open(save_path) as tar:
            tar.extractall(path=os.path.dirname(save_path))
            print("Done!")