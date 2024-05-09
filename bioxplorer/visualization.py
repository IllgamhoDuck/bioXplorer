
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd


def draw_line_plot(
    data,
    methods=['uce', 'geneformer', 'scgpt'],
    scrna=True,
    pseudo_bulk=True,
    count=True,
):
    """
    Draw a line plot comparing F1 scores across multiple datasets for each method.

    data: dict
    {
        'dataset_id': {
            'scrna': {
                'uce': {
                    'f1': 0.5,
                    'total_data': 1000
                }
                'geneformer': {
                    'f1': 0.5,
                    'total_data': 1000
                }
            },
            'pseudo_bulk': {
                'uce': {
                    'f1': 0.5,
                    'total_data': 1000
                }
            }
        }
        '...': {...},
        '...': {...},
    }

    Args:
        data (dict): Dictionary containing the data for each dataset
        methods (list): List of methods to compare
        scrna (bool): Whether to include scRNA data in the plot
        pseudo_bulk (bool): Whether to include pseudo-bulk data in the plot
        count (bool): Whether to include total data count as annotation
    """
    assert scrna or pseudo_bulk, 'at least one of scrna or pseudo_bulk must be True'

    # Prepare figure
    fig, ax = plt.subplots(figsize=(12, 8))

    # Process each method
    for method_index, method in enumerate(methods):
        # Collect scores and total data across datasets
        scrna_scores = []
        scrna_total_data = []
        pseudobulk_scores = []
        pseudobulk_total_data = []

        for dataset_id, dataset in data.items():
            if scrna:
                if 'scrna' in dataset and method in dataset['scrna']:
                    scrna_scores.append(dataset['scrna'][method]['f1'])
                    scrna_total_data.append(dataset['scrna'][method]['total_data'])
                else:
                    scrna_scores.append(None)
                    scrna_total_data.append(None)
            if pseudo_bulk:
                if 'pseudo_bulk' in dataset and method in dataset['pseudo_bulk']:
                    pseudobulk_scores.append(dataset['pseudo_bulk'][method]['f1'])
                    pseudobulk_total_data.append(dataset['pseudo_bulk'][method]['total_data'])
                else:
                    pseudobulk_scores.append(None)
                    pseudobulk_total_data.append(None)

        # X positions for each dataset point
        x_positions = np.arange(len(data))

        # Plot scRNA scores for the method
        if scrna:
            scRNA_line, = ax.plot(x_positions, scrna_scores, '-o', label=f'scRNA {method}', alpha=0.7)

        # Plot pseudo-bulk scores for the method
        if pseudo_bulk:
            if any(score is not None for score in pseudobulk_scores):
                pseudoBulk_line, = ax.plot(x_positions, pseudobulk_scores, '-o', label=f'Pseudo-bulk {method}', alpha=0.7)

        # Annotate total data for each point
        if count:
            if scrna:
                for i, (score, txt) in enumerate(zip(scrna_scores, scrna_total_data)):
                    if score is not None and txt is not None:
                        ax.annotate(txt, (x_positions[i], score), textcoords="offset points", xytext=(0,10), ha='center')
            if pseudo_bulk:
                for i, (score, txt) in enumerate(zip(pseudobulk_scores, pseudobulk_total_data)):
                    if score is not None and txt is not None:
                        ax.annotate(txt, (x_positions[i], score), textcoords="offset points", xytext=(0,10), ha='center')

    # Customize plot
    ax.set_xlabel('Dataset Index')
    ax.set_ylabel('F1 Scores')
    ax.set_title('Comparison of F1 Scores Across Multiple Datasets by Method')
    ax.set_xticks(x_positions)
    ax.set_xticklabels([f'Dataset {i+1}' for i in x_positions])
    ax.legend()

    plt.show()


def draw_violin_plot(data, methods=['uce', 'geneformer', 'scgpt'], scrna=True, pseudo_bulk=True):
    """
    Draw a violin plot comparing F1 scores across multiple datasets for each method.
    """
    assert scrna or pseudo_bulk, 'At least one of scRNA or pseudo bulk must be true.'

    # Preparing data for the plot
    scores_list = []
    for dataset_id, types in data.items():
        for data_type in types:
            if (data_type == 'scrna' and scrna) or (data_type == 'pseudo_bulk' and pseudo_bulk):
                for method in methods:
                    if method in types[data_type]:
                        score = types[data_type][method]['f1']
                        scores_list.append({'Dataset ID': dataset_id, 'Method': method, 'Data Type': data_type, 'F1 Score': score})

    # Converting to DataFrame
    df = pd.DataFrame(scores_list)

    # Drawing the violin plot
    plt.figure(figsize=(12, 8))
    sns.violinplot(x='Method', y='F1 Score', hue='Data Type', data=df, split=True, inner='quart', scale='count')
    plt.title('F1 Score Distribution by Method and Data Type')
    plt.xlabel('Method')
    plt.ylabel('F1 Score')
    plt.legend(title='Data Type')
    plt.grid(True)
    plt.show()