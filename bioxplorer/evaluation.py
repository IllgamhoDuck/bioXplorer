
import numpy as np
from copy import deepcopy
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import confusion_matrix



def evaluate_knn_f1(
    adata,
    embedding_name,
    classify_col='cell_type'
):
    """
    Calculate the performance of the classifier using KNN

    Args:
        adata (AnnData object): Anndata object containing the data
        embedding_name (str): Name of the embedding in adata.obsm
        classify_col (str): Column name in adata.obs containing the labels
    """
    adata = deepcopy(adata)
    adata.obs.reset_index(inplace=True)

    # Extracting labels and groups
    X = adata.obsm[embedding_name]
    y = adata.obs[classify_col].to_numpy() # Cell type labels

    # Splitting the dataset into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, random_state=42)
    knn = KNeighborsClassifier(n_neighbors=5)
    knn.fit(X_train, y_train)
    y_pred = knn.predict(X_test)

    cm = confusion_matrix(y_test, y_pred)
    if cm.size == 1:  # If there is only one class present in both y_test and y_pred
        TN = FP = FN = 0
        TP = cm[0, 0] if len(np.unique(y_test)) == 1 else 0
    elif cm.size == 4:
        TN, FP, FN, TP = cm.ravel()
    else:
        raise ValueError("Unexpected confusion matrix size")

    return {'TN': TN, 'FP': FP, 'FN': FN, 'TP': TP}

def calculate_weighted_f1(tn, fp, fn, tp):
    # Calculate precision, recall, and F1 score for each class
    precision_class_1 = tp / (tp + fp) if (tp + fp) != 0 else 0
    recall_class_1 = tp / (tp + fn) if (tp + fn) != 0 else 0
    f1_class_1 = 2 * (precision_class_1 * recall_class_1) / (precision_class_1 + recall_class_1) if (precision_class_1 + recall_class_1) != 0 else 0

    precision_class_0 = tn / (tn + fn) if (tn + fn) != 0 else 0
    recall_class_0 = tn / (tn + fp) if (tn + fp) != 0 else 0
    f1_class_0 = 2 * (precision_class_0 * recall_class_0) / (precision_class_0 + recall_class_0) if (precision_class_0 + recall_class_0) != 0 else 0

    # Calculate support for each class
    support_class_1 = tp + fn
    support_class_0 = tn + fp

    # Calculate the total number of instances
    total_instances = support_class_1 + support_class_0

    # Calculate the weighted F1 score
    weighted_f1 = (f1_class_1 * support_class_1 + f1_class_0 * support_class_0) / total_instances

    return weighted_f1

def confusion_matrix_to_f1(confusion_matrix):
    """
    confusion_matrix = {
        'TN': 0,
        'FP': 0,
        'FN': 0,
        'TP': 0
    }
    """
    tn = confusion_matrix['TN']
    fp = confusion_matrix['FP']
    fn = confusion_matrix['FN']
    tp = confusion_matrix['TP']

    weighted_f1 = calculate_weighted_f1(tn, fp, fn, tp)
    total_data = tn + fp + fn + tp
    return weighted_f1, total_data