'''
---------------------
sceptic functions
author: Gang Li
e-mail:gangliuw@uw.edu
MIT LICENSE
---------------------
'''
from sklearn.model_selection import KFold, GridSearchCV, StratifiedKFold
from sklearn import svm
import numpy as np
import sklearn
import xgboost as xgb

eFold=3
iFold=4
def run_sceptic_and_evaluate(data, labels, label_list, parameters, method="svm", use_gpu=False):
    """
    Run pseudotime estimation using SVM or XGBoost.
    
    Args:
        data (np.ndarray): Cell-by-feature matrix.
        labels (np.ndarray): Ground-truth time labels.
        label_list (List[float or int]): Ordered time labels (e.g., [0, 1, 2]).
        parameters (dict): Grid search parameters for the classifier.
        method (str): "svm" or "xgboost".
        use_gpu (bool): Only applies if method="xgboost".
        
    Returns:
        cm, label_predicted, pseudotime, sceptic_prob
    """
    # Set default parameters if none provided
    if not parameters:
        if method == "svm":
            parameters = {
                "C": [1, 10],
                "kernel": ["linear", "rbf"],
                "gamma": ["scale"]
            }
        elif method == "xgboost":
            parameters = {
                "max_depth": [3, 5],
                "learning_rate": [0.1, 0.3],
                "n_estimators": [100],
                "subsample": [0.8]
            }
        else:
            raise ValueError(f"Unsupported method '{method}' and no parameters provided.")
        
    data = np.asarray(data)
    labels = np.asarray(labels).astype(np.int64) # Ensure labels are correct type

    cm = np.zeros((len(label_list), len(label_list)))
    label_predicted = np.zeros(len(labels))
    sceptic_prob = np.zeros((len(labels), len(label_list)))
    pseudotime = np.zeros(len(labels))
    
    # Outer CV: Must be StratifiedKFold and requires (data, labels) for splitting
    kf = StratifiedKFold(n_splits=eFold, shuffle=True, random_state=42)
    
    for i, (train_index, test_index) in enumerate(kf.split(data, labels)):
        X_train, X_test = data[train_index], data[test_index]
        y_train, y_test = labels[train_index], labels[test_index]

        # Inner CV: Must be a StratifiedKFold object, not just the number iFold
        inner_cv = StratifiedKFold(n_splits=iFold, shuffle=True, random_state=42)

        if method == "xgboost":
            xgb_model = xgb.XGBClassifier(
                tree_method='auto',
                objective='multi:softprob',
                num_class=len(label_list), 
                eval_metric='mlogloss',
                use_label_encoder=False,
                random_state=42
            )
            clf = GridSearchCV(xgb_model, parameters, cv=inner_cv, scoring='roc_auc', n_jobs=-1)
            
        elif method == "svm":
            svc = svm.SVC(probability=True, random_state=42)
            clf = GridSearchCV(svc, parameters, cv=inner_cv, scoring='roc_auc', n_jobs=-1)
            
        else:
            raise ValueError(f"Unsupported method '{method}'. Choose 'svm' or 'xgboost'.")

        clf.fit(X_train, y_train)
        
        # Use predict_proba and argmax to ensure a 1D label array is returned.
        # This resolves the shape mismatch (334, 2) vs (334,) when assigning to label_predicted.
        prob = clf.predict_proba(X_test)
        predicted = np.argmax(prob, axis=1)
        
        label_predicted[test_index] = predicted
        cm += sklearn.metrics.confusion_matrix(y_test, predicted)
        
        # The 'prob' variable is already defined and correct, no need for the try/except block 
        # that was previously just assigning 'prob' again. We use the existing 'prob' variable.
        # try:
        #     prob = clf.predict_proba(X_test)
        # except Exception as e:
        #     print(f"Warning: predict_proba failed on fold {i}: {e}")
        #     prob = np.zeros((len(X_test), len(label_list)))

        sceptic_prob[test_index, :] = prob
        # Pseudotime calculation: weighted average of probabilities and time labels (0, 1)
        pseudotime[test_index] = np.sum(prob * label_list, axis=1)

    return cm, label_predicted, pseudotime, sceptic_prob
