import numpy as np

from recurscan.data.dataset import encode_clinical, make_synthetic_survival


def test_encode_clinical_shapes_and_imputation():
    recs = [
        {"AGE_AT_DIAGNOSIS": "55.0", "GRADE": "2", "ER_STATUS": "Positive",
         "CLAUDIN_SUBTYPE": "LumA", "CHEMOTHERAPY": "YES"},
        {"AGE_AT_DIAGNOSIS": "", "GRADE": "3", "ER_STATUS": "Negative",
         "CLAUDIN_SUBTYPE": "Basal"},
    ]
    X, names = encode_clinical(recs)
    assert X.shape[0] == 2
    assert not np.isnan(X).any()  # median imputed
    j_age = names.index("AGE_AT_DIAGNOSIS")
    assert X[0, j_age] == 55.0
    assert X[1, j_age] == 55.0  # imputed to median of present values
    assert X[0, names.index("ER_STATUS=Positive")] == 1.0
    assert X[1, names.index("CLAUDIN_SUBTYPE=Basal")] == 1.0


def test_synthetic_survival_effect_direction():
    X, t, e, beta = make_synthetic_survival(n=500, seed=9, effect=2.0,
                                            censor_rate=0.2)
    hi = X[:, 0] > 0
    assert t[hi].mean() < t[~hi].mean()  # higher risk -> shorter survival
    assert 0.5 < e.mean() < 1.0
