"""G2 hard-call gate under amendment A1: train-fold-only thresholds, candidate ladder R1-R5.
Comparator: FBA minimal with the SAME train-only MCC threshold rule (symmetric);
fba_min > 0.5 reported as sensitivity. Pooled exact McNemar; G2 passes iff p>=0.05 and b>=c."""
import json, hashlib, numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score, f1_score, matthews_corrcoef
from scipy.stats import binomtest

al = pd.read_csv('results/aligned_predictions_all_models.csv')
oof = pd.read_csv('results/ensemble_v2_oof.csv')
ext = pd.read_csv('results/external_features.csv')
df = al.merge(oof.drop(columns=['essential']), on='bnumber')
y = df.essential.values.astype(int)
folds = np.zeros(len(y), int)
for k, (_, te) in enumerate(StratifiedKFold(3, shuffle=True, random_state=7).split(np.zeros(len(y)), y)):
    folds[te] = k

def pick_t(s_tr, y_tr):
    grid = np.quantile(s_tr, np.linspace(0.005, 0.995, 199))
    best_m, best_t = -2, 0.5
    for t in grid:
        m = matthews_corrcoef(y_tr, s_tr >= t)
        if m > best_m + 1e-12 or (abs(m - best_m) <= 1e-12 and abs(t - 0.5) < abs(best_t - 0.5)):
            best_m, best_t = m, t
    return best_t

def hard_calls(score, transform=None):
    """Per fold: optional calibration fit on train folds, threshold from train folds."""
    pred = np.zeros(len(y), int); ts = {}
    for k in range(3):
        tr, te = folds != k, folds == k
        s_tr, s_te = score[tr], score[te]
        if transform == 'isotonic':
            f = IsotonicRegression(out_of_bounds='clip').fit(s_tr, y[tr])
            s_tr, s_te = f.predict(s_tr), f.predict(s_te)
        elif transform == 'platt':
            f = LogisticRegression(max_iter=4000).fit(s_tr.reshape(-1, 1), y[tr])
            s_tr = f.predict_proba(s_tr.reshape(-1, 1))[:, 1]
            s_te = f.predict_proba(s_te.reshape(-1, 1))[:, 1]
        t = pick_t(s_tr, y[tr]); ts[f'fold{k}'] = float(t)
        pred[te] = (s_te >= t).astype(int)
    return pred, ts

def oversampled_lr_oof():
    base = ['fba_min', 'fba_rich', 'gnn', 'cnn', 'kmer']
    NEW3 = ['pax_log_ppm', 'cai', 'gc3']
    extc = [c for c in ext.columns if c != 'bnumber' and not c.startswith('strctx') and c not in NEW3]
    d2 = al.merge(ext, on='bnumber', how='left').fillna(0)
    X = d2[base + extc].values
    out = np.zeros(len(y))
    for k, (tr, te) in enumerate(StratifiedKFold(3, shuffle=True, random_state=7).split(X, y)):
        rng = np.random.default_rng(7)
        Xtr, ytr = X[tr], y[tr]
        mi, ma = (ytr == 1), (ytr == 0)
        extra = rng.choice(np.where(mi)[0], size=ma.sum() - mi.sum(), replace=True)
        Xb = np.vstack([Xtr, Xtr[extra]]); yb = np.concatenate([ytr, ytr[extra]])
        m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=4000, C=0.5))
        out[te] = m.fit(Xb, yb).predict_proba(X[te])[:, 1]
    return out

def mcnemar(pred_ens, pred_fba):
    eok, fok = pred_ens == y, pred_fba == y
    b, c = int(((~fok) & eok).sum()), int((fok & (~eok)).sum())
    return {'b_ensemble_wins': b, 'c_fba_wins': c, 'p_exact_2sided': binomtest(b, b + c).pvalue,
            'net': b - c, 'g2_pass': bool(binomtest(b, b + c).pvalue >= 0.05 and b >= c)}

# comparator: FBA minimal, symmetric train-only MCC threshold
fba_pred, fba_ts = hard_calls(df.fba_min.values)
fba_fixed = (df.fba_min > 0.5).astype(int).values
res = {'comparator': {'rule': 'fba_min, per-fold train-only MCC threshold (symmetric)',
                      'thresholds': fba_ts,
                      'metrics': {'f1': f1_score(y, fba_pred), 'mcc': matthews_corrcoef(y, fba_pred),
                                  'acc': float((fba_pred == y).mean())},
                      'sensitivity_fba_gt_0.5': {'f1': f1_score(y, fba_fixed),
                                                 'acc': float((fba_fixed == y).mean())}},
       'baseline_reproduction_check': {'rule': 'v1 ensemble best-F1 (99-pt grid, eval-inclusive) vs fba_min best-F1',
                                       'note': 'amendment A1 harness check'}}
# harness check: reproduce locked current-state 27:50
from sklearn.metrics import f1_score as _f1
def best_f1(s):
    g = np.quantile(s, np.linspace(0.01, 0.99, 99)); bm, bt = 0, .5
    for t in g:
        f = _f1(y, s >= t)
        if f > bm: bm, bt = f, t
    return bm, bt
v1 = pd.read_csv('results/ensemble_oof.csv')
d1 = al.merge(v1[['bnumber', 'ensemble']], on='bnumber')
fe, te = best_f1(d1.ensemble.values); ff, tf = best_f1(d1.fba_min.values)
eok = (d1.ensemble.values >= te) == y; fok = (d1.fba_min.values >= tf) == y
res['baseline_reproduction_check'].update(
    {'b': int(((~fok) & eok).sum()), 'c': int((fok & (~eok)).sum()),
     'ens_f1': fe, 'fba_f1': ff, 'matches_locked_27_50': bool(int(((~fok) & eok).sum()) == 27 and int((fok & (~eok)).sum()) == 50)})

rungs = [('R1_v2_lr_trainonly_mcc', df.v2_lr.values, None),
         ('R2_v2_lr_isotonic', df.v2_lr.values, 'isotonic'),
         ('R3_v2_lr_platt', df.v2_lr.values, 'platt'),
         ('R4_v2_lr_oversampled_refit', oversampled_lr_oof(), None),
         ('R5_v2_lr_gbm_prob_avg', (df.v2_lr.values + df.v2_gbm.values) / 2, None)]
res['ladder'] = []
first_pass = None
for name, s, tr in rungs:
    pred, ts = hard_calls(s, tr)
    mc = mcnemar(pred, fba_pred)
    mc_sens = mcnemar(pred, fba_fixed)
    row = {'rung': name, 'auroc': roc_auc_score(y, s), 'f1': f1_score(y, pred),
           'mcc': matthews_corrcoef(y, pred), 'thresholds': ts, 'mcnemar_vs_comparator': mc,
           'mcnemar_vs_fba_fixed_0.5_sensitivity': mc_sens,
           'g1_pass': bool(roc_auc_score(y, s) > 0.666)}
    res['ladder'].append(row)
    if first_pass is None and mc['g2_pass'] and row['g1_pass']:
        first_pass = name
res['verdict'] = {'first_passing_rung': first_pass,
                  'g2_overall': 'PASS' if first_pass else 'FAIL - pivot per base prereg (rule 6)'}
res['input_sha256'] = {p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
                       for p in ['results/aligned_predictions_all_models.csv', 'results/ensemble_v2_oof.csv',
                                 'results/ensemble_oof.csv', 'results/external_features.csv']}
json.dump(res, open('results/calibration_gate.json', 'w'), indent=1)
print(json.dumps({'baseline_check': res['baseline_reproduction_check'],
                  'verdict': res['verdict'],
                  'ladder_summary': [{'rung': r['rung'], 'auroc': round(r['auroc'], 4),
                                      'f1': round(r['f1'], 4), 'b': r['mcnemar_vs_comparator']['b_ensemble_wins'],
                                      'c': r['mcnemar_vs_comparator']['c_fba_wins'],
                                      'p': round(r['mcnemar_vs_comparator']['p_exact_2sided'], 4),
                                      'g2': r['mcnemar_vs_comparator']['g2_pass']} for r in res['ladder']]}, indent=1))
