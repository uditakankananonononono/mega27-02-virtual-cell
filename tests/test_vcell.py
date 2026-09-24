"""Hermetic test suite: runs entirely on cached in-repo data, no network."""
import json
import os

import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")

import sys
sys.path.insert(0, ROOT)


@pytest.fixture(scope="module")
def core_model():
    import warnings
    warnings.filterwarnings("ignore")
    from vcell import metabolism as vm
    return vm.load_model(os.path.join(DATA, "e_coli_core.json"))


# ---------------- metabolism ----------------
class TestMetabolism:
    def test_model_loads_real(self, core_model):
        assert len(core_model.reactions) == 95
        assert len(core_model.genes) == 137

    def test_glucose_growth_matches_published_core_value(self, core_model):
        from vcell import metabolism as vm
        cond = vm.Condition("EX_glc__D_e", aerobic=True)
        gr = vm.growth_rate(core_model, cond, vm.carbon_exchange_ids(core_model))
        assert gr == pytest.approx(0.873922, abs=1e-4)

    def test_medium_restored_after_condition(self, core_model):
        from vcell import metabolism as vm
        cids = vm.carbon_exchange_ids(core_model)
        vm.growth_rate(core_model, vm.Condition("EX_succ_e", aerobic=False), cids)
        gr = vm.growth_rate(core_model, vm.Condition("EX_glc__D_e", aerobic=True), cids)
        assert gr == pytest.approx(0.873922, abs=1e-4)

    def test_anaerobic_glc_lower_than_aerobic(self, core_model):
        from vcell import metabolism as vm
        cids = vm.carbon_exchange_ids(core_model)
        aer = vm.growth_rate(core_model, vm.Condition("EX_glc__D_e", True), cids)
        ana = vm.growth_rate(core_model, vm.Condition("EX_glc__D_e", False), cids)
        assert 0 < ana < aer

    def test_essentiality_scan_regression(self, core_model):
        from vcell import metabolism as vm
        ess = vm.gene_essentiality_scan(core_model)
        assert int(ess.predicted_essential.sum()) == 7
        assert set(ess[ess.predicted_essential].gene) == {"b0720", "b1136", "b1779", "b2415", "b2416", "b2779", "b2926"}
        assert float(ess.wt_growth_h.iloc[0]) == pytest.approx(0.873922, abs=1e-4)


# ---------------- data ----------------
class TestData:
    def test_gerdes_parse(self):
        from vcell import data as vd
        g = vd.parse_gerdes_s1(os.path.join(DATA, "raw", "gerdes_table_s1.txt"))
        assert len(g) == 4226
        lab = vd.essentiality_labels(g)
        assert int(lab.essential.sum()) == 617
        assert len(lab) == 3689

    def test_labels_overlap_model(self, core_model):
        from vcell import data as vd
        lab = vd.essentiality_labels(vd.parse_gerdes_s1(os.path.join(DATA, "raw", "gerdes_table_s1.txt")))
        model_genes = {g.id for g in core_model.genes}
        assert len(set(lab.bnumber) & model_genes) == 125


# ---------------- dynamics ----------------
class TestDynamics:
    def test_unregulated_coutilization(self, core_model):
        from vcell import dynamics as vdyn
        traj = vdyn.dynamic_fba(core_model, {"EX_glc__D_e": 10.0, "EX_ac_e": 2.0},
                                dt=0.2, t_end=8.0)
        met = vdyn.diauxie_metrics(traj, "EX_glc__D_e", "EX_ac_e")
        assert met["co_utilization_primary_phase"] is True

    def test_regulated_diauxie(self, core_model):
        from vcell import dynamics as vdyn
        traj = vdyn.dynamic_fba(core_model, {"EX_glc__D_e": 10.0, "EX_ac_e": 2.0},
                                dt=0.2, t_end=8.0,
                                regulators=[("EX_glc__D_e", "EX_ac_e", 0.0)])
        met = vdyn.diauxie_metrics(traj, "EX_glc__D_e", "EX_ac_e")
        assert met["co_utilization_primary_phase"] is False
        assert met["t_primary_exhausted_h"] == pytest.approx(5.0, abs=1.0)
        assert np.all(np.diff(traj["biomass"]) >= -1e-9)  # biomass never declines


# ---------------- sequence CNN pieces ----------------
class TestSeqCNN:
    def test_one_hot(self):
        from vcell import seqcnn as vs
        x = vs.one_hot("ACGTN", 8)
        assert x.shape == (4, 8)
        assert x[:, 0].tolist() == [1, 0, 0, 0]
        assert x[:, 4].sum() == 0.0  # N -> all zeros

    def test_kmer_normalized(self):
        from vcell import seqcnn as vs
        v = vs.kmer_freqs("ACGTACGTACGT", 3)
        assert v.shape == (64,)
        assert v.sum() == pytest.approx(1.0)

    def test_cds_extraction(self):
        from vcell import seqcnn as vs
        cds = vs.extract_cds(os.path.join(DATA, "U00096.3.gb"))
        assert len(cds) == 4218
        assert {"bnumber", "gene_name", "cds_len", "cds"} <= set(cds.columns)

    def test_cnn_forward_shape(self):
        import torch
        from vcell import seqcnn as vs
        m = vs.make_cnn(1200)
        out = m(torch.zeros(2, 4, 1200))
        assert out.shape == (2, 1)


# ---------------- graph GNN pieces ----------------
class TestGraphGNN:
    def test_gcn_norm_properties(self):
        from vcell import graphgnn as vg
        A = np.array([[0, 1], [1, 0]], dtype=np.float32)
        An = vg.gcn_norm(A)
        assert An.shape == (2, 2)
        assert np.allclose(An, An.T)  # symmetric

    def test_graph_symmetric_binary(self, core_model):
        from vcell import graphgnn as vg
        genes = sorted(g.id for g in core_model.genes)[:40]
        A = vg.build_gene_graph(core_model, genes)
        assert np.allclose(A, A.T)
        assert set(np.unique(A)) <= {0.0, 1.0}

    def test_gnn_toy_training(self):
        from vcell import graphgnn as vg
        rng = np.random.default_rng(0)
        A = (rng.random((30, 30)) > 0.8).astype(np.float32)
        A = np.maximum(A, A.T)
        X = rng.normal(size=(30, 4)).astype(np.float32)
        y = (X[:, 0] > 0).astype(int)
        res = vg.train_eval_gnn(A, X, y, folds=2, epochs=5)
        assert 0.0 <= res["oof_auroc"] <= 1.0
        assert len(res["oof_scores"]) == 30


# ---------------- benchmark ----------------
class TestBenchmark:
    def _df(self):
        return pd.DataFrame({"essential": [0, 0, 1, 1, 0, 1, 0, 1],
                             "score": [0.1, 0.2, 0.9, 0.8, 0.3, 0.7, 0.2, 0.6],
                             "pred": [0, 0, 1, 1, 0, 1, 1, 0]})

    def test_metrics(self):
        from vcell import benchmark as vb
        rep = vb.classification_report(self._df(), "score", pred_col="pred")
        assert rep["auroc"] == pytest.approx(1.0)
        assert 0 <= rep["f1"] <= 1

    def test_bootstrap_ci_ordered(self):
        from vcell import benchmark as vb
        ci = vb.bootstrap_ci(self._df(), "score", n_boot=200)
        assert ci["auroc_ci"][0] <= ci["auroc_ci"][1]


def test_cli_audit_biomass_core(tmp_path):
    """CLI end-to-end on the bundled core model for a known pair of genes."""
    import pandas as pd
    from pathlib import Path
    from vcell.cli import main
    model = str(Path(__file__).resolve().parents[1] / "data" / "e_coli_core.json")
    out = tmp_path / "audit.csv"
    rc = main(["audit-biomass", model, "--genes", "b2415", "b1136", "--out", str(out)])
    assert rc == 0
    df = pd.read_csv(out)
    assert set(df.gene) == {"b2415", "b1136"}
    assert set(df.verdict) <= {"biomass_forced", "network_forced", "network_forced_no_biomass_link", "not_essential"}


def test_bernstein_prauc_perfect_and_bootstrap():
    import numpy as np
    from vcell.bernstein import bernstein_prauc, paired_gene_bootstrap, greedy_select
    fit = np.array([[-3.0, 0.1], [-2.5, 0.0], [0.2, -0.1], [0.1, 0.3]])
    sim = np.array([[0.0, 1.0], [0.0, 1.0], [1.0, 1.0], [1.0, 1.0]])  # no-growth exactly where fitness lowest
    assert abs(bernstein_prauc(sim, fit) - 1.0) < 1e-9
    # a supplement that wrongly rescues gene 0 on carbon 0 must be rejected by greedy selection
    pairs = np.array([[0, 0], [1, 0]]); R = np.array([[1, 0], [0, 0]], dtype=np.int8)
    S, trace, apply = greedy_select(sim, fit, pairs, R, ["bad", "inert"], [0, 1])
    assert S == []
    bs = paired_gene_bootstrap((sim > 0.001).astype(int), (sim > 0.001).astype(int), fit, [0, 1], n=50)
    assert bs["diff"] == 0.0


def test_rescue_audit_on_off_pathway(tmp_path):
    """Hermetic: icd (b1136) knockout rescued by 2-oxoglutarate is on-pathway (TCA);
    enolase (b2779) rescued by pyruvate is off-pathway when pyruvate is tagged TCA only."""
    import warnings
    warnings.filterwarnings("ignore")
    from pathlib import Path
    from vcell import metabolism as vm
    from vcell.rescue import rescue_audit, load_kegg_links
    from vcell.cli import main
    mp = Path(__file__).resolve().parents[1] / "data" / "e_coli_core.json"
    m = vm.load_model(str(mp))
    kg = tmp_path / "links.tsv"
    kg.write_text("eco:b1136\tpath:eco00020\neco:b2779\tpath:eco00010\n")
    gp = load_kegg_links(str(kg))
    assert gp == {"b1136": {"00020"}, "b2779": {"00010"}}
    rows = rescue_audit(m, {"akg_c": {"00020"}, "pyr_c": {"00020"}}, gp, genes=["b1136", "b2779", "b0008"])
    lab = {(r["gene"], r["supplement"]): r["label"] for r in rows}
    assert lab[("b1136", "akg_c")] == "on_pathway"
    assert lab[("b2779", "pyr_c")] == "off_pathway"
    assert not any(r["gene"] == "b0008" for r in rows)  # non-essential gene is never a rescue
    out = tmp_path / "r.csv"
    assert main(["rescue-audit", str(mp), "--supplement", "akg_c=00020", "--gene-pathways", str(kg),
                 "--genes", "b1136", "--out", str(out)]) == 0
    assert "on_pathway" in out.read_text()


def test_rescue_audit_production_label():
    """Hermetic: icd (b1136) knockout blocks 2-oxoglutarate synthesis, so its akg rescue is on-pathway."""
    import warnings
    warnings.filterwarnings("ignore")
    from pathlib import Path
    from vcell import metabolism as vm
    from vcell.rescue import rescue_audit_production, can_produce
    m = vm.load_model(str(Path(__file__).resolve().parents[1] / "data" / "e_coli_core.json"))
    assert can_produce(m, "akg_c")
    rows = rescue_audit_production(m, ["akg_c"], genes=["b1136", "b0008"])
    assert [(r["gene"], r["label"]) for r in rows] == [("b1136", "on_pathway")]
