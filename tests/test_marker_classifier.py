from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCORER_PATH = (
    PROJECT_ROOT / "workflow" / "scripts" / "cls" / "marker_classifier.py"
)


@pytest.fixture(scope="module")
def scorer():
    spec = importlib.util.spec_from_file_location(
        "sencat_marker_classifier",
        SCORER_PATH,
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {SCORER_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_adata(matrix, *, labeled: bool = False) -> ad.AnnData:
    obs = pd.DataFrame(index=["sample-b", "sample-a"])
    if labeled:
        obs["is_sen"] = [1, 0]
    return ad.AnnData(
        X=matrix,
        obs=obs,
        var=pd.DataFrame(index=["gene_a", "gene_b"]),
    )


def markers() -> pd.DataFrame:
    return pd.DataFrame(
        {"coef": [0.25, -0.5]},
        index=["gene_b", "gene_a"],
    )


def expected_scores(matrix) -> np.ndarray:
    dense = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)
    return dense[:, 0] * -0.5 + dense[:, 1] * 0.25


@pytest.mark.parametrize(
    "matrix",
    [
        np.array([[1.0, 4.0], [3.0, 2.0]]),
        sparse.csr_matrix([[1.0, 4.0], [3.0, 2.0]]),
    ],
)
def test_unlabeled_dense_and_sparse_input(scorer, matrix):
    result = scorer.classify_samples(make_adata(matrix), markers())

    assert list(result.columns) == ["score"]
    assert list(result.index) == ["sample-b", "sample-a"]
    np.testing.assert_allclose(result["score"], expected_scores(matrix))


def test_labeled_input_remains_compatible(scorer):
    matrix = np.array([[1.0, 4.0], [3.0, 2.0]])
    result = scorer.classify_samples(
        make_adata(matrix, labeled=True),
        markers(),
    )

    assert list(result.columns) == ["label", "score"]
    assert result["label"].tolist() == [1, 0]
    np.testing.assert_allclose(result["score"], expected_scores(matrix))


def test_zero_marker_overlap_is_a_clear_error(scorer):
    missing = pd.DataFrame({"coef": [1.0]}, index=["missing_gene"])

    with pytest.raises(ValueError, match="No markers"):
        scorer.classify_samples(make_adata(np.ones((2, 2))), missing)


def test_cli_preserves_order_and_formula(tmp_path: Path):
    matrix = np.array([[0.0, 3.0], [1.0, 7.0]])
    adata_path = tmp_path / "input.h5ad"
    markers_path = tmp_path / "markers.csv"
    output_path = tmp_path / "scores.csv"
    log_path = tmp_path / "scorer.log"
    make_adata(matrix).write_h5ad(adata_path)
    markers().to_csv(markers_path)

    completed = subprocess.run(
        [
            sys.executable,
            str(SCORER_PATH),
            "--input-h5ad",
            str(adata_path),
            "--markers",
            str(markers_path),
            "--output-results-csv",
            str(output_path),
            "--log",
            str(log_path),
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr
    result = pd.read_csv(output_path, index_col=0)
    assert list(result.index) == ["sample-b", "sample-a"]
    np.testing.assert_allclose(
        result["score"],
        expected_scores(np.log1p(matrix)),
    )


def test_cli_zero_overlap_is_nonzero_without_output(tmp_path: Path):
    adata_path = tmp_path / "input.h5ad"
    markers_path = tmp_path / "markers.csv"
    output_path = tmp_path / "scores.csv"
    log_path = tmp_path / "scorer.log"
    make_adata(np.ones((2, 2))).write_h5ad(adata_path)
    pd.DataFrame({"coef": [1.0]}, index=["missing_gene"]).to_csv(markers_path)

    completed = subprocess.run(
        [
            sys.executable,
            str(SCORER_PATH),
            "--input-h5ad",
            str(adata_path),
            "--markers",
            str(markers_path),
            "--output-results-csv",
            str(output_path),
            "--log",
            str(log_path),
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode != 0
    assert not output_path.exists()
    assert "No markers" in completed.stderr + log_path.read_text()
