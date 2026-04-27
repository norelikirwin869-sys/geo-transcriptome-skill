from pathlib import Path


def test_enrichment_script_uses_base_subsetting_for_dynamic_columns():
    content = Path('scripts/03_enrichment_gsea.R').read_text(encoding='utf-8')
    assert '.data[[' not in content
    assert 'deg[[p_col]]' in content
    assert 'deg[[fc_col]]' in content


def test_ml_template_contains_documented_svm_rfe_step():
    content = Path('scripts/05_ml_biomarker_template.R').read_text(encoding='utf-8')
    assert 'svm_rfe <- rfe(' in content
    assert 'svm_rfe_genes.txt' in content
    assert 'Reduce(intersect, list(lasso_genes, top_rf, svm_genes))' in content


def test_readme_example_title_uses_mitochondrial_term():
    content = Path('README.md').read_text(encoding='utf-8')
    assert 'mitochondrial-related biomarkers' in content


def test_ml_template_has_small_feature_guard_and_dynamic_svm_sizes():
    content = Path('scripts/05_ml_biomarker_template.R').read_text(encoding='utf-8')
    assert 'Need at least 2 candidate genes present in expression matrix for ML modeling.' in content
    assert 'svm_sizes <- unique(sort(pmin(c(5, 10, 20, 30), n_features)))' in content
    assert 'sizes = svm_sizes' in content


def test_config_uses_canonical_mcp_counter_name():
    content = Path('templates/project_config.yaml').read_text(encoding='utf-8')
    assert 'MCP-counter' in content
