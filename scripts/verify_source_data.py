"""Run the standalone IA analytical reproductions and aggregate checks."""
import argparse
import json

from . import reproduce_gene_claims, reproduce_rupture_decomposition
from .common import check, read, validate_inputs


def verify_aggregates(collective):
    strict = read('strict58_outcomes.tsv').set_index('endpoint', verify_integrity=True)
    breakdown = strict.drop(index='ALL')
    if not breakdown.drop(columns='n').sum(axis=1).equals(breakdown.n):
        raise ValueError('Strict-set endpoint counts failed')
    if not breakdown.sum().equals(strict.loc['ALL']):
        raise ValueError('Strict-set totals failed')
    check(strict.loc['ALL', 'n'], 58)
    check(strict.loc['ALL', 'REPLICATED_CALIBRATED'], 9)

    perm = read('q1_label_draws.tsv')
    if perm.draw.tolist() != list(range(1000)):
        raise ValueError('Joint label-permutation draw IDs failed')
    check(round(perm.T_perm.max(), 3), 0.687)
    check(int((perm.T_perm >= collective['mean_percentile']).sum()), 0)

    genes = read('fig1_A_genes.tsv')
    check(genes.gene.nunique(), 1345)
    check(genes.n_families.sum(), 2792)
    share = genes.n_hw_families.sum() / genes.n_families.sum()
    check(round(share, 2), 0.46)
    cohorts = read('fig1_B_cohorts.tsv')
    check(cohorts.cohort.nunique(), 6)
    if set(cohorts.families_exposed_total) != {87}:
        raise ValueError('Cohort reuse tally failed')
    models = read('fig1_D_matrix.tsv').set_index('model_id', verify_integrity=True)
    if sorted(models.group.value_counts()) != [7, 11, 36, 39]:
        raise ValueError('Model class counts failed')
    if list(models.loc['IAOBJ_X_008', ['specimen_compartment', 's1']]) != ['PERIPHERAL_BLOOD', 'missing']:
        raise ValueError('November 7 MLDGS correction failed')
    annotation = read('fig1_E_effect_summary.tsv').set_index('contrast').loc['IA-CTRL']
    check(annotation.n_genes, 18918)
    check(round(annotation.pearson_r, 2), -0.01)
    check(round(annotation.sign_agreement, 2), 0.50)

    programmes = read('fig3_B_PROGRAMME_OBS_PRED.tsv').set_index('program')
    blood = programmes.loc[['PLATELET', 'NEUTROPHIL_ENRICHED', 'ERYTHROID_HEMOGLOBIN']]
    if not ((blood.obs - blood.pred - blood.resid).abs() < 1e-8).all() or not (blood.resid_lo > 0).all():
        raise ValueError('Fig. 3B blood programme residuals failed')
    check(round(blood.loc['PLATELET', 'resid'], 2), 1.23)

    mapped = read('fig4_A_MODELS.tsv')
    missing = read('fig4_A_NOT_MAPPABLE.tsv')
    check(len(set(mapped.label) | set(missing.label)), 9)
    if dict(zip(missing.label, missing.missing_member)) != {'ICDR-8': 'HMGB1', 'FRG-4': 'ABCC1'}:
        raise ValueError('Non-projectable score identities failed')
    if sorted(mapped.dimension_class.value_counts()) != [1, 2, 4]:
        raise ValueError('Score projection class tally failed')
    transport = read('fig4_B_TRANSPORT.tsv')
    passing = transport[(transport.decision == 'beyond matched null') & (transport.calibration_gate == 'PASS')]
    if sorted(zip(passing.label, passing.cohort)) != [('m6A-4', 'GSE54083'), ('m6A-4', 'GSE75436')]:
        raise ValueError('Unexposed score decision tally failed')

    atlas = read('fig5_atlas_class_counts.tsv').set_index('atlas_class').n_genes
    atlas_classes = ['Compartment-compatible', 'Agrees within expected cell type',
                     'Disagrees (opposite within cell type)', 'Unresolved in cells']
    if atlas.loc[atlas_classes].tolist() != [139, 8, 13, 74]:
        raise ValueError('Supplementary Table S5 concordance counts failed')
    check(atlas.sum(), 1345)
    bounds = read('fig5_SourceData_Fig5D_contractile_survival.tsv').drop_duplicates()
    check(len(bounds), 6)
    rounded_bounds = sorted((100 * bounds.upper_bound_control_like).round(1))
    if rounded_bounds != [1.5, 9.2, 11.3, 11.8, 13.2, 17.0]:
        raise ValueError('November 7 Fig. 5D bound labels failed')

    worlds = read('fig6_FIG6A_WORLDS.tsv')
    sc = worlds[(worlds.category == 'STATE') & (worlds.assay_class == 'sc') & worlds.matchable]
    check(len(sc), 6)
    sc_range = [round(sc.matched_factor.min(), 2), round(sc.matched_factor.max(), 2)]
    if sc_range != [0.34, 0.50]:
        raise ValueError('Fig. 6A scRNA range failed')
    auc = read('fig6_FIG6D1_DISCRIMINATION_AUROC.tsv')
    nnls = auc[(auc.task == 'T1_loss_vs_state') & (auc.regime == 'ORIG') & (auc.condition == 'C1')]
    check(round(nnls.auroc_groupcv_mean.item(), 2), 0.60)
    return {'strict_claims': 58, 'strict_calibrated': 9, 'label_permutation_max': float(perm.T_perm.max()),
            'label_permutation_n': 1000, 'fig1_literature_genes': 1345,
            'fig1_human_wall_fraction': float(share), 'fig4_score_count': 9,
            'fig3_blood_programme_residuals': blood.resid.to_dict(),
            's5_clear_bulk_genes': int(atlas.loc[atlas_classes].sum()),
            's5_compartment_compatible_genes': int(atlas.loc['Compartment-compatible']),
            'fig5_upper_bounds_percent': rounded_bounds, 'fig6_sc_equivalent_scaling_range': sc_range,
            'fig6_nnls_only_discrimination_AUROC': nnls.auroc_groupcv_mean.item()}


def run(workflow='all'):
    validate_inputs()
    report = {}
    if workflow in ('all', 'claims'):
        report['fig2_level_B'] = reproduce_gene_claims.run()
    if workflow in ('all', 'rupture'):
        report['fig3'] = reproduce_rupture_decomposition.run()
    if workflow == 'all':
        report['frozen_checks_level_C'] = verify_aggregates(report['fig2_level_B'])
    return {'status': 'PASS', **report}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workflow', choices=['all', 'claims', 'rupture'], default='all')
    args = parser.parse_args()
    print(json.dumps(run(args.workflow), indent=2, sort_keys=True))
