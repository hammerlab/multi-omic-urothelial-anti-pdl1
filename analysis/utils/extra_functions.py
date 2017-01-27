from varcode.effects.effect_classes import Exonic
from cohorts.functions import first_not_none_param, get_patient_to_mb, no_filter, snv_count, effect_expressed_filter

def exonic_snv_count(row, cohort, filter_fn=None,
              normalized_per_mb=True, **kwargs):
    filter_fn = first_not_none_param([filter_fn, cohort.filter_fn], no_filter)
    normalized_per_mb = first_not_none_param([normalized_per_mb, cohort.normalized_per_mb], False)
    patient_id = row["patient_id"]
    def exonic_filter_fn(filterable_effect):
        if filter_fn is not None:
            return (isinstance(filterable_effect.effect, Exonic) and
                    filter_fn(filterable_effect))
        return isinstance(filterable_effect.effect, Exonic)
    patient_effects = cohort.load_effects(
        patients=[cohort.patient_from_id(patient_id)],
        filter_fn=exonic_filter_fn,
        **kwargs)
    if patient_id in patient_effects:
        count = len(patient_effects[patient_id])
        if normalized_per_mb:
            count /= float(get_patient_to_mb(cohort)[patient_id])
        return count
    return np.nan


def expressed_snv_count(row, cohort, filter_fn=None,
                                 normalized_per_mb=None):
    filter_fn = first_not_none_param([filter_fn, cohort.filter_fn], no_filter)
    normalized_per_mb = first_not_none_param([normalized_per_mb, cohort.normalized_per_mb], False)
    def expressed_filter_fn(filterable_effect):
        assert filter_fn is not None, "filter_fn should never be None, but it is."
        return filter_fn(filterable_effect) and effect_expressed_filter(filterable_effect)
    return snv_count(row, cohort,
                     filter_fn=expressed_filter_fn,
                     normalized_per_mb=normalized_per_mb)

## alternate name of expressed_snv_count
## used in cases when it is useful to match expressed to non-expressed analogs
## (since non-exonic snvs cannot be expressed, this is more correct comparison)
expressed_exonic_snv_count = expressed_snv_count
expressed_exonic_snv_count.__name__ = 'expressed_exonic_snv_count'
