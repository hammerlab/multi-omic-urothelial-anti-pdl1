from os import getcwd, path
import numpy as np
from cohorts.rounding import float_str
from PIL import Image

def range_formatter(series, use_mean=False, round_to_int=False):
    agg = series.mean() if use_mean else series.median()
    def round_func(val):
        if round_to_int:
            return format(int(round(val)), ",d")
        else:
            return float_str(val)
    return "%s (range %s-%s)" % (round_func(agg), round_func(series.min()), round_func(series.max()))

def percent_formatter(fraction):
    return "%d%%" % round(fraction * 100)

def mann_whitney_formatter(results):
    return "n=%d, Mann-Whitney p=%s" % (len(results.with_condition_series) + len(results.without_condition_series), float_str(results.p_value))

def fishers_exact_formatter(results):
    return "n=%d, Fisher's Exact p=%s" % (len(results.with_condition1_series) + len(results.without_condition1_series), float_str(results.p_value))

def logrank_formatter(results):
    return "n=%d, log-rank p=%s" % (
        len(results.with_condition_series) + len(results.without_condition_series), float_str(results.p_value))

def pearsonr_formatter(results):
    return "n=%d, Pearson r=%s p=%s" % (len(results.series_x), float_str(results.coeff), float_str(results.p_value))

def spearmanr_formatter(results):
    return "n=%d, Spearman rho=%s p=%s" % (len(results.series_x), float_str(results.coeff), float_str(results.p_value))

def bootstrap_mean_formatter(series, q=0.95):
    if q < 0 or q > 1:
        raise ValueError("Invalid q %0.2f, needs to be within [0, 1]" % q)
    q = q * 100
    value = series.mean()
    low = np.percentile(series, (100 - q) / 2.0)
    high = np.percentile(series, (q + ((100 - q) / 2.0)))
    return ("%s, %d%% CI (%s, %s)" % (float_str(value), int(q), float_str(low), float_str(high)))

def hr_posterior_formatter(series, q=0.95, n=None,
                           summary='mean',
                           q_format="{q}% CI ({low}, {high})",
                           stat='HR',
                           stat_format="{stat}={value}",
                           include_p=False, p_compare=0):
    # construct resulting string
    string_parts = list()
    if n:
        string_parts.append("n=%d" % (int(n)))
    # summarize stat if summary is given as a list
    if isinstance(summary, list):
        for summary_type in summary:
            if summary_type == 'mean': 
                value = series.mean()
            elif summary_type == 'median':
                value = series.median()
            else:
                raise ValueError("invalid summary %s; needs to be either None, 'mean', or 'median'" % summary_type)
            string_parts.append(stat_format.format(stat='{} ({})'.format(stat, summary_type),
                                                   value=float_str(value))
                               )
    # summarize stat if summary is given as a single str
    elif isinstance(summary, str):
        if summary == 'mean':
            value = series.mean()
        elif summary == 'median':
            value = series.median()
        else:
            raise ValueError("invalid summary %s; needs to be either None, 'mean', or 'median'" % summary)
        string_parts.append(stat_format.format(stat=stat, value=float_str(value)))
    # add posterior interval if q is given
    if q:
        if q < 0 or q > 1:
            raise ValueError("Invalid q %0.2f, needs to be within [0, 1]" % q)
        q = q * 100
        low = np.percentile(series, (100 - q) / 2.0)
        high = np.percentile(series, (q + ((100 - q) / 2.0)))
        string_parts.append(q_format.format(q=int(q), low=float_str(low), high=float_str(high)))
    # compute 'bayesian p-value', ie probability that HR <> 0
    if include_p:
        prob = 1-np.mean(series >= p_compare)
        direction = '>'
        if prob >= 0.5:
            prob = 1-prob
            direction = '<'
        string_parts.append('p({}{}{})={}'.format(stat, direction, p_compare, float_str(prob)))
    if len(string_parts) == 0:
        raise ValueError('No summary components specified')
    return ', '.join(string_parts)

def compare_posterior_dist(df, value='exp(beta)', by='expressed', stat='HR',
                           **kwargs):
    among_true = df.loc[df[by] == True, value]
    among_false = df.loc[df[by] == False, value]
    default_args = dict(summary='median',
                        q_format='[{low}, {high}]',
                        q=None,
                        stat=stat,
                        stat_format='{value}',
                        **kwargs)
    return 'median {stat} of {true} vs. {false} for {by} versus total'.format(
        stat=stat,
        by=by,
        true=hr_posterior_formatter(series=among_true, **default_args),
        false=hr_posterior_formatter(series=among_false, **default_args)) 

def null_formatter(value):
    return value

def hyper_label_printer(formatter, label, **kwargs):
    if not formatter:
        formatter = null_formatter
    print("{{{%s:%s}}}" % (label, formatter(**kwargs)))

def hyper_figure_label_printer(label):
    print("{{{%s}}}" % label)

def mann_whitney_hyper_label_printer(results, label, split_col="benefit"):
    hyper_figure_label_printer("%s_plot" % label)
    hyper_label_printer(range_formatter, "%s_%s" % (label, split_col),
                        series=results.with_condition_series)
    hyper_label_printer(range_formatter, "%s_no_%s" % (label, split_col),
                        series=results.without_condition_series)
    hyper_label_printer(mann_whitney_formatter, "%s_mw" % label, results=results)

def fishers_exact_hyper_label_printer(results, label, split_col="benefit"):
    hyper_figure_label_printer("%s_plot" % label)
    hyper_label_printer(percent_formatter, "%s_%s" % (label, split_col),
                        fraction=results.with_condition1_series.mean())
    hyper_label_printer(percent_formatter, "%s_no_%s" % (label, split_col),
                        fraction=results.without_condition1_series.mean())
    hyper_label_printer(fishers_exact_formatter, "%s_fishers" % label, results=results)

def survival_hyper_label_printer(results, label):
    hyper_figure_label_printer("%s_plot" % label)
    hyper_label_printer(logrank_formatter, "%s_logrank" % label, results=results)

def pearsonr_hyper_label_printer(results, label):
    hyper_figure_label_printer("%s_plot" % label)
    hyper_label_printer(pearsonr_formatter, "%s_pearsonr" % label, results=results)

def spearmanr_hyper_label_printer(results, label):
    hyper_figure_label_printer("%s_plot" % label)
    hyper_label_printer(spearmanr_formatter, "%s_spearmanr" % label, results=results)

def resize_image_plos(figure_name, new_figure_name, figure_path=getcwd(), height=2625, min_width=789, max_width=2250, dpi=300):
    image = Image.open(path.join(figure_path, figure_name))
    image_width, image_height = image.size
    height_percent = (height / float(image_height))
    width = int((float(image_width) * float(height_percent)))
    assert width >= min_width, "Image is not wide enough"
    assert width < max_width, "Image is too wide"
    new_image = image.resize((width, height), Image.ANTIALIAS)
    new_image.save(path.join(figure_path, new_figure_name), dpi=(dpi, dpi))
