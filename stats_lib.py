# --- R / rpy2 bridge ---
import os
import pingouin as pg
os.environ["R_HOME"] = r"C:\Program Files\R\R-4.3.1"
os.environ["PATH"]  = r"C:\Program Files\R\R-4.3.1\bin\x64;" + os.environ["PATH"]
#import rpy2.robjects.packages as rpackages
#utils = rpackages.importr('utils')
#utils.install_packages('BayesFactor', repos="https://cloud.r-project.org")
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import pandas as pd
import numpy as np
from datetime import datetime

pandas2ri.activate()
_BF_READY = True

try:
    BayesFactor = importr('BayesFactor')
except Exception as e:
    raise RuntimeError(
        "R package 'BayesFactor' is required. In an R console run: install.packages('BayesFactor')"
    ) from e

# ---------- helpers ----------
def _safe_sphericity_and_eps(wide):
    """Return (W, pval, eps_dict) across Pingouin versions."""
    res = pg.sphericity(wide)

    # Any tuple/list length: first is W; pick the last number in [0,1] as p.
    if isinstance(res, (tuple, list)):
        W = float(res[0])
        nums = []
        for x in res[1:]:
            try:
                nums.append(float(x))
            except Exception:
                pass
        p_candidates = [x for x in nums if 0.0 <= x <= 1.0]
        pval = float(p_candidates[-1]) if p_candidates else np.nan
        eps = pg.epsilon(wide)
        return W, pval, dict(eps)

    # Dict-like
    if isinstance(res, dict):
        W = float(res.get('W', np.nan))
        pval = float(res.get('pval', np.nan))
        eps = pg.epsilon(wide)
        return W, pval, dict(eps)

    # DataFrame-like
    if isinstance(res, pd.DataFrame):
        row = res.iloc[0]
        W = float(row.get('W', np.nan))
        pval = float(row.get('pval', row.get('p', np.nan)))
        eps = pg.epsilon(wide)
        return W, pval, dict(eps)

    # Fallback
    raise RuntimeError(f"Unhandled sphericity return type: {type(res)}")

def _bf_rm_anova(df_long, dv, within, subject):
    """BF10 for RM-ANOVA: BF(Subject+Within) / BF(Subject). Robust for your dataset."""
    import numpy as np, pandas as pd
    if not _BF_READY:
        return np.nan
    try:
        use = df_long[[subject, within, dv]].copy()
        use[dv] = pd.to_numeric(use[dv], errors='coerce')
        use = use.dropna()
        # complete subjects across all within levels
        wide = use.pivot_table(index=subject, columns=within, values=dv)
        wide = wide.dropna(axis=0, how='any')
        # drop zero-variance levels (safety)
        keep_cols = [c for c in wide.columns if np.nanvar(wide[c].astype(float)) > 0]
        if len(keep_cols) < 2 or wide.shape[0] < 2:
            return np.nan
        wide = wide[keep_cols]
        # back to long with simple R-safe names
        long = wide.reset_index().melt(id_vars=subject, var_name=within, value_name=dv)
        long = long.rename(columns={subject:'s', within:'w', dv:'y'})
        ro.globalenv['rdf'] = pandas2ri.py2rpy(long)
        ro.r('rdf$s <- factor(rdf$s); rdf$w <- factor(rdf$w)')
        # explicit models; each returns a single BF vs intercept
        ro.r('bf_full <- anovaBF(y ~ w + s, data=rdf, whichRandom="s", progress=FALSE)')
        ro.r('bf_null <- anovaBF(y ~ s,     data=rdf, whichRandom="s", progress=FALSE)')
        ro.r('bf10 <- as.numeric(extractBF(bf_full, onlybf=TRUE)[1]) / as.numeric(extractBF(bf_null, onlybf=TRUE)[1])')
        val = float(ro.r('bf10')[0])
        return val if np.isfinite(val) and val>0 else np.nan
    except Exception:
        return np.nan

def _bf_between_anova(df, dv, between):
    if not _BF_READY:
        return np.nan
    try:
        use = df[[between, dv]].dropna().copy()
        use[between] = use[between].astype(str)
        ro.globalenv['rdf'] = pandas2ri.py2rpy(use)
        ro.r(f'''
            rdf${between} <- factor(rdf${between})
            bf_mod <- anovaBF({dv} ~ {between}, data=rdf, progress=FALSE)
            bf10 <- as.numeric(extractBF(bf_mod, onlybf=TRUE)[1])
        ''')
        return float(ro.r('bf10')[0])
    except Exception:
        return np.nan

def _bf_paired(df_long, dv, within, subject, A, B):
    if not _BF_READY:
        return np.nan
    try:
        sub = df_long[[subject, within, dv]].dropna().copy()
        sub[within] = sub[within].astype(str)
        A = str(A); B = str(B)
        sub = sub[sub[within].isin([A, B])]
        wide = sub.pivot_table(index=subject, columns=within, values=dv)
        if A not in wide.columns or B not in wide.columns:
            return np.nan
        pair = wide.dropna(subset=[A, B])
        if pair.shape[0] < 2:
            return np.nan
        ro.globalenv['xa'] = ro.FloatVector(pair[A].astype(float).values)
        ro.globalenv['xb'] = ro.FloatVector(pair[B].astype(float).values)
        ro.r('bf <- ttestBF(x=xa, y=xb, paired=TRUE, progress=FALSE)')
        ro.r('bf10 <- extractBF(bf, onlybf=TRUE)')
        return float(ro.r('bf10')[0])
    except Exception:
        return np.nan

def _bf_indep(df, dv, group, A, B):
    if not _BF_READY:
        return np.nan
    try:
        sub = df[[group, dv]].dropna().copy()
        sub[group] = sub[group].astype(str)
        A = str(A); B = str(B)
        xa = sub.loc[sub[group]==A, dv].astype(float).values
        xb = sub.loc[sub[group]==B, dv].astype(float).values
        if xa.size < 2 or xb.size < 2:
            return np.nan
        ro.globalenv['xa'] = ro.FloatVector(xa)
        ro.globalenv['xb'] = ro.FloatVector(xb)
        ro.r('bf <- ttestBF(x=xa, y=xb, paired=FALSE, progress=FALSE)')
        ro.r('bf10 <- extractBF(bf, onlybf=TRUE)')
        return float(ro.r('bf10')[0])
    except Exception:
        return np.nan

def _is_df(x): return isinstance(x, pd.DataFrame) and not x.empty

# ---------- main ----------
def anova_analysis(
        data,
        dv='value',
        within=None,                 # repeated-measures factor (optional)
        between=None,                # between-subjects factor (optional)
        subject='sub',
        padjust='fdr_bh',
        alpha=0.05,
        add_signi_stars=True
):
    df = data.copy()

    # validate design
    has_within  = (within  is not None and within  in df.columns)
    has_between = (between is not None and between in df.columns)
    if not has_within and not has_between:
        raise ValueError("Specify at least one factor: within=... or between=...")

    # ---------- bookkeeping ----------
    N = df[subject].nunique() if subject in df.columns else df.shape[0]

    # WITHIN specifics
    if has_within:
        df[within] = df[within].astype(str)
        levels_w = pd.Index(df[within].unique()).sort_values()
        k = len(levels_w)
        if k < 2:
            raise ValueError(f"Within factor '{within}' must have >=2 levels.")
        wide = df.pivot_table(index=subject, columns=within, values=dv)
    else:
        levels_w, k, wide = pd.Index([]), 1, None

    # BETWEEN specifics
    if has_between:
        df[between] = df[between].astype(str)
        levels_b = pd.Index(df[between].unique()).sort_values()
        g = len(levels_b)
        if g < 2:
            raise ValueError(f"Between factor '{between}' must have >=2 groups.")
    else:
        levels_b, g = pd.Index([]), 1

    # ---------- assumptions ----------
    ass_rows = []

    # Normality by WITHIN level
    if has_within:
        try:
            n_tbl = pg.normality(df, dv=dv, group=within)
            pvals = n_tbl['pval'].astype(float).tolist()
            viol_norm_w = any((p < alpha) for p in pvals if pd.notna(p))
            labels = n_tbl['group'].astype(str).tolist() if 'group' in n_tbl else n_tbl.index.astype(str).tolist()
            ass_rows.append(['Normality (by WITHIN level)', bool(viol_norm_w),
                             '; '.join([f"{g}: p={p:.3g}" for g, p in zip(labels, pvals)])])
        except Exception as e:
            ass_rows.append(['Normality (by WITHIN level)', True, f'check failed ({e})'])

    # Sphericity
    if has_within and k >= 3 and wide is not None and wide.shape[1] >= 3:
        try:
            W, p_spher, eps = _safe_sphericity_and_eps(wide)
            ass_rows.append(['Sphericity (Mauchly, WITHIN)', bool(p_spher < alpha),
                             f"W={W:.3g}, p={p_spher:.3g}, GG={eps.get('gg', np.nan):.3g}, HF={eps.get('hf', np.nan):.3g}"])
        except Exception as e:
            ass_rows.append(['Sphericity (Mauchly, WITHIN)', True, f'sphericity computation failed ({e})'])
    elif has_within:
        ass_rows.append(['Sphericity (Mauchly, WITHIN)', False, 'not applicable (<3 levels)'])

    # Normality by BETWEEN group
    if has_between:
        try:
            n_tblb = pg.normality(df, dv=dv, group=between)
            pvalsb = n_tblb['pval'].astype(float).tolist()
            viol_norm_b = any((p < alpha) for p in pvalsb if pd.notna(p))
            labelsb = n_tblb['group'].astype(str).tolist() if 'group' in n_tblb else n_tblb.index.astype(str).tolist()
            ass_rows.append(['Normality (by BETWEEN group)', bool(viol_norm_b),
                             '; '.join([f"{g}: p={p:.3g}" for g, p in zip(labelsb, pvalsb)])])
        except Exception as e:
            ass_rows.append(['Normality (by BETWEEN group)', True, f'check failed ({e})'])

    # Homogeneity (Levene)
    viol_homo = False
    if has_between:
        try:
            lev = pg.homoscedasticity(data=df[[between, dv]].dropna(), dv=dv, group=between, method='levene')
            p_h = float(lev['pval'].iat[0])
            viol_homo = (p_h < alpha)
            ass_rows.append(['Homogeneity (Levene, BETWEEN)', bool(viol_homo), f'p={p_h:.3g}'])
        except Exception as e:
            ass_rows.append(['Homogeneity (Levene, BETWEEN)', True, f'check failed ({e})'])

    assumptions = pd.DataFrame(ass_rows, columns=['Assumption', 'Violated', 'Details']) if ass_rows else pd.DataFrame(
        {'Assumption':[], 'Violated':[], 'Details':[]})

    # convenience flags
    viol_norm_w = any(assumptions.loc[assumptions['Assumption'].str.contains('WITHIN level'), 'Violated'].tolist()) if has_within else False
    viol_spher  = any(assumptions.loc[assumptions['Assumption'].str.contains('Sphericity'), 'Violated'].tolist()) if has_within else False
    viol_norm_b = any(assumptions.loc[assumptions['Assumption'].str.contains('BETWEEN group'), 'Violated'].tolist()) if has_between else False

    omnibus = None
    posthoc = None
    path = None

    # ---------- WITHIN only ----------
    if has_within and not has_between:
        use_nonparam = (viol_norm_w or (viol_spher and k >= 3))

        # Bayes omnibus
        bf_omni = _bf_rm_anova(df, dv, within, subject)
        bayes_row = pd.DataFrame({
            'Test':['bayes_rm_anova'], 'N':[N], 'Factor':[within], 'BF10':[bf_omni],
            'Note':[f"Model: {dv} ~ {subject} + {within}  vs  {dv} ~ {subject}"]
        })

        if k >= 3:
            if not use_nonparam:
                aov = pg.rm_anova(data=df, dv=dv, within=within, subject=subject, detailed=True, correction=True).copy()
                aov.insert(0, 'Test', 'anova_rm')
                aov.insert(1, 'N', N)
                aov['Factor'] = within
                omnibus = pd.concat([bayes_row, aov], ignore_index=True)

                ph = pg.pairwise_tests(data=df, dv=dv, within=within, subject=subject,
                                       padjust=padjust, parametric=True).copy()
                if add_signi_stars and 'p-corr' in ph:
                    ph['signi'] = ph['p-corr'].apply(lambda x: '***' if x < 0.001 else '**' if x < 0.01 else '*' if x < 0.05 else '')
                ph.insert(0, 'Test', 'anova_rm_post_hoc')
                ph['BF10'] = [ _bf_paired(df, dv, within, subject, A, B) for A, B in zip(ph['A'], ph['B']) ]
                posthoc = ph
                path = 'Parametric RM-ANOVA (WITHIN) + BF10 (omnibus & pairs)'
            else:
                fr = pg.friedman(data=df, dv=dv, within=within, subject=subject).copy()
                fr = fr.rename(columns={'Q':'stat'})  # keep 'dof' and 'p-unc'
                fr.insert(0, 'Test', 'friedman')
                fr.insert(1, 'N', N)
                fr['Factor'] = within
                fr['BF10'] = np.nan
                omnibus = pd.concat([bayes_row, fr], ignore_index=True)

                ph = pg.pairwise_tests(data=df, dv=dv, within=within, subject=subject,
                                       padjust=padjust, parametric=False).copy()
                if add_signi_stars and 'p-corr' in ph:
                    ph['signi'] = ph['p-corr'].apply(lambda x: '***' if x < 0.001 else '**' if x < 0.01 else '*' if x < 0.05 else '')
                ph.insert(0, 'Test', 'friedman_post_hoc')
                ph['BF10'] = [ _bf_paired(df, dv, within, subject, A, B) for A, B in zip(ph['A'], ph['B']) ]
                posthoc = ph
                path = 'Nonparametric (WITHIN): Friedman + Wilcoxon; BF10 for pairwise'

        else:  # k == 2
            a, b = list(levels_w)
            if not use_nonparam:
                t = pg.ttest(x=df.loc[df[within]==a, dv], y=df.loc[df[within]==b, dv], paired=True)
                t_row = pd.DataFrame({
                    'Test':['ttest_paired'], 'N':[N], 'A':[a], 'B':[b],
                    'T':[t['T'].iat[0]], 'dof':[t['dof'].iat[0]], 'p-val':[t['p-val'].iat[0]],
                    'CI95%':[str((t['CI95%'][0][0], t['CI95%'][0][1]))], 'cohen-d':[t['cohen-d'].iat[0]],
                    'BF10':[_bf_paired(df, dv, within, subject, a, b)], 'Factor':[within]
                })
                omnibus = pd.concat([bayes_row, t_row], ignore_index=True)
                path = 'Parametric (WITHIN, 2 levels): paired t-test + BF10'
            else:
                w = pg.wilcoxon(x=df.loc[df[within]==a, dv].values, y=df.loc[df[within]==b, dv].values)
                w_row = pd.DataFrame({
                    'Test':['wilcoxon_signed_rank'], 'N':[N], 'A':[a], 'B':[b],
                    'W-val':[w['W-val'].iat[0]], 'p-val':[w['p-val'].iat[0]],
                    'RBC':[w['RBC'].iat[0]] if 'RBC' in w.columns else [np.nan],
                    'effsize':[w['RBC'].iat[0]] if 'RBC' in w.columns else [np.nan],
                    'BF10':[_bf_paired(df, dv, within, subject, a, b)], 'Factor':[within]
                })
                omnibus = pd.concat([bayes_row, w_row], ignore_index=True)
                path = 'Nonparametric (WITHIN, 2 levels): Wilcoxon + BF10'

    # ---------- BETWEEN only ----------
    elif has_between and not has_within:
        viol_any = viol_norm_b or (viol_homo and g >= 3)
        # Bayes omnibus
        bf_omni = _bf_between_anova(df, dv, between)
        bayes_row = pd.DataFrame({
            'Test':['bayes_between_anova'], 'N':[N], 'Factor':[between], 'BF10':[bf_omni],
            'Note':[f"Model: {dv} ~ {between}  vs  intercept"]
        })

        if g >= 3:
            if not viol_any:
                aov = pg.anova(data=df, dv=dv, between=between, detailed=True).copy()
                aov.insert(0, 'Test', 'anova_between')
                aov.insert(1, 'N', N)
                aov['Factor'] = between
                omnibus = pd.concat([bayes_row, aov], ignore_index=True)

                ph = pg.pairwise_tests(data=df, dv=dv, between=between, padjust=padjust, parametric=True).copy()
                if add_signi_stars and 'p-corr' in ph:
                    ph['signi'] = ph['p-corr'].apply(lambda x: '***' if x < 0.001 else '**' if x < 0.01 else '*' if x < 0.05 else '')
                ph.insert(0, 'Test', 'anova_between_post_hoc')
                ph['BF10'] = [ _bf_indep(df, dv, between, A, B) for A, B in zip(ph['A'], ph['B']) ]
                posthoc = ph
                path = 'Parametric one-way ANOVA (BETWEEN) + BF10 (omnibus & pairs)'
            else:
                kw = pg.kruskal(data=df, dv=dv, between=between).copy()
                kw = kw.rename(columns={'H':'stat'})  # keep 'p-unc' and 'ddof1'
                kw.insert(0, 'Test', 'kruskal_wallis')
                kw.insert(1, 'N', N)
                kw['Factor'] = between
                kw['BF10'] = np.nan
                omnibus = pd.concat([bayes_row, kw], ignore_index=True)

                ph = pg.pairwise_tests(data=df, dv=dv, between=between, padjust=padjust, parametric=False).copy()
                if add_signi_stars and 'p-corr' in ph:
                    ph['signi'] = ph['p-corr'].apply(lambda x: '***' if x < 0.001 else '**' if x < 0.01 else '*' if x < 0.05 else '')
                ph.insert(0, 'Test', 'kruskal_between_post_hoc')
                ph['BF10'] = [ _bf_indep(df, dv, between, A, B) for A, B in zip(ph['A'], ph['B']) ]
                posthoc = ph
                path = 'Nonparametric (BETWEEN): Kruskal–Wallis + Mann–Whitney; BF10 for pairwise'
        else:  # g == 2
            a, b = list(levels_b)
            if not viol_any and not viol_homo:
                t = pg.ttest(x=df.loc[df[between]==a, dv], y=df.loc[df[between]==b, dv], paired=False)
                t_row = pd.DataFrame({
                    'Test':['ttest_independent'], 'N':[N], 'A':[a], 'B':[b],
                    'T':[t['T'].iat[0]], 'dof':[t['dof'].iat[0]], 'p-val':[t['p-val'].iat[0]],
                    'CI95%':[str((t['CI95%'][0][0], t['CI95%'][0][1]))], 'hedges':[t['hedges'].iat[0]],
                    'BF10':[_bf_indep(df, dv, between, a, b)], 'Factor':[between]
                })
                omnibus = pd.concat([bayes_row, t_row], ignore_index=True)
                path = 'Parametric (BETWEEN, 2 groups): independent t-test + BF10'
            else:
                w = pg.mwu(x=df.loc[df[between]==a, dv], y=df.loc[df[between]==b, dv])
                w_row = pd.DataFrame({
                    'Test':['mann_whitney'], 'N':[N], 'A':[a], 'B':[b],
                    'U-val':[w['U-val'].iat[0]], 'p-val':[w['p-val'].iat[0]],
                    'RBC':[w['RBC'].iat[0]] if 'RBC' in w.columns else [np.nan],
                    'effsize':[w['RBC'].iat[0]] if 'RBC' in w.columns else [np.nan],
                    'BF10':[_bf_indep(df, dv, between, a, b)], 'Factor':[between]
                })
                omnibus = pd.concat([bayes_row, w_row], ignore_index=True)
                path = 'Nonparametric (BETWEEN, 2 groups): Mann–Whitney + BF10'

    # ---------- MIXED (WITHIN + BETWEEN) ----------
    else:
        ma = pg.mixed_anova(data=df, dv=dv, within=within, between=between, subject=subject).copy()
        ma.insert(0, 'Test', 'mixed_anova')
        ma.insert(1, 'N', N)
        # try to present factors consistently
        if 'Within' in ma.columns:
            ma.rename(columns={'Within':'Within_Factor', 'Between':'Between_Factor'}, inplace=True)
        omnibus = ma

        # Simple-effects post hocs
        ph = pg.pairwise_tests(data=df, dv=dv, within=within, subject=subject,
                               between=between, padjust=padjust, parametric=True).copy()
        if add_signi_stars and 'p-corr' in ph:
            ph['signi'] = ph['p-corr'].apply(lambda x: '***' if x < 0.001 else '**' if x < 0.01 else '*' if x < 0.05 else '')
        ph.insert(0, 'Test', 'mixed_post_hoc')

        # Add BF10 for pure within or pure between contrasts; leave NaN for interactions
        bf_list = []
        for _, r in ph.iterrows():
            wlab = r.get('within', np.nan)
            blab = r.get('between', np.nan)
            if pd.notna(wlab) and pd.isna(blab):
                bf_list.append(_bf_paired(df, dv, within, subject, r['A'], r['B']))
            elif pd.notna(blab) and pd.isna(wlab):
                bf_list.append(_bf_indep(df, dv, between, r['A'], r['B']))
            else:
                bf_list.append(np.nan)
        ph['BF10'] = bf_list
        posthoc = ph
        path = 'Mixed ANOVA (WITHIN+BETWEEN): omnibus via pg.mixed_anova; simple-effects post hocs with BF10 where applicable'

    # ---------- combined + BF note ----------
    base = omnibus if _is_df(omnibus) else pd.DataFrame()
    phdf = posthoc if _is_df(posthoc) else pd.DataFrame()
    if _is_df(base) and _is_df(phdf):
        combined = pd.concat([base, phdf.drop(columns=['signi']) if 'signi' in phdf.columns else phdf],
                             ignore_index=True)
    elif _is_df(base):
        combined = base.copy()
    elif _is_df(phdf):
        combined = phdf.drop(columns=['signi']) if 'signi' in phdf.columns else phdf.copy()
    else:
        combined = pd.DataFrame()

    if not _BF_READY:
        for tbl in (omnibus, posthoc, combined):
            if _is_df(tbl):
                if 'Note' not in tbl.columns:
                    tbl['Note'] = np.nan
                tbl['Note'] = tbl['Note'].fillna('BF10 unavailable (R BayesFactor not loaded)')

    return {
        'path': path,
        'assumptions': assumptions,
        'omnibus': omnibus,
        'posthoc': posthoc,
        'combined': combined
    }


def _fmt_p(p):
    try: p = float(p)
    except: return str(p)
    if np.isnan(p): return "NA"
    if p < 1e-4: return "<1e-4"
    if p < 0.001: return f"{p:.1e}"
    return f"{p:.3f}"

def _stars(p):
    try: p = float(p)
    except: return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""

def _bf_strength(x):
    # returns label and bucket for 1-sided BF
    if x >= 100:   return "decisive", 5
    if x >= 30:    return "very strong", 4
    if x >= 10:    return "strong", 3
    if x >= 3:     return "moderate", 2
    if x >= 1:     return "anecdotal", 1
    return "NA", 0

def _fmt_bf_dual(bf10):
    """Return a compact string showing evidence for H1 and H0.
       Examples:
         113  -> 'BF10=113 (strong for H1), BF01=0.00885'
         0.477-> 'BF10=0.477 (moderate for H0; BF01=2.10)'
         ~1   -> 'BF10=1.05 (inconclusive; BF01=0.952)'
    """
    try:
        bf10 = float(bf10)
    except Exception:
        return "BF10=NA"

    if not np.isfinite(bf10) or bf10 <= 0:
        return "BF10=NA"

    # cap for printing
    bf10_txt = ">=1e4" if bf10 >= 1e4 else f"{bf10:.3g}"
    bf01 = 1.0 / bf10
    bf01_txt = "<=1e-4" if bf01 <= 1e-4 else f"{bf01:.5g}"

    if 1/3 < bf10 < 3:
        return f"BF10={bf10_txt} (inconclusive; BF01={bf01_txt})"

    if bf10 >= 3:
        lab, _ = _bf_strength(bf10)
        return f"BF10={bf10_txt} ({lab} for H1), BF01={bf01_txt}"

    # bf10 < 1 -> H0 favored
    lab, _ = _bf_strength(1.0/bf10)
    return f"BF10={bf10_txt} ({lab} for H0; BF01={bf01_txt})"

def _table_as_text(df: pd.DataFrame, maxcolwidth=40):
    if df is None or len(df)==0: return "(no data)"
    df2 = df.copy()
    for c in df2.columns: df2[c] = df2[c].apply(lambda x: str(x)[:maxcolwidth])
    return df2.to_string(index=False, line_width=2000)

def _pick_col(df, names, case_insensitive=True):
    cols = list(df.columns)
    if case_insensitive:
        lower = {str(c).lower(): c for c in cols}
        for n in names:
            if n.lower() in lower: return lower[n.lower()]
    else:
        for n in names:
            if n in cols: return n
    return None

def _num(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else np.nan
    except:
        return np.nan

def _tidy_combined_text(res):
    omni = res.get('omnibus')
    ph   = res.get('posthoc')

    lines = []
    # Omnibus block (handles rm_anova, friedman, bayes row)
    if isinstance(omni, pd.DataFrame) and len(omni):
        factor_col = _pick_col(omni, ['Factor'])
        p_cols = ['p-gg-corr','p-hf-corr','p-unc','p-val','Pr(>F)','p']
        stat_cols = ['F','F Value','stat','T','W-val','Q']
        df1_cols = ['ddof1','df1','num df','df']
        df2_cols = ['ddof2','df2','den df']
        pcol  = next((c for c in omni.columns if str(c) in p_cols), None)
        scol  = next((c for c in omni.columns if str(c) in stat_cols), None)
        df1c  = next((c for c in omni.columns if str(c).lower() in [x.lower() for x in df1_cols]), None)
        df2c  = next((c for c in omni.columns if str(c).lower() in [x.lower() for x in df2_cols]), None)
        bfcol = _pick_col(omni, ['BF10'])

        for _, r in omni.iterrows():
            test = str(r.get('Test','')).strip()
            if not test and isinstance(omni.columns[0], str):
                test = str(r.get(omni.columns[0], ''))
            if test == 'bayes_rm_anova':
                bf = r.get(bfcol, np.nan) if bfcol in omni.columns else np.nan
                fac = str(r.get(factor_col,'')) if factor_col else ''
                lines.append(f"OMNIBUS  bayes_rm_anova [{fac}]: BF10={_fmt_bf_dual(bf)}")
                continue

            fac = str(r.get(factor_col,'')) if factor_col else ''
            stat = r.get(scol, np.nan) if scol else np.nan
            pval = r.get(pcol, np.nan) if pcol else np.nan
            part = f"OMNIBUS  {test}" + (f" [{fac}]" if fac else "")
            if not pd.isna(stat):
                if df1c and df2c and (df1c in omni.columns and df2c in omni.columns):
                    part += f" | stat={_num(stat):.3f} (df={int(_num(r[df1c]))},{int(_num(r[df2c]))})"
                elif df1c and (df1c in omni.columns):
                    part += f" | stat={_num(stat):.3f} (df={int(_num(r[df1c]))})"
                else:
                    part += f" | stat={_num(stat):.3f}"
            if not pd.isna(pval): part += f" | p={_fmt_p(pval)}{_stars(pval)}"
            if bfcol and bfcol in omni.columns and not pd.isna(r[bfcol]):
                part += f" | BF10={_fmt_bf_dual(r[bfcol])}"
            lines.append(part)

    # Posthoc block
    if isinstance(ph, pd.DataFrame) and len(ph):
        Acol = _pick_col(ph, ['A'])
        Bcol = _pick_col(ph, ['B'])
        pcol = _pick_col(ph, ['p-corr','p-unc','p-val','p.adjusted'], case_insensitive=False)
        bfcol = _pick_col(ph, ['BF10'])
        effcol = _pick_col(ph, ['cohen-d','hedges','RBC','effsize'])
        statcol = _pick_col(ph, ['T','U-val','W-val','Z'])

        # Safe sort
        try:
            ph_sorted = ph.sort_values(pcol) if pcol else ph
        except Exception:
            ph_sorted = ph

        for _, r in ph_sorted.iterrows():
            A = r.get(Acol, '')
            B = r.get(Bcol, '')
            p = r.get(pcol, np.nan) if pcol else np.nan
            bf = r.get(bfcol, np.nan) if bfcol else np.nan
            eff= r.get(effcol, np.nan) if effcol else np.nan
            stat = r.get(statcol, np.nan) if statcol else np.nan

            # Fallback inference for A,B when columns are missing but row contains two level-like numbers
            if (not Acol or not Bcol) and isinstance(r.values, np.ndarray):
                nums = [v for v in r.values if isinstance(v, (int,float,np.floating)) and np.isfinite(v)]
                if len(nums)>=2 and (np.isnan(_num(A)) or np.isnan(_num(B)) or A=='' or B==''):
                    A, B = nums[0], nums[1]

            line = f"CONTRAST  {A} vs {B}"
            if not pd.isna(stat): line += f" | stat={_num(stat):.3f}"
            if not pd.isna(p):    line += f" | p={_fmt_p(p)}{_stars(p)}"
            if not pd.isna(eff):  line += f" | effect={_num(eff):.3f}"
            if not pd.isna(bf):   line += f" | BF10={_fmt_bf_dual(bf)}"
            lines.append(line)

    return "\n".join(lines) if lines else "(no combined summary)"

def _ass_line(df, name_pattern, label):
    if not isinstance(df, pd.DataFrame) or df.empty:
        return f"{label}: NA"
    m = df['Assumption'].str.contains(name_pattern, case=False, na=False)
    if not m.any():
        return f"{label}: NA"
    r = df.loc[m].iloc[0]
    status = 'violated' if bool(r.get('Violated', False)) else 'ok'
    det = r.get('Details', '')
    det = '' if (pd.isna(det) or str(det).strip().lower() in ('nan','none','')) else str(det)
    return f"{label}: {status}" + (f" | {det}" if det else "")

def print_within_report(res, show_full=False, return_text=False):
    import pandas as pd, numpy as np
    from datetime import datetime

    def _is_df(x): return isinstance(x, pd.DataFrame) and not x.empty
    def _val(r, cols):
        for c in cols:
            if c in r and pd.notna(r[c]): return r[c]
        return np.nan

    def _fmt_p(p):
        try: p=float(p)
        except: return "p=NA"
        stars = "***" if p<1e-4 else "**" if p<1e-2 else "*" if p<5e-2 else ""
        if p<1e-4: return "p=<1e-4"+stars
        if p<1e-3: return f"p={p:.1e}{stars}"
        if p<1e-2: return f"p={p:.3f}{stars}"
        return f"p={p:.3f}{stars}"

    def _bf_dual(bf10):
        try:
            bf10 = float(bf10)
            if not np.isfinite(bf10) or bf10<=0: return "BF10=NA"
        except: return "BF10=NA"
        bf10_txt = ">=1e4" if bf10>=1e4 else f"{bf10:.3g}"
        bf01 = 1.0/bf10
        bf01_txt = "<=1e-4" if bf01<=1e-4 else f"{bf01:.5g}"
        if 1/3<bf10<3: return f"BF10={bf10_txt} (inconclusive; BF01={bf01_txt})"
        if bf10>=3:   return f"BF10={bf10_txt} (for H1), BF01={bf01_txt}"
        # bf10<1
        lab = "moderate" if bf01>=3 and bf01<10 else "strong" if bf01>=10 and bf01<30 else "very strong" if bf01>=30 and bf01<100 else "decisive" if bf01>=100 else "anecdotal"
        return f"BF10={bf10_txt} ({lab} for H0; BF01={bf01_txt})"

    def _ass_line(df, patt, label):
        if not _is_df(df): return f"{label}: NA"
        m = df['Assumption'].str.contains(patt, case=False, na=False)
        if not m.any(): return f"{label}: NA"
        r = df.loc[m].iloc[0]
        status = 'violated' if bool(r.get('Violated', False)) else 'ok'
        det = r.get('Details', '')
        det = '' if (pd.isna(det) or str(det).strip().lower() in ('nan','none','')) else str(det)
        return f"{label}: {status}" + (f" | {det}" if det else "")

    lines = []
    lines.append("WITHIN/BETWEEN ANOVA REPORT")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    path = res.get('path') or 'NA'
    lines.append(f"Analysis path: {path}")
    lines.append("")

    # Assumptions
    ass = res.get('assumptions')
    lines.append("Assumptions:")
    lines.append("  " + _ass_line(ass, r'Normality.*WITHIN',  'Normality by level'))
    lines.append("  " + _ass_line(ass, r'Sphericity',         'Sphericity (Mauchly)'))
    lines.append("  " + _ass_line(ass, r'Normality.*BETWEEN', 'Normality by group'))
    lines.append("  " + _ass_line(ass, r'Homogeneity',        'Homogeneity (Levene)'))
    lines.append("")

    # Omnibus
    omni = res.get('omnibus')
    if _is_df(omni):
        lines.append("Omnibus test:")
        for _, r in omni.iterrows():
            test = str(r.get('Test','')).strip()
            fac  = r.get('Factor', '')
            fac  = '' if (pd.isna(fac) or str(fac).lower() in ('nan','none','')) else f" [{fac}]"
            if test.startswith('bayes'):
                lines.append(f"OMNIBUS  {test}{fac}: {_bf_dual(r.get('BF10', np.nan))}")
            elif test in ('anova_rm','anova_between','mixed_anova'):
                F = _val(r, ['F','F-value','F_val'])
                p = _val(r, ['p-gg-corr','p-unc','p-GG-corr','p-val','p'])
                d1 = _val(r, ['ddof1','DF1'])
                d2 = _val(r, ['ddof2','DF2'])
                df_part = f"(df={int(d1)},{int(d2)})" if np.isfinite(d1) and np.isfinite(d2) else ""
                lines.append(f"OMNIBUS  {test}{fac} | F={F:.3g} {df_part} | {_fmt_p(p)}")
            elif test=='friedman':
                stat = _val(r, ['stat','Q'])
                p = _val(r, ['p-unc','p-val','p'])
                df1 = _val(r, ['ddof1','dof'])
                df_part = f"(df={int(df1)})" if np.isfinite(df1) else ""
                lines.append(f"OMNIBUS  friedman{fac} | stat={stat:.3g} {df_part} | {_fmt_p(p)}")
            elif test=='kruskal_wallis':
                stat = _val(r, ['stat','H'])
                p = _val(r, ['p-unc','p-val','p'])
                df1 = _val(r, ['ddof1'])
                df_part = f"(df={int(df1)})" if np.isfinite(df1) else ""
                lines.append(f"OMNIBUS  kruskal{fac} | stat={stat:.3g} {df_part} | {_fmt_p(p)}")
            elif test in ('ttest_paired','ttest_independent'):
                T = _val(r, ['T'])
                p = _val(r, ['p-val','p'])
                d = _val(r, ['cohen-d','hedges'])
                A = r.get('A','?'); B = r.get('B','?')
                bf = r.get('BF10', np.nan)
                parts = [f"OMNIBUS  {test}{fac} {A} vs {B} | T={T:.3g}",
                         f"dof={_val(r,['dof'])}",
                         _fmt_p(p)]
                if np.isfinite(d): parts.append(f"effect={d:.3g}")
                parts.append(_bf_dual(bf))
                lines.append(" | ".join(map(str, parts)))
            elif test in ('wilcoxon_signed_rank','mann_whitney'):
                stat = _val(r, ['W-val','U-val','Z','stat'])
                p = _val(r, ['p-val','p'])
                eff = _val(r, ['RBC','effsize'])
                A = r.get('A','?'); B = r.get('B','?')
                bf = r.get('BF10', np.nan)
                parts = [f"OMNIBUS  {test}{fac} {A} vs {B} | stat={stat:.3g}", _fmt_p(p)]
                if np.isfinite(eff): parts.append(f"effect={eff:.3g}")
                parts.append(_bf_dual(bf))
                lines.append(" | ".join(parts))
        lines.append("")
    else:
        lines.append("Omnibus test: NA\n")

    # Posthocs
    ph = res.get('posthoc')
    if _is_df(ph):
        lines.append("Pairwise contrasts:")
        for _, r in ph.iterrows():
            A = r.get('A','?'); B = r.get('B','?')
            stat = _val(r, ['T','W-val','U-val','Z','stat'])
            p = _val(r, ['p-corr','p-unc','p-val','p'])
            eff = _val(r, ['RBC','cohen-d','hedges','effsize','CLES'])
            bf = r.get('BF10', np.nan)
            parts = [f"CONTRAST  {A} vs {B}",
                     f"stat={stat:.3g}" if np.isfinite(stat) else "stat=NA",
                     _fmt_p(p)]
            if np.isfinite(eff): parts.append(f"effect={eff:.3g}")
            parts.append(_bf_dual(bf))
            lines.append(" | ".join(parts))
        lines.append("")
    else:
        lines.append("Pairwise contrasts: NA\n")

    # Combined (brief)
    comb = res.get('combined')
    if _is_df(comb):
        lines.append("Combined summary:")
        # show Bayes rows + first 10 other rows
        show = comb.copy()
        rows = []
        for _, r in show.iterrows():
            test = str(r.get('Test',''))
            if test.startswith('bayes'):
                fac = r.get('Factor','')
                fac = '' if (pd.isna(fac) or str(fac).lower() in ('nan','none','')) else f" [{fac}]"
                rows.append(f"OMNIBUS  {test}{fac}: {_bf_dual(r.get('BF10', np.nan))}")
        # add up to 10 more condensed entries
        cnt = 0
        for _, r in show.iterrows():
            if cnt>=10: break
            test = str(r.get('Test',''))
            if test.startswith('bayes'): continue
            A = r.get('A'); B = r.get('B')
            if pd.notna(A) and pd.notna(B):
                stat = _val(r, ['T','W-val','U-val','Z','stat'])
                p = _val(r, ['p-corr','p-unc','p-val','p'])
                bf = r.get('BF10', np.nan)
                rows.append(f"CONTRAST  {A} vs {B} | stat={stat if np.isfinite(stat) else 'NA'} | {_fmt_p(p)} | {_bf_dual(bf)}")
                cnt += 1
            else:
                # omnibus line
                p = _val(r, ['p-gg-corr','p-unc','p-val','p'])
                if pd.notna(p):
                    rows.append(f"OMNIBUS  {test} | {_fmt_p(p)}")
                    cnt += 1
        lines.extend(rows)
        lines.append("")
    else:
        lines.append("Combined summary: NA\n")

    # Footer
    lines.append("Notes:")
    lines.append("  p-values are adjusted where available.")
    lines.append("  BF10 prints support for H1 or H0 with BF01.")
    txt = "\n".join(lines)

    if show_full:
        # raw dumps appended after a separator
        dumps = []
        for key in ('assumptions','omnibus','posthoc','combined'):
            df = res.get(key)
            if _is_df(df):
                dumps.append(f"\n--- {key.upper()} ---\n{df.to_string(index=False)}")
        txt += "".join(dumps)

    if return_text: return txt
    print(txt)

def _bf10_paired(x, y):
    if not _BF_READY: return np.nan
    try:
        xa = pd.Series(x, dtype=float).dropna().values
        xb = pd.Series(y, dtype=float).dropna().values
        n = min(len(xa), len(xb))
        if n < 2: return np.nan
        ro.globalenv['xa'] = ro.FloatVector(xa[:n])
        ro.globalenv['xb'] = ro.FloatVector(xb[:n])
        ro.r('bf <- ttestBF(x=xa, y=xb, paired=TRUE, progress=FALSE)')
        ro.r('bf10 <- extractBF(bf, onlybf=TRUE)')
        return float(ro.r('bf10')[0])
    except Exception:
        return np.nan

def _bf10_indep(x, y):
    if not _BF_READY: return np.nan
    try:
        xa = pd.Series(x, dtype=float).dropna().values
        xb = pd.Series(y, dtype=float).dropna().values
        if len(xa) < 2 or len(xb) < 2: return np.nan
        ro.globalenv['xa'] = ro.FloatVector(xa)
        ro.globalenv['xb'] = ro.FloatVector(xb)
        ro.r('bf <- ttestBF(x=xa, y=xb, paired=FALSE, progress=FALSE)')
        ro.r('bf10 <- extractBF(bf, onlybf=TRUE)')
        return float(ro.r('bf10')[0])
    except Exception:
        return np.nan

def auto_t(x, y, paired=None, alpha=0.05, alternative='two-sided', add_bayes=True):
    """
    Decide parametric vs nonparametric and run the appropriate test.
    Returns a 1-row DataFrame with: test, paired, parametric, n or n1/n2, stat, dof, p, effect, effect_name, bf10.
    """
    x = pd.Series(x, dtype=float)
    y = pd.Series(y, dtype=float)

    # Pairing heuristic if not specified: equal length -> paired
    if paired is None:
        paired = (len(x) == len(y))

    if paired:
        df = pd.concat([x.rename('x'), y.rename('y')], axis=1).dropna()
        if df.empty or len(df) < 3:
            raise ValueError("Not enough paired observations after dropping NaNs.")
        d = df['x'] - df['y']
        try:
            pnorm = float(pg.normality(d)['pval'].iat[0])
        except Exception:
            pnorm = np.nan

        if not np.isnan(pnorm) and pnorm >= alpha:
            # Paired t-test
            res = pg.ttest(df['x'], df['y'], paired=True, alternative=alternative)
            out = {
                'test': 'ttest_paired',
                'paired': True,
                'parametric': True,
                'n': len(df),
                'stat': float(res['T'].iat[0]),
                'dof': float(res['dof'].iat[0]),
                'p': float(res['p-val'].iat[0]),
                'effect': float(res['cohen-d'].iat[0]),
                'effect_name': 'cohen-d',
                'bf10': _bf10_paired(df['x'], df['y']) if add_bayes else np.nan
            }
        else:
            # Wilcoxon signed-rank
            w = pg.wilcoxon(df['x'].values, df['y'].values, alternative=alternative)
            eff = w['RBC'].iat[0] if 'RBC' in w.columns else np.nan
            out = {
                'test': 'wilcoxon_signed_rank',
                'paired': True,
                'parametric': False,
                'n': len(df),
                'stat': float(w[['W-val','Z']].iloc[0].dropna().values[0]),
                'dof': np.nan,
                'p': float(w['p-val'].iat[0]),
                'effect': float(eff) if pd.notna(eff) else np.nan,
                'effect_name': 'RBC',
                'bf10': _bf10_paired(df['x'], df['y']) if add_bayes else np.nan
            }

    else:
        # Independent samples
        xa = x.dropna()
        xb = y.dropna()
        if len(xa) < 3 or len(xb) < 3:
            raise ValueError("Not enough observations per group after dropping NaNs.")
        try:
            pnx = float(pg.normality(xa)['pval'].iat[0])
            pny = float(pg.normality(xb)['pval'].iat[0])
        except Exception:
            pnx = pny = np.nan
        both_normal = (not np.isnan(pnx) and not np.isnan(pny) and pnx >= alpha and pny >= alpha)

        # Levene on two groups
        try:
            lev = pg.homoscedasticity([xa.values, xb.values], method='levene')
            p_lev = float(lev['pval'].iat[0])
        except Exception:
            p_lev = np.nan
        equal_var = (not np.isnan(p_lev) and p_lev >= alpha)

        if both_normal:
            # Student or Welch
            welch = not equal_var
            res = pg.ttest(xa, xb, paired=False, alternative=alternative, correction=welch)
            out = {
                'test': 'welch_t' if welch else 'ttest_independent',
                'paired': False,
                'parametric': True,
                'n1': len(xa),
                'n2': len(xb),
                'stat': float(res['T'].iat[0]),
                'dof': float(res['dof'].iat[0]),
                'p': float(res['p-val'].iat[0]),
                'effect': float(res['hedges'].iat[0]),
                'effect_name': 'hedges-g',
                'bf10': _bf10_indep(xa, xb) if add_bayes else np.nan
            }
        else:
            # Mann–Whitney U
            w = pg.mwu(xa, xb, alternative=alternative)
            eff = w['RBC'].iat[0] if 'RBC' in w.columns else np.nan
            out = {
                'test': 'mann_whitney',
                'paired': False,
                'parametric': False,
                'n1': len(xa),
                'n2': len(xb),
                'stat': float(w[['U-val','Z']].iloc[0].dropna().values[0]),
                'dof': np.nan,
                'p': float(w['p-val'].iat[0]),
                'effect': float(eff) if pd.notna(eff) else np.nan,
                'effect_name': 'RBC',
                'bf10': _bf10_indep(xa, xb) if add_bayes else np.nan
            }

    return pd.DataFrame([out])


