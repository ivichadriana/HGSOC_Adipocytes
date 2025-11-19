"""
The following script has functions used in the notebooks of the repository.
Functions can be imported with import src.hp as hp (make sure python path includes folder).
"""

import gzip
import os
import zipfile
import pandas as pd
from statsmodels.othermod.betareg import BetaModel
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatter
import matplotlib
from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter
from lifelines.utils import k_fold_cross_validation
import statsmodels.api as sm
from sklearn.model_selection import ParameterGrid


def plot_forest_subset(
    logsumm, variables, title, x_label="Hazard ratio (95% CI)", xlim=None
):
    """
    logsumm: DataFrame from build_loghr_summary
    variables: list of covariate names to plot
    title: plot title
    """
    # subset in the provided order, keep only those that exist
    variables = [v for v in variables if v in logsumm.index]
    if len(variables) == 0:
        print(f"(No variables to plot for: {title})")
        return

    sub = logsumm.loc[variables].copy()
    ypos = np.arange(len(sub))

    sns.set_style("ticks", {"xtick.major.size": 5, "ytick.major.size": 5})
    fig, ax = plt.subplots(figsize=(15, 1.5 * len(sub)))

    # error bars around log(HR)
    ax.errorbar(
        x=sub["HR"],
        y=ypos,
        xerr=[sub["HR"] - sub["lower95"], sub["upper95"] - sub["HR"]],
        fmt="o",
        capsize=10,
        linewidth=4.5,
        elinewidth=4.5,
        capthick=4.5,
        markersize=12,
    )
    # star p-values (same thresholds as your code)
    for i, p in enumerate(sub["p"]):
        if p < 0.05:
            ax.text(
                sub["upper95"].iloc[i] + 0.01,
                ypos[i],
                "*",
                va="center",
                fontsize=50,
                color="red",
            )
        if p < 0.005:
            ax.text(
                sub["upper95"].iloc[i] + 0.01,
                ypos[i],
                "**",
                va="center",
                fontsize=50,
                color="red",
            )
        if p < 0.0005:
            ax.text(
                sub["upper95"].iloc[i] + 0.01,
                ypos[i],
                "***",
                va="center",
                fontsize=50,
                color="red",
            )

    # tiny y padding
    padding = 0.35
    ax.set_ylim(-padding, len(sub) - 1 + padding)
    if xlim:
        ax.set_xlim(xlim)
    ax.set_yticks(ypos)
    ax.set_xscale("log")
    ax.set_yticklabels(sub.index, fontsize=35)

    # Major ticks at 1–2–5 per decade
    ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
    # Label all major ticks, not just powers of 10
    ax.xaxis.set_major_formatter(LogFormatter(base=10, labelOnlyBase=False))
    ax.xaxis.set_minor_locator(LogLocator(base=10, subs=tuple(range(1, 10))))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.tick_params(axis="x", which="major", labelsize=25)
    ax.axvline(1, ls="--", lw=1)
    ax.set_ylabel("", fontsize=35)
    ax.set_xlabel(x_label, fontsize=35, fontweight="bold")
    ax.set_title(title, fontsize=45, fontweight="bold", y=1.05, x=0.4)
    ax.invert_yaxis()
    plt.tight_layout()
    plt.show()


def build_loghr_summary(cph_model, desired_order=None):
    """
    Returns a summary DataFrame with:
      index = covariate names
      columns = ['logHR', 'lower95_log', 'upper95_log', 'p']
    using the model's coefficient (log HR) and its CI.
    """
    summ = cph_model.summary.copy()

    # optional ordering if provided
    if desired_order is not None:
        keep = [v for v in desired_order if v in summ.index]
        # include anything not in desired_order at the end
        keep += [v for v in summ.index if v not in keep]
        summ = summ.loc[keep]

    # find CI columns robustly (lifelines changes names across versions)
    ci_lower_col = [
        c
        for c in summ.columns
        if c.startswith("coef lower") or c.startswith("lower 95%")
    ][0]
    ci_upper_col = [
        c
        for c in summ.columns
        if c.startswith("coef upper") or c.startswith("upper 95%")
    ][0]

    out = pd.DataFrame(index=summ.index)
    out["logHR"] = summ["coef"]
    out["lower95_log"] = summ[ci_lower_col]
    out["upper95_log"] = summ[ci_upper_col]
    # lifelines may use 'p' or another column name for p-values—fall back to any column whose lowercase name is 'p' if needed
    pcol = (
        "p" if "p" in summ.columns else [c for c in summ.columns if c.lower() == "p"][0]
    )
    out["p"] = summ[pcol]
    return out


def meta_summary_table(df):
    rows = []

    # Race
    for cat, val in df["Race"].value_counts(dropna=False).items():
        label = f"Race = {cat}" if not pd.isna(cat) else "Race = Unknown"
        rows.append(["Race", label, val])

    # Age
    rows.extend(
        [
            ["Age at Diagnosis", "N", df["Age"].count()],
            ["Age at Diagnosis", "Mean", round(df["Age"].mean(), 2)],
            ["Age at Diagnosis", "Std", round(df["Age"].std(), 2)],
            ["Age at Diagnosis", "Min", df["Age"].min()],
            ["Age at Diagnosis", "Max", df["Age"].max()],
        ]
    )

    # Vital Status
    vital_counts = df["Event"].map({1: "Deceased", 0: "Alive/Censored"}).value_counts()
    for cat, val in vital_counts.items():
        rows.append(["Vital Status", cat, val])

    # Survival Time
    rows.extend(
        [
            ["Years from diagnosis to last follow up", "N", df["Time_Yrs"].count()],
            [
                "Years from diagnosis to last follow up",
                "Mean",
                round(df["Time_Yrs"].mean(), 2),
            ],
            [
                "Years from diagnosis to last follow up",
                "Std",
                round(df["Time_Yrs"].std(), 2),
            ],
            [
                "Years from diagnosis to last follow up",
                "Min",
                round(df["Time_Yrs"].min(), 2),
            ],
            [
                "Years from diagnosis to last follow up",
                "Max",
                round(df["Time_Yrs"].max(), 2),
            ],
        ]
    )

    # FIGO Stage
    for cat, val in df["Stage"].value_counts(dropna=False).items():
        label = f"Stage {cat}" if not pd.isna(cat) else "Stage Unknown"
        rows.append(["FIGO Stage", label, val])

    # BMI
    rows.extend(
        [
            ["BMI", "N", df["BMI"].count()],
            ["BMI", "Mean", round(df["BMI"].mean(), 2)],
            ["BMI", "Std", round(df["BMI"].std(), 2)],
            ["BMI", "Min", round(df["BMI"].min(), 2)],
            ["BMI", "Max", round(df["BMI"].max(), 2)],
            ["BMI", "Unknown", df["BMI"].isna().sum()],
        ]
    )

    # Residual Disease
    for cat, val in df["Residual"].value_counts(dropna=False).items():
        label = f"Residual {cat}" if not pd.isna(cat) else "Residual Unknown"
        rows.append(["Residual Status", label, val])

    # Adjuvant Therapy
    for cat, val in df["AdjTx"].value_counts(dropna=False).items():
        label = f"Adjuvant {cat}" if not pd.isna(cat) else "Adjuvant Unknown"
        rows.append(["Adjuvant Status", label, val])

    # Build DataFrame
    summary_df = pd.DataFrame(
        rows, columns=["Clinical Variable", "Category/Statistic", "Value"]
    )
    return summary_df


def count_decimal_places(x):
    """
    Returns the number of decimal places in a given number.
    It formats the input with high precision, removes trailing zeros,
    and counts the digits after the decimal point.
    If the number has no decimal part, it returns 0.
    """
    s = format(x, ".16f").rstrip("0")  # high precision, strip trailing zeros
    if "." in s:
        return len(s.split(".")[1])
    else:
        return 0


def boxplots_by_cutoff(
    df: pd.DataFrame,
    feature: str,  # e.g., "BMI" or "Age"
    threshold: float,  # e.g., 30 or 50
    macro_cols: list,  # e.g., fractions
    *,
    flag_name: str | None = None,  # default: f"{feature}_high"
    ylabel: str = "Proportion",
    palette: str = "Set2",
    figsize=(6, 4),
    rotate_x=30,
    title: str | None = None,
    legend_title: str | None = None,
    ax=None,
    add_flag_to_df: bool = False,  # set True if you want the 0/1 column added to df
):
    """
    Creates a binary flag (>= threshold), melts macro_cols for tidy plotting, and draws a boxplot.
    Returns the long (melted) DataFrame and the matplotlib Axes used.
    """
    # Determine flag/labels
    if flag_name is None:
        flag_name = f"{feature}_high"
    if legend_title is None:
        legend_title = f"{feature} \u2265 {threshold}"
    if title is None:
        title = f"Cell group proportions by {feature} class"

    # Build a temporary frame with the flag and macro columns (no mutation unless requested)
    flag_series = (df[feature] >= threshold).astype(int)
    tmp = pd.concat([flag_series.rename(flag_name), df[macro_cols]], axis=1)

    # Optionally persist the flag on the original df
    if add_flag_to_df:
        df[flag_name] = flag_series

    # Tidy (long) table
    long = tmp.melt(id_vars=flag_name, var_name="Group", value_name="Prop")

    # Plot
    sns.set(style="whitegrid", palette=palette)
    created_ax = False
    if ax is None:
        plt.figure(figsize=figsize)
        ax = plt.gca()
        created_ax = True

    sns.boxplot(
        data=long,
        x="Group",
        y="Prop",
        hue=flag_name,
        showcaps=False,
        fliersize=3,
        width=0.6,
        ax=ax,
    )
    ax.set_ylabel(ylabel)
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=rotate_x)
    ax.legend(title=legend_title)
    ax.set_title(title)
    plt.tight_layout()
    if created_ax:
        plt.show()

    return long, ax


def beta_diagnostics(model, frac_name):
    """
    Diagnostics for a statsmodels.othermod.betareg.BetaResults object.
    """
    # ---- pull residuals, fitted, influence ------------------------------
    pearson = model.resid_pearson
    fitted = model.fittedvalues

    infl = model.get_influence()
    leverage = infl.hat_matrix_diag  # h_ii
    cooks_d = infl.cooks_distance[0]

    # studentised Pearson residuals
    stud_res = pearson / np.sqrt(1 - leverage)

    # ---- plotting --------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), dpi=300)

    # 1 ─ Residuals vs. fitted -------------------------------------------
    axes[0].scatter(fitted, pearson, alpha=0.7)
    axes[0].axhline(0, color="gray", lw=1)
    axes[0].set_xlabel("Fitted value", weight="bold")
    axes[0].set_ylabel("Pearson residual", weight="bold")
    axes[0].set_title(f"{frac_name}: Residuals vs. Fitted")

    # 2 ─ QQ-plot of studentised residuals --------------------------------
    sm.qqplot(stud_res, line="45", ax=axes[1])
    axes[1].set_title(f"{frac_name}: QQ-plot (studentised res.)")
    axes[1].set_xlabel("Theoretical Quantiles", weight="bold")
    axes[1].set_ylabel("Sample Quantiles", weight="bold")

    # 3 ─ Leverage vs. residuals with Cook’s D ----------------------------
    sc = axes[2].scatter(leverage, pearson, c=cooks_d, cmap="viridis", alpha=0.7)
    axes[2].set_xlabel("Leverage", weight="bold")
    axes[2].set_ylabel("Pearson residual", weight="bold")
    axes[2].set_title(f"{frac_name}: Influence")
    cb = fig.colorbar(sc, ax=axes[2])
    cb.set_label("Cook’s D")

    plt.tight_layout()
    plt.show()


def tune_cox_penalty(
    data: pd.DataFrame,
    duration_col: str,
    event_col: str,
    search_grid: dict | None = None,
    k: int = 5,
    strata: str | None = None,
    scoring: str = "concordance_index",
) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    """
    Grid-search cross-validation for lifelines.CoxPHFitter penalties.

    Parameters
    ----------
    data : pd.DataFrame
        Input table (rows = subjects, cols = covariates + time + event).
    duration_col, event_col : str
        Column names used by lifelines.
    search_grid : dict
        Keys = hyper-parameter names, values = list/array of candidate values.
        Defaults to a modest ridge grid.
    k : int
        Number of CV folds (default 5).
    strata : str | None
        Optional column name to stratify by (passed to fitter_kwargs).
    scoring : str
        lifelines scoring method (default "concordance_index").

    Returns
    -------
    cv_results : pd.DataFrame
        One row per grid point with mean and SD of the CV metric.
    best_params : dict
        Row (as dict) with the highest mean CV score.
    failures : pd.DataFrame
        Grid points that raised exceptions, with error messages.
    """
    if search_grid is None:
        search_grid = {
            "penalizer": [1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.2, 0.3, 0.5],
            "l1_ratio": [0.0],
        }

    successes, failures = [], []

    for params in ParameterGrid(search_grid):
        cph = CoxPHFitter(**{k: params[k] for k in ("penalizer", "l1_ratio")})
        try:
            cidx = k_fold_cross_validation(
                cph,
                data,
                duration_col=duration_col,
                event_col=event_col,
                k=k,
                scoring_method=scoring,
                fitter_kwargs={"strata": strata} if strata else None,
            )
            successes.append(
                {
                    **params,
                    "c_index_mean": np.mean(cidx),
                    "c_index_sd": np.std(cidx),
                }
            )
        except Exception as e:
            failures.append({**params, "error": str(e)})

    cv_df = (
        pd.DataFrame(successes).sort_values("c_index_mean", ascending=False)
        if successes
        else pd.DataFrame()
    )
    best = cv_df.iloc[0].to_dict() if not cv_df.empty else {}
    fail_df = pd.DataFrame(failures)

    return cv_df, best, fail_df


def get_variable_renaming():
    """
    Provides a dictionary mapping variable names to more readable labels.
    Used to rename dataset columns into standardized, human-friendly terms.
    Returns the renaming dictionary for consistent variable handling.
    """
    renaming = {
        "suid": "ID",
        "refage": "Age",
        "vital_status_fin": "Event",
        "years_extend": "Time_Yrs",
        "tissue": "Tissue",
        "stage": "Stage",
        "race": "Race",
        "dblk_treat": "Debulk",
        "hispanic": "Hispanic",
        "bmi_recent": "BMI",
        "neoadj_treat": "NeoTx",
        "adj_treat": "AdjTx",
        "resdis_treat": "Residual",
    }
    return renaming


def get_tissue_dictionary():
    """
    Provides a dictionary that standardizes tissue descriptions.
    Maps various raw tissue labels to consistent categories like Ovary, Fallopian Tube, Omentum, or Other.
    Helps unify terminology for downstream data analysis or reporting.
    """
    tissue_dictionary = {
        "right ovary": "Ovary",
        "left ovary": "Ovary",
        "ovary": "Ovary",
        "left ovarian mass": "Ovary",
        "right fallopian tube": "Fallopian Tube",
        "left fallopian tube": "Fallopian Tube",
        "fallopian tube": "Fallopian Tube",
        "left fallopian tube and ovary": "Fallopian Tube and Ovary",
        "right fallopian tube and ovary": "Fallopian Tube and Ovary",
        "right ovary and fallopian tube": "Fallopian Tube and Ovary",
        "left ovary and fallopian tube": "Fallopian Tube and Ovary",
        "bilateral tubes and ovaries: tumor including possible ovarian tissue": "Fallopian Tube and Ovary",
        "tubes and ovaries/cancer": "Fallopian Tube and Ovary",
        "fallopian tube and ovary": "Fallopian Tube and Ovary",
        "omentum": "Omentum",
        "omental tumor": "Omentum",
        "omentum (note: ovarian primary)": "Omentum",
        "omentum or peritoneum": "Omentum",
        "peritoneum or omentum": "Omentum",
        "omentum or cul-de-sac implant": "Omentum",
        "umbilicus": "Other",
        "left ovary with adherent omentum": "Other",
        "representative section of mesenteric nodule": "Other",
        "representative sections of tumor": "Other",
        "posterior wall of myometrium": "Other",
        "left ovary and peritoneum": "Other",
        "omentum or ovary": "Other",
        "fallopian tube or ovary": "Other",
        "cervix and colon": "Other",
        "probable adnexal structure with papillary mass": "Other",
        "peritoneum": "Other",
    }
    return tissue_dictionary


def p_to_star(p):
    """
    Converts a p-value into a significance star annotation.
    Returns "***" for p < 0.001, "**" for p < 0.01, "*" for p < 0.05, and an empty string otherwise.
    Useful for quickly marking statistical significance in results.
    """
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def ungzip_file(gzip_path, dest_path):
    """
    Extracts the contents of a gzip-compressed file to a specified destination path.
    Opens the gzip file in binary mode, reads its contents, and writes them to the output file.
    Prints a confirmation message showing the source and destination paths after extraction.
    """
    with gzip.open(gzip_path, "rb") as f_in, open(dest_path, "wb") as f_out:
        f_out.write(f_in.read())
    print(f"Extracted gzip: {gzip_path} -> {dest_path}")


def unzip_file(zip_path, dest_dir):
    """
    Extracts the contents of a zip archive into a destination directory.
    Iterates through each file in the archive, skipping files that already exist.
    Prints whether each file was extracted or skipped.
    """
    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        for member in zip_ref.namelist():
            target_path = os.path.join(dest_dir, member)
            if not os.path.exists(target_path):
                zip_ref.extract(member, dest_dir)
                print(f"Extracted zip: {member} -> {target_path}")
            else:
                print(f"Skipped existing file: {target_path}")


def plot_distributions(
    props,
    *,
    height=10,  # overall figure height (inch)
    aspect=1.4,  # width-to-height ratio
    base_font=18,
    title,
):  # master font size
    """
    Box-plots of all cell-type proportions with large fonts / figure.

    Parameters
    ----------
    props : DataFrame
        Wide table of samples × cell-type proportions (0-1).
    height : float
        Figure height in inches.  Width = height * aspect.
    aspect : float
        Width / height for the catplot.
    base_font : int
        Base font size used everywhere (title, labels, ticks).
    """
    # ------------------------------------------------------------
    # 1.  Re-shape to long format
    # ------------------------------------------------------------
    long_props = props.reset_index(drop=True).melt(
        var_name="Cell type", value_name="Proportion"
    )

    # ------------------------------------------------------------
    # 2.  Global font sizes (rcParams)
    # ------------------------------------------------------------
    plt.rcParams.update(
        {
            "font.size": base_font,
            "axes.labelsize": base_font * 2,
            "xtick.labelsize": base_font * 2,
            "ytick.labelsize": base_font * 2,
            "legend.fontsize": base_font * 2,
        }
    )

    # ------------------------------------------------------------
    # 3.  Plot
    # ------------------------------------------------------------
    g = sns.catplot(
        data=long_props,
        x="Cell type",
        y="Proportion",
        kind="box",
        height=height,
        aspect=aspect,
        showcaps=True,
        fliersize=4,
        palette="Set2",
    )
    g.set_xticklabels(rotation=45, ha="right")

    g.fig.subplots_adjust(top=0.88)  # room for title
    g.fig.suptitle(
        f"Distribution of Proportions: {title}", weight="bold", fontsize=base_font * 2.5
    )

    plt.show()


def group_cts(cell_types_to_use_grouped, name, props):
    """
    Function to group cell types and add them to the props DataFrame.
    """
    props[name] = props[cell_types_to_use_grouped].sum(axis=1).copy()
    return props


def corr_mat_sp(props, title="Spearman correlations"):
    """
    Generates and displays a heatmap of Spearman correlation coefficients for a given dataset.
    Uses Seaborn to create a color-coded matrix with annotations, ranging from -1 to 1.
    Includes customizable title, bold formatting, and large font sizes for readability.
    """
    corr_all = props.corr(method="spearman")

    plt.figure(figsize=(16, 13))
    sns.heatmap(
        corr_all,
        cmap="vlag",
        center=0,
        square=True,
        cbar_kws=dict(label="Spearman ρ"),
        linewidths=0.5,
        annot=True,
        vmin=-1,
        vmax=1,
    )
    plt.title(title, weight="bold", fontsize=36)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)

    plt.tight_layout()
    plt.show()


def corr_mat_pe(props, title="Pearson correlations"):
    """
    Generates and displays a heatmap of Pearson correlation coefficients for a dataset.
    Adjusts figure size based on the number of variables to improve readability.
    Uses Seaborn to plot a color-coded, annotated correlation matrix with values from -1 to 1.
    """
    corr_all = props.corr(method="pearson")
    if props.shape[1] > 4:
        plt.figure(figsize=(19, 16))
    else:
        plt.figure(figsize=(16, 13))

    sns.heatmap(
        corr_all,
        cmap="vlag",
        center=0,
        square=True,
        cbar_kws=dict(label="Pearson r"),
        linewidths=0.5,
        annot=True,
        vmin=-1,
        vmax=1,
    )
    plt.title(title, weight="bold", fontsize=36)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.tight_layout()
    plt.show()


def open_and_clean_meta(meta_path, renaming, tissue_dictionary):
    """
    Loads, cleans, and transforms clinical metadata from an Excel file.
    Parameters
    ----------
    meta_path : str
        Path to the Excel file containing the metadata.
    renaming : dict
        Dictionary mapping original column names to new column names.
    tissue_dictionary : dict
        Dictionary mapping tissue codes to descriptive tissue names.
    Returns
    -------
    pandas.DataFrame
        Cleaned and transformed metadata DataFrame with selected columns, standardized types,
        and mapped tissue names. Unnecessary columns are removed.
    """
    # ------------------------------ clinical columns -----------------
    meta_full = (
        pd.read_excel(meta_path, sheet_name=0)
        .rename(columns=str.strip)
        .rename(columns=renaming)
    )

    meta_full = meta_full[renaming.values()]

    meta_full["Stage"] = pd.to_numeric(meta_full["Stage"], errors="coerce")

    meta_full["Event"] = meta_full["Event"].astype(int)

    meta_full["Tissue"] = meta_full["Tissue"].map(tissue_dictionary)

    # In no case we'd use Hispanic variable, and
    # we are removing debulking treatment that includes CA125.
    meta_full.drop(columns=["Hispanic", "Debulk", "NeoTx"], inplace=True)

    return meta_full


def plot_km(has_info, missing_info, title, duration_col, event_col):
    """
    Plots Kaplan-Meier survival curves comparing two groups: those with treatment information and those missing it.
    Uses the KaplanMeierFitter to estimate survival probabilities over time.
    Customizes labels, titles, and axis formatting for clarity and readability.
    """

    # Initialize Kaplan-Meier fitter
    kmf = KaplanMeierFitter()

    plt.figure(figsize=(14, 12))

    # Plot: patients with treatment info
    kmf.fit(
        durations=has_info[duration_col],
        event_observed=has_info[event_col],
        label="Treatment Data Present",
    )
    kmf.plot_survival_function()

    # Plot: patients missing all treatment info
    kmf.fit(
        durations=missing_info[duration_col],
        event_observed=missing_info[event_col],
        label="Treatment Data Missing",
    )
    kmf.plot_survival_function()

    # Customize plot
    plt.title(title, fontsize=35, fontweight="bold")
    plt.yticks(fontsize=30)
    plt.xticks(fontsize=30)
    plt.xlabel("Time (Days)", fontsize=30)
    plt.ylabel("Survival Probability", fontsize=30)
    plt.grid(True)
    plt.ylim(0, 1)
    plt.legend(fontsize=30, markerscale=2)
    plt.tight_layout()
    plt.show()
