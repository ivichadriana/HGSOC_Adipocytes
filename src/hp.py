import gzip
import os
import zipfile
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def ungzip_file(gzip_path, dest_path):
    with gzip.open(gzip_path, 'rb') as f_in, open(dest_path, 'wb') as f_out:
        f_out.write(f_in.read())
    print(f"Extracted gzip: {gzip_path} -> {dest_path}")

def unzip_file(zip_path, dest_dir):
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        for member in zip_ref.namelist():
            target_path = os.path.join(dest_dir, member)
            if not os.path.exists(target_path):
                zip_ref.extract(member, dest_dir)
                print(f"Extracted zip: {member} -> {target_path}")
            else:
                print(f"Skipped existing file: {target_path}")

def plot_distributions(props, *,
                       height=8,           # overall figure height (inch)
                       aspect=1.4,         # width-to-height ratio
                       base_font=18, title):      # master font size
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
    long_props = (
        props
        .reset_index(drop=True)
        .melt(var_name="Cell type", value_name="Proportion")
    )

    # ------------------------------------------------------------
    # 2.  Global font sizes (rcParams)
    # ------------------------------------------------------------
    plt.rcParams.update({
        "font.size"       : base_font,
        "axes.titlesize"  : base_font * 1.2,
        "axes.labelsize"  : base_font,
        "xtick.labelsize" : base_font * 0.9,
        "ytick.labelsize" : base_font * 0.9,
        "legend.fontsize" : base_font,
    })

    # ------------------------------------------------------------
    # 3.  Plot
    # ------------------------------------------------------------
    g = sns.catplot(
        data      = long_props,
        x         = "Cell type",
        y         = "Proportion",
        kind      = "box",
        height    = height,
        aspect    = aspect,
        showcaps  = False,
        fliersize = 3,
        palette   = "Set2"
    )
    g.set_xticklabels(rotation=45, ha="right")

    g.fig.subplots_adjust(top=0.88)                     # room for title
    g.fig.suptitle(f"Distribution of Proportions: {title}", weight="bold")

    plt.show()


def group_cts(cell_types_to_use_grouped, name, props):
    """
    Function to group cell types and add them to the props DataFrame.
    """
    props[name] = props[cell_types_to_use_grouped].sum(axis=1).copy()
    return props
    

def corr_mat_sp(props, title="Spearman correlations"):
    corr_all = props.corr(method="spearman")

    plt.figure(figsize=(16,13))
    sns.heatmap(corr_all, cmap="vlag", center=0, square=True,
                cbar_kws=dict(label="Spearman ρ"), linewidths=.5, annot=True, vmin=-1, vmax=1)
    plt.title(title, weight="bold", fontsize=36)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)

    plt.tight_layout()
    plt.show()

def corr_mat_pe(props, title="Pearson correlations"):
    corr_all = props.corr(method="pearson")
    if props.shape[1] > 4:
        plt.figure(figsize=(19,16))
    else:
        plt.figure(figsize=(16,13))

    sns.heatmap(corr_all, cmap="vlag", center=0, square=True,
                cbar_kws=dict(label="Pearson r"), linewidths=.5, annot=True, vmin=-1, vmax=1)
    plt.title(title, weight="bold", fontsize=36)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.tight_layout()
    plt.show()

def open_and_clean_meta(meta_path,  renaming, tissue_dictionary):
    # ------------------------------ clinical columns -----------------
    meta_full = (pd.read_excel(meta_path, sheet_name=0)
                .rename(columns=str.strip)
                .rename(columns=renaming))

    meta_full = meta_full[renaming.values()]

    meta_full["Stage"] = pd.to_numeric(meta_full["Stage"], errors="coerce")

    meta_full["Event"] = meta_full["Event"].astype(int)

    meta_full["Tissue"] = meta_full["Tissue"].map(tissue_dictionary)

    # In no case we'd use Hispanic variable:
    meta_full.drop(columns=["Hispanic"], inplace=True)
    
    # We are removing debulking treatment that includes CA125.
    meta_full.drop(columns=["Debulk", "NeoTx"], inplace=True)

    return meta_full