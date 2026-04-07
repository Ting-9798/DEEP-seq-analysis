import scvelo as scv
import argparse
from pkg_resources import parse_version
import os
import scanpy as sc
import matplotlib.pyplot as plt

if parse_version(scv.__version__) < parse_version("0.3.2"):
    raise ValueError("scvelo version needs to be greater than 0.3.2")


def creat_cli(**kwargs):
    parser = argparse.ArgumentParser(description="scvelo analysis")
    parser.add_argument(
        "-m",
        "--mode",
        dest="mode",
        action="store",
        type=str,
        help="Velocity estimation model algorithm, https://scvelo.readthedocs.io/about",
        default="stochastic",
        choices=["stochastic", "dynamical", "deterministic"],
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        action="store",
        type=str,
        help="H5ad file containing clustered data",
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        action="store",
        type=str,
        help="Output directory"
    )

    parser.add_argument(
        "--column",
        required=False,
        action="store",
        type=str,
        help="Column name for cell type annotation",
        default="celltype"
    )

    parser.add_argument(
        "--root-cell",
        required=False,
        type=str,
        help="Specify root cell (can use cell name or cell type)"
    )




    return parser.parse_args(**kwargs)


def check_anndata(adata):
    """
    Check if input meets requirements
    """
    assert (
        adata.layers["spliced"] != adata.X
    ).nnz == 0, "Input anndata object X should be consistent with layers['spliced']"
    return


def get_args():
    """
    For jupyter interactive testing
    """
    args = creat_cli(
        args=[
            "-i","/data/velocyto.h5ad",
            "--column","celltype",
            "-o","./"
        ]
    )
    return args


def run_scvelo():
    args = creat_cli()
    adata = sc.read(args.input)
    check_anndata(adata)
    
    # 
    os.chdir(args.output)

    # 1. Observe the proportion of spliced/unspliced reads
    if args.column in adata.obs.columns:
        scv.pl.proportions(adata, groupby=args.column, save="proportions.pdf")
        scv.pl.proportions(adata, groupby=args.column, save="proportions.svg")

    else:
        raise ValueError("Parameter column needs to be in `adata.obs` columns")

    # 2. Preprocess data
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    # https://github.com/theislab/scvelo/issues/1212
    # Using scv.pp.moments directly reports an error
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    # 3. Estimate RNA velocity
    if args.mode == "dynamical":
        # https://scvelo.readthedocs.io/en/stable/DynamicalModeling.html
        print("Using dynamical model for velocity estimation!")
        scv.tl.recover_dynamics(adata, n_jobs=8)
    
    # Calculated velocity is stored in adata.layers, similar to count matrices
    scv.tl.velocity(adata, mode=args.mode)
    # Generated velocity graph has dimensions n_obs × n_obs
    scv.tl.velocity_graph(adata, n_jobs=8)
    
    if "X_umap" not in adata.obsm_keys():
        sc.tl.umap(adata)
    # Set figure dimensions
    #plt.rcParams["figure.figsize"] = (10, 8)  # width 10, height 8

    custom_palette = [
        "#31CDEE", "#D0F199", "#79BC98", "#3C8487", "#094867", '#E59CC4', "#6666CC",
        "#FEDD81", "#FF9A84", "#9B6194", "#43457B", "#1965B0", "#CCFFCC", "#CCCCFF",
        "#390A4A", "#443880", "#39558B", "#31668D", "#28818E", "#20958B", "#20A487",
        "#48C06E", "#6ECC5A", "#A6DC36", "#F5E24B", "#82E1F6", "#E2F8C3", "#ADD8C0",
        "#89B5B2", "#6C92A0", "#32CBF1", "#FEDA84", "#FF9B84", "#966392", "#094869"
    ]

    # For continuous variables (e.g., pseudotime), custom cmap can be used instead
    # custom_cmap = "plasma"  # Options include "viridis", "inferno", "coolwarm", "Spectral", etc.

    # Plot and save velocity embedding (UMAP)
    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color=args.column,
        palette=custom_palette,   # For categorical color variables
        cmap='RdBu_r',        # For continuous color variables
        legend_loc="right margin",
        legend_fontsize=12,
        legend_fontoutline=1,
        show=False,
        dpi=600,
        save="velocity_umap.pdf"
    )

    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color=args.column,
        palette=custom_palette,
        cmap='RdBu_r',
        legend_loc="right margin",
        legend_fontsize=12,
        legend_fontoutline=1,
        show=False,
        dpi=600,
        save="velocity_umap.svg"
    )

    # 4. Calculate pseudotime based on velocity graph
    scv.tl.velocity_pseudotime(adata)

    # Plot pseudotime
    scv.pl.scatter(
        adata,
        color="velocity_pseudotime",
        cmap='RdBu_r',  # Can change to "viridis", "cool", "Spectral", etc.
        legend_loc="right margin",
        legend_fontsize=12,
        legend_fontoutline=1,
        show=False,
        dpi=600,
        save="velocity_pseudotime.pdf"
    )

    scv.pl.scatter(
        adata,
        color="velocity_pseudotime",
        cmap='RdBu_r',
        legend_loc="right margin",
        legend_fontsize=12,
        legend_fontoutline=1,
        show=False,
        dpi=600,
        save="velocity_pseudotime.svg"
    )

    # 4. Diffusion map and root cell visualization
    sc.tl.diffmap(adata)
    # Select a stable component (use the next best if fewer than 4 components)
    comps = adata.obsm['X_diffmap'].shape[1]
    comp_pair = "1,2" if comps >= 2 else "1,1"
    scv.pl.scatter(
        adata, basis="diffmap", c=[args.column],
        legend_loc="right", components=[comp_pair],
        save="diffmap_by_group.pdf"
    )

    # 6. Cell cycle velocities analysis
    scv.tl.score_genes_cell_cycle(adata)
    scv.pl.scatter(
        adata, color_gradients=['S_score', 'G2M_score'], 
        smooth=True, perc=[5, 95],
        save="velocity_score.pdf"
    )

    # 7. Dynamical model
    scv.tl.recover_dynamics(adata, n_jobs=8)
    # Calculated velocity is stored in adata.layers, similar to count matrices
    scv.tl.velocity(adata, mode='dynamical')
    # Generated velocity graph has dimensions n_obs × n_obs
    scv.tl.velocity_graph(adata, n_jobs=8)
    scv.pl.velocity_embedding_stream(
        adata, basis="umap", save="velocity_umap_dynamical.pdf", color=args.column , dpi=600
    )

    scv.pl.velocity_embedding_stream(
        adata, basis="diffmap", save="velocity_diffmap_dynamical.pdf", color=args.column , dpi=600
    )


    # 5. PAGA plot
    # Error reported
    adata.uns['neighbors']['distances'] = adata.obsp['distances']
    adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']


    scv.tl.paga(adata, groups=args.column,vkey='velocity')
    scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5,save="paga_umap.pdf")



    adata.write_h5ad("scvelo.h5ad")


    return


if __name__ == "__main__":
    run_scvelo()