import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path


def create_matplotlib_plot(
        rmsx_csv,
        rmsd_csv,
        rmsf_csv,
        interpolate=False,
        triple=False,
        palette="plasma",
        min_val=None,
        max_val=None,
        log_transform=True,
        custom_fill_label="",
        verbose=True,
        window_check=False,
):
    """
    Generate RMSX plots using matplotlib.
    Replaces the previous R/ggplot2 implementation.
    """
    try:
        from IPython.display import display, Image
        has_ipython = True
    except ImportError:
        has_ipython = False

    # Set font to bold for labels and tickers
    plt.rcParams['font.weight'] = 'bold'
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.titleweight'] = 'bold'

    # Verify palette exists
    try:
        plt.get_cmap(palette)
    except ValueError:
        if verbose:
            print(f"Warning: Palette '{palette}' not found. Falling back to 'plasma'.")
        palette = "plasma"

    # 1. Read the main RMSX CSV
    if not os.path.exists(rmsx_csv):
        if verbose:
            print(f"Error: RMSX CSV not found: {rmsx_csv}")
        return False
        
    rmsx_raw = pd.read_csv(rmsx_csv)
    
    if log_transform and verbose:
        print("Log transform requested: assuming CSV already contains log-scaled values.")
        
    chain_ids = rmsx_raw['ChainID'].unique() if 'ChainID' in rmsx_raw.columns else [None]
    
    # Process each chain
    for chain_id in chain_ids:
        if chain_id is not None:
            rmsx_chain = rmsx_raw[rmsx_raw['ChainID'] == chain_id].copy()
            rmsx_chain = rmsx_chain.drop(columns=['ChainID'])
        else:
            rmsx_chain = rmsx_raw.copy()
            
        if 'ResidueID' not in rmsx_chain.columns:
            if verbose:
                print(f"Error: 'ResidueID' column missing in {rmsx_csv}")
            continue
            
        residues = rmsx_chain['ResidueID'].values
        
        # The remaining columns should be our time slices
        slice_cols = [c for c in rmsx_chain.columns if c != 'ResidueID']
        if not slice_cols:
            if verbose:
                print(f"Error: No slice columns found in CSV for ChainID {chain_id}")
            continue
            
        # Try to infer simulation time from filename: e.g. prefix_simname_0.500_ns.csv
        filename = os.path.basename(rmsx_csv)
        try:
            parts = filename.replace('.csv', '').split('_')
            # Look for '_ns' and assume the value before it is the simulation time
            if parts[-1] == 'ns':
                sim_len = float(parts[-2])
            else:
                sim_len = float(parts[-1])
                
            # Snap to the nearest whole integer if it's extremely close (e.g., 99.99 -> 100)
            if abs(sim_len - round(sim_len)) < 0.1:
                sim_len = float(round(sim_len))

        except (ValueError, IndexError):
            if verbose:
                print(f"Warning: Could not parse simulation time from {filename}. Defaulting to number of slices.")
            sim_len = len(slice_cols)
            
        if verbose:
            print(f"Simulation time = {sim_len}")
            
        num_time_points = len(slice_cols)
        step_size = sim_len / num_time_points if num_time_points > 0 else 1.0
        
        # Time centers for the x-axis
        centers = np.linspace(step_size / 2, sim_len - step_size / 2, num_time_points)
        
        # Determine Fill Label
        if custom_fill_label:
            fill_label = custom_fill_label
        elif log_transform:
            fill_label = "Log-\nScaled\nRMSX"
        else:
            fill_label = "RMSX (Å)"
            
        # Data for Heatmap
        # Z has shape (len(residues), len(centers))
        try:
            Z = rmsx_chain[slice_cols].values
        except Exception as e:
            if verbose:
                print(f"Error converting slice columns to values: {e}")
            continue
            
        # Determine vmin, vmax
        vmin = min_val if min_val is not None else np.nanmin(Z)
        vmax = max_val if max_val is not None else np.nanmax(Z)
        
        # Read RMSD and RMSF if needed
        rmsd_data = None
        rmsf_data = None
        if triple or window_check:
            if rmsd_csv and os.path.exists(rmsd_csv):
                rmsd_data = pd.read_csv(rmsd_csv)
            if rmsf_csv and os.path.exists(rmsf_csv):
                rmsf_data = pd.read_csv(rmsf_csv)

        fig = None
        if window_check and rmsd_data is not None:
            fig = plt.figure(figsize=(10, 8))
            gs = GridSpec(4, 1, height_ratios=[1, 1, 1, 3], hspace=0.4)
            
            # 1. RMSD (top)
            ax_rmsd = fig.add_subplot(gs[0])
            if 'Time' in rmsd_data.columns:
                t = rmsd_data['Time'] / 1000.0 if rmsd_data['Time'].max() > 100 else rmsd_data['Time']
            else:
                t = np.linspace(0, sim_len, len(rmsd_data))
            ax_rmsd.plot(t, rmsd_data['RMSD'], color='black')
            ax_rmsd.set_ylabel("RMSD (Å)")
            ax_rmsd.set_xlim(0, sim_len)
            
            # 2. Slice-mean RMSD
            ax_slice_rmsd = fig.add_subplot(gs[1], sharex=ax_rmsd)
            # compute mean RMSD per slice
            slice_indices = np.floor((t / sim_len) * num_time_points).astype(int)
            slice_indices = np.clip(slice_indices, 0, num_time_points - 1)
            mean_rmsd = [rmsd_data['RMSD'][slice_indices == i].mean() for i in range(num_time_points)]
            ax_slice_rmsd.plot(centers, mean_rmsd, color='black')
            ax_slice_rmsd.set_ylabel("Mean RMSD")
            
            # 3. Mean RMSX
            ax_mean_rmsx = fig.add_subplot(gs[2], sharex=ax_rmsd)
            mean_rmsx = np.nanmean(Z, axis=0)
            ax_mean_rmsx.plot(centers, mean_rmsx, color='black')
            ax_mean_rmsx.set_ylabel("Mean RMSX")
            
            # 4. RMSX Heatmap
            ax_heatmap = fig.add_subplot(gs[3], sharex=ax_rmsd)
            if interpolate:
                X, Y = np.meshgrid(centers, residues)
                shading = 'gouraud'
            else:
                X, Y = np.meshgrid(
                    np.linspace(0, sim_len, num_time_points + 1),
                    np.append(residues, residues[-1] + 1) - 0.5
                )
                shading = 'flat'
            # Use pcolormesh
            c = ax_heatmap.pcolormesh(X, Y, Z, cmap=palette, vmin=vmin, vmax=vmax, shading=shading)
            ax_heatmap.set_ylabel("Residue (Index)")
            ax_heatmap.set_xlabel("Time (ns)")
            # Add colorbar
            cbar = plt.colorbar(c, ax=ax_heatmap, aspect=20, pad=0.02)
            cbar.set_label(fill_label)

        elif triple and rmsd_data is not None and rmsf_data is not None:
            fig = plt.figure(figsize=(12, 6))
            gs = GridSpec(2, 2, height_ratios=[1, 3], width_ratios=[4, 1], wspace=0.1, hspace=0.1)
            
            # Top Left: RMSD
            ax_rmsd = fig.add_subplot(gs[0, 0])
            if 'Time' in rmsd_data.columns:
                t = rmsd_data['Time'] / 1000.0 if rmsd_data['Time'].max() > 100 else rmsd_data['Time']
            else:
                t = np.linspace(0, sim_len, len(rmsd_data))
            ax_rmsd.plot(t, rmsd_data['RMSD'], color='black')
            ax_rmsd.set_ylabel("RMSD (Å)")
            ax_rmsd.set_xlim(0, sim_len)
            ax_rmsd.tick_params(labelbottom=False)  # hide x labels
            
            # Bottom Left: Heatmap
            ax_heatmap = fig.add_subplot(gs[1, 0], sharex=ax_rmsd)
            if interpolate:
                X, Y = np.meshgrid(centers, residues)
                shading = 'gouraud'
            else:
                X, Y = np.meshgrid(
                    np.linspace(0, sim_len, num_time_points + 1),
                    np.append(residues, residues[-1] + 1) - 0.5
                )
                shading = 'flat'
            c = ax_heatmap.pcolormesh(X, Y, Z, cmap=palette, vmin=vmin, vmax=vmax, shading=shading)
            ax_heatmap.set_ylabel("Residue (Index)")
            ax_heatmap.set_xlabel("Time (ns)")
            
            # Add colorbar inside the heatmap grid area but specifically for it
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(ax_heatmap)
            cax = divider.append_axes("right", size="2%", pad=0.1)
            cbar = plt.colorbar(c, cax=cax)
            cbar.set_label(fill_label)
            
            # Bottom Right: RMSF
            ax_rmsf = fig.add_subplot(gs[1, 1], sharey=ax_heatmap)
            if 'ChainID' in rmsf_data.columns and chain_id is not None:
                rmsf_chain_data = rmsf_data[rmsf_data['ChainID'] == chain_id]
            else:
                rmsf_chain_data = rmsf_data
            if 'ResidueID' in rmsf_chain_data.columns and 'RMSF' in rmsf_chain_data.columns:
                ax_rmsf.plot(rmsf_chain_data['RMSF'], rmsf_chain_data['ResidueID'], color='black')
                ax_rmsf.set_xlabel("RMSF (Å)")
                ax_rmsf.tick_params(labelleft=False)
            
        else:
            # Heatmap only
            fig, ax_heatmap = plt.subplots(figsize=(10, 5))
            if interpolate:
                X, Y = np.meshgrid(centers, residues)
                shading = 'gouraud'
            else:
                X, Y = np.meshgrid(
                    np.linspace(0, sim_len, num_time_points + 1),
                    np.append(residues, residues[-1] + 1) - 0.5
                )
                shading = 'flat'
            c = ax_heatmap.pcolormesh(X, Y, Z, cmap=palette, vmin=vmin, vmax=vmax, shading=shading)
            ax_heatmap.set_ylabel("Residue (Index)")
            ax_heatmap.set_xlabel("Time (ns)")
            ax_heatmap.set_xlim(0, sim_len)
            cbar = plt.colorbar(c, ax=ax_heatmap)
            cbar.set_label(fill_label)
            
        # Save figure
        if fig is not None:
            cid_str = str(chain_id) if chain_id is not None else "ALL"
            outfile = str(rmsx_csv).replace('.csv', f'_rmsx_plot_chain_{cid_str}.png')
            outfile_transparent = str(rmsx_csv).replace('.csv', f'_rmsx_plot_chain_{cid_str}_transparent.png')
            
            fig.savefig(outfile, dpi=600, bbox_inches='tight')
            fig.savefig(outfile_transparent, dpi=600, bbox_inches='tight', transparent=True)
            
            plt.close(fig)
            if verbose:
                print(f"Saved: {outfile}")
                print(f"Saved: {outfile_transparent}")
                
            # Try to display in IPython if available
            if has_ipython:
                try:
                    display(Image(filename=outfile))
                except Exception:
                    pass
            
    return True
