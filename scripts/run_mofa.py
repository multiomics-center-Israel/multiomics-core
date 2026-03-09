#!/usr/bin/env python3
"""
Run MOFA+ (Multi-Omics Factor Analysis) on multiple omics views.

This script provides a command-line interface to run MOFA+ analysis on
multi-omics data provided as CSV files (features x samples format).

Author: Generated for multiomics pipeline
Date: 2026-02-06
"""

import argparse
import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path

try:
    from mofapy2.run.entry_point import entry_point
except ImportError:
    print("Error: mofapy2 package not found. Please install it:")
    print("  pip install mofapy2")
    sys.exit(1)


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run MOFA+ on multi-omics data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python run_mofa.py \\
    --views rna=data/rna.csv prot=data/prot.csv \\
    --outfile results/mofa_model.hdf5 \\
    --outdir results/mofa_outputs \\
    --factors 10 \\
    --seed 42

Input CSV format:
  - Features as rows, samples as columns
  - First row: sample IDs (header)
  - First column: feature IDs
  - Data should be normalized but not centered/scaled
        """
    )

    parser.add_argument(
        "--views",
        nargs="+",
        required=True,
        metavar="VIEW=FILE",
        help="List of view=filepath pairs (e.g., rna=rna.csv prot=prot.csv)"
    )

    parser.add_argument(
        "--outfile",
        type=str,
        required=True,
        help="Output path for trained MOFA model (.hdf5)"
    )

    parser.add_argument(
        "--outdir",
        type=str,
        required=True,
        help="Output directory for CSV files (factors, weights, R2)"
    )

    parser.add_argument(
        "--factors",
        type=int,
        default=10,
        help="Number of factors to train (default: 10)"
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)"
    )

    parser.add_argument(
        "--max_iter",
        type=int,
        default=1000,
        help="Maximum number of iterations (default: 1000)"
    )

    parser.add_argument(
        "--convergence_mode",
        type=str,
        default="fast",
        choices=["fast", "medium", "slow"],
        help="Convergence mode (default: fast)"
    )

    parser.add_argument(
        "--scale_views",
        action="store_true",
        help="Scale views to unit variance (default: False, only center)"
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output"
    )

    return parser.parse_args()


def parse_view_arguments(view_args):
    """
    Parse view=filepath arguments into a dictionary.

    Args:
        view_args: List of strings in format "view_name=filepath"

    Returns:
        Dictionary mapping view names to file paths
    """
    views = {}
    for arg in view_args:
        if "=" not in arg:
            raise ValueError(f"Invalid view argument format: {arg}. Expected format: view_name=filepath")

        view_name, filepath = arg.split("=", 1)

        if not os.path.exists(filepath):
            raise FileNotFoundError(f"View file not found: {filepath}")

        views[view_name] = filepath

    return views


def load_and_transpose_view(filepath, view_name):
    """
    Load a CSV file and transpose it from (Features x Samples) to (Samples x Features).

    Args:
        filepath: Path to CSV file
        view_name: Name of the view (for error messages)

    Returns:
        Transposed pandas DataFrame with samples as rows, features as columns
    """
    print(f"  Loading {view_name} from {filepath}...")

    try:
        # Read CSV with first column as index (feature IDs)
        df = pd.read_csv(filepath, index_col=0)

        # Check for empty dataframe
        if df.empty:
            raise ValueError(f"View {view_name} is empty")

        # Transpose: features (rows) -> columns, samples (columns) -> rows
        df_t = df.T

        print(f"    Shape after transpose: {df_t.shape[0]} samples x {df_t.shape[1]} features")

        # Check for missing values
        n_missing = df_t.isnull().sum().sum()
        if n_missing > 0:
            print(f"    Warning: {n_missing} missing values detected")

        return df_t

    except Exception as e:
        raise RuntimeError(f"Error loading view {view_name}: {str(e)}")


def align_samples(views_dict):
    """
    Find common samples across all views and subset data.

    Args:
        views_dict: Dictionary of view_name -> DataFrame

    Returns:
        Dictionary of view_name -> aligned DataFrame
        List of common sample IDs
    """
    print("\nAligning samples across views...")

    # Get sample IDs from each view
    sample_sets = {name: set(df.index) for name, df in views_dict.items()}

    # Print sample counts
    for name, samples in sample_sets.items():
        print(f"  {name}: {len(samples)} samples")

    # Find intersection
    common_samples = set.intersection(*sample_sets.values())

    if len(common_samples) == 0:
        raise ValueError("No common samples found across all views!")

    print(f"  Common samples: {len(common_samples)}")

    # Check if any view lost all samples
    for name, samples in sample_sets.items():
        lost = len(samples) - len(common_samples)
        if lost == len(samples):
            raise ValueError(f"View {name} has no samples in common with other views!")
        if lost > 0:
            print(f"    Warning: {name} lost {lost} samples not present in other views")

    # Sort sample IDs for consistency
    common_samples = sorted(list(common_samples))

    # Subset all views to common samples in same order
    aligned_views = {}
    for name, df in views_dict.items():
        aligned_views[name] = df.loc[common_samples, :]

    return aligned_views, common_samples


def prepare_mofa_data(aligned_views):
    """
    Prepare data in the format required by MOFA.

    Args:
        aligned_views: Dictionary of view_name -> aligned DataFrame

    Returns:
        List of numpy arrays (one per view)
        List of view names
        List of feature names per view
        Sample names
    """
    print("\nPreparing data for MOFA...")

    data_list = []
    view_names = []
    feature_names = []

    for view_name in sorted(aligned_views.keys()):
        df = aligned_views[view_name]

        # Convert to numpy array and replace NaN with actual NaN (MOFA handles missing data)
        data_array = df.values.astype(np.float32)

        data_list.append(data_array)
        view_names.append(view_name)
        feature_names.append(df.columns.tolist())

        print(f"  {view_name}: {data_array.shape}")

    sample_names = aligned_views[list(aligned_views.keys())[0]].index.tolist()

    return data_list, view_names, feature_names, sample_names


def run_mofa(data_list, view_names, feature_names, sample_names, args):
    """
    Initialize and train MOFA model.

    Args:
        data_list: List of numpy arrays
        view_names: List of view names
        feature_names: List of feature name lists
        sample_names: List of sample names
        args: Command-line arguments

    Returns:
        Trained MOFA model
    """
    print("\n" + "="*60)
    print("Initializing MOFA model")
    print("="*60)

    # Initialize entry point
    ent = entry_point()

    # Set data - MOFA expects a list of matrices, one per view
    # Each matrix should be (samples x features)
    ent.set_data_matrix(
        data_list,
        views_names=view_names,
        groups_names=["group1"],  # Single group
        samples_names=sample_names,
        features_names=feature_names
    )

    # Set data options
    print("\nSetting data options...")
    ent.set_data_options(
        scale_groups=False,  # Don't scale across groups (only one group)
        scale_views=args.scale_views,  # Scale views to unit variance if requested
        center_groups=True  # Center features within each group
    )

    # Set model options
    print("Setting model options...")
    ent.set_model_options(
        factors=args.factors,
        likelihoods=["gaussian"] * len(view_names),  # Gaussian likelihood for continuous data
        spikeslab_weights=True,  # Use spike-and-slab prior for automatic relevance determination
        ard_factors=True,  # Automatic relevance determination for factors
        ard_weights=True  # Automatic relevance determination for weights
    )

    # Set training options
    print("Setting training options...")
    ent.set_train_options(
        iter=args.max_iter,
        convergence_mode=args.convergence_mode,
        dropR2=0.001,  # Drop factors with R2 < 0.1%
        seed=args.seed,
        verbose=args.verbose
    )

    # Build and train
    print("\n" + "="*60)
    print("Building MOFA model")
    print("="*60)
    ent.build()

    print("\n" + "="*60)
    print("Training MOFA model")
    print("="*60)
    print(f"Factors: {args.factors}")
    print(f"Max iterations: {args.max_iter}")
    print(f"Convergence mode: {args.convergence_mode}")
    print(f"Seed: {args.seed}")
    print()

    ent.run()

    print("\n" + "="*60)
    print("Training completed!")
    print("="*60)

    return ent


def save_outputs(model, view_names, feature_names, sample_names, args):
    """
    Save MOFA model and extract outputs to CSV files.

    Args:
        model: Trained MOFA model
        view_names: List of view names
        feature_names: List of feature name lists
        sample_names: List of sample names
        args: Command-line arguments
    """
    print("\n" + "="*60)
    print("Saving outputs")
    print("="*60)

    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)

    # Save model
    print(f"\nSaving model to {args.outfile}...")
    model.save(args.outfile)

    # Get factors (latent variables)
    print("\nExtracting factors...")
    factors = model.model.getExpectations()["Z"]["group1"]  # Shape: (samples, factors)
    factors_df = pd.DataFrame(
        factors,
        index=sample_names,
        columns=[f"Factor{i+1}" for i in range(factors.shape[1])]
    )
    factors_file = os.path.join(args.outdir, "factors.csv")
    factors_df.to_csv(factors_file)
    print(f"  Saved: {factors_file}")
    print(f"  Shape: {factors_df.shape}")

    # Get weights (loadings) for each view
    print("\nExtracting weights...")
    weights = model.model.getExpectations()["W"]

    for i, view_name in enumerate(view_names):
        view_weights = weights[view_name]  # Shape: (features, factors)
        weights_df = pd.DataFrame(
            view_weights,
            index=feature_names[i],
            columns=[f"Factor{j+1}" for j in range(view_weights.shape[1])]
        )
        weights_file = os.path.join(args.outdir, f"weights_{view_name}.csv")
        weights_df.to_csv(weights_file)
        print(f"  Saved: {weights_file}")
        print(f"    Shape: {weights_df.shape}")

    # Get variance explained (R2)
    print("\nExtracting variance explained...")
    try:
        r2_total = model.model.calculate_variance_explained(total=True)
        r2_per_factor = model.model.calculate_variance_explained(total=False)

        # R2 per view (total)
        r2_total_df = pd.DataFrame(r2_total).T
        r2_total_df.index = ["Total_R2"]
        r2_total_file = os.path.join(args.outdir, "variance_explained_total.csv")
        r2_total_df.to_csv(r2_total_file)
        print(f"  Saved: {r2_total_file}")

        # R2 per factor per view
        r2_per_factor_df = pd.DataFrame(r2_per_factor).T
        r2_per_factor_df.index = [f"Factor{i+1}" for i in range(len(r2_per_factor_df))]
        r2_per_factor_file = os.path.join(args.outdir, "variance_explained_per_factor.csv")
        r2_per_factor_df.to_csv(r2_per_factor_file)
        print(f"  Saved: {r2_per_factor_file}")

        # Print summary
        print("\n  Variance explained summary:")
        for view_name in view_names:
            total_r2 = r2_total_df.loc["Total_R2", view_name]
            print(f"    {view_name}: {total_r2:.2f}%")

    except Exception as e:
        print(f"  Warning: Could not extract variance explained: {str(e)}")

    print("\n" + "="*60)
    print("All outputs saved successfully!")
    print("="*60)


def main():
    """Main function."""
    print("="*60)
    print("MOFA+ Multi-Omics Factor Analysis")
    print("="*60)

    # Parse arguments
    args = parse_arguments()

    try:
        # Parse view arguments
        print("\nParsing input views...")
        views_dict_paths = parse_view_arguments(args.views)
        print(f"  Found {len(views_dict_paths)} views: {', '.join(views_dict_paths.keys())}")

        # Load and transpose data
        print("\nLoading data...")
        views_dict = {}
        for view_name, filepath in views_dict_paths.items():
            views_dict[view_name] = load_and_transpose_view(filepath, view_name)

        # Align samples
        aligned_views, common_samples = align_samples(views_dict)

        # Prepare data for MOFA
        data_list, view_names, feature_names, sample_names = prepare_mofa_data(aligned_views)

        # Run MOFA
        model = run_mofa(data_list, view_names, feature_names, sample_names, args)

        # Save outputs
        save_outputs(model, view_names, feature_names, sample_names, args)

        print("\n" + "="*60)
        print("MOFA analysis completed successfully!")
        print("="*60)
        print(f"\nOutputs saved to:")
        print(f"  Model: {args.outfile}")
        print(f"  CSVs:  {args.outdir}/")
        print()

        return 0

    except Exception as e:
        print(f"\nError: {str(e)}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
