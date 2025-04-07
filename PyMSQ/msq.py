"""Genetic evaluation of individuals."""
import os
import sys
import math
import numpy as np
import pandas as pd
import pkg_resources
from scipy.stats import norm
from numba import njit



def load_package_data():
    """
    Load genetic data embedded within the package and return them in a dictionary.
    Returns
    -------
    dict of pandas.DataFrame
        Keys:
            - 'chromosome' -> chromosome map data
            - 'marker_effects' -> marker effects for traits
            - 'genotype' -> phased genotype matrix
            - 'group' -> group/classification data
            - 'pedigree' -> pedigree information
    Raises
    ------
    FileNotFoundError
        If any of the data files are not found in the package resources.
    """
    data_files = {
        'chromosome':     'data/chr.txt',
        'marker_effects': 'data/effects.txt',
        'genotype':       'data/phase.txt',
        'group':          'data/group.txt',
        'pedigree':       'data/pedigree.txt'
    }
    data_frames = {}
    for key, file_path in data_files.items():
        try:
            # Read data from the package resource as a stream:
            with pkg_resources.resource_stream(__name__, file_path) as stream:
                # For genotype/pedigree, assume no header; for others, default to standard usage
                if key in ['genotype', 'pedigree']:
                    df = pd.read_table(stream, sep=" ", header=None)
                else:
                    df = pd.read_table(stream, sep=" ")

                data_frames[key] = df
        except FileNotFoundError:
            # If pkg_resources fails to open, it raises a FileNotFoundError
            raise FileNotFoundError(
                f"Data file '{file_path}' not found in package resources. "
                "Check that it exists in 'PyMSQ/data/' and is listed under package_data in setup.py."
            )
    return data_frames


if __name__ == "__main__":
    load_package_data()


def expldmat(gmap, group, **kwargs):
    """
    Calculate the expected within-family LD matrix (R).

    Args:
        gmap (pd.DataFrame): Genetic map information.
        group (pd.DataFrame): Group information.
        **kwargs: Additional keyword arguments.
            mposunit (str): Marker position unit ("cM" or "reco"). 
                Default is "cM".
            method (int): LD calculation method (Bonk = 1, Santos = 2). 
                Default is 1.
            threshold (float): threshold of cM or reco at which independence 
                is assumed. Default is None.
    Returns:
        dict or list: Expected within-family LD matrices for each chromosome
            or group.
    Raises:
        ValueError: If the marker unit is invalid.
    """
    # Retrieve optional keyword arguments
    mposunit = kwargs.get('mposunit', "cM")
    method = kwargs.get('method', 1)
    threshold = kwargs.get('threshold', None)
    # Valid marker position units
    valid_mposunits = {"cM", "cm", "CM", "Cm", "reco", "RECO", "Reco"}
    # Validate marker unit
    if mposunit not in valid_mposunits:
        raise ValueError("Marker unit should be either cM or reco")
    # Get unique group names and chromosome numbers
    probn = group.iloc[:, 1].astype(str).unique().tolist()
    chromos = gmap.iloc[:, 0].unique()
    # Calculate the number of groups
    no_grp = gmap.shape[1] - 3
    
    n_unique_groups = len(probn)
    # checks
    if n_unique_groups > 1 and no_grp == 1:
        print(f"Assuming a group-average map: {n_unique_groups} groups but only one map")
    
    if no_grp > 1 and no_grp != n_unique_groups:
        raise ValueError(
            f"{no_grp} maps!= {n_unique_groups} groups; they must match for distinct group maps."
        )
    # Initialize list to store LD matrices for each group
    mylist = [[] for _ in range(no_grp)]
    # Iterate over each group
    for ngp in range(no_grp):
        grouprecodist = gmap.iloc[:, 3 + ngp]
        # Iterate over each chromosome
        for chrm in chromos:
            # Get marker positions for this chromosome and cast to float32
            mpo = grouprecodist[gmap.iloc[:, 0] == chrm].to_numpy().astype(np.float32)
            if mposunit in ("cM", "cm", "CM", "Cm"):
                # Calculate distance differences using float32 arithmetic
                dist_diff = np.abs(mpo - mpo[:, None])
                if method == 1:  # Bonk et al.'s approach
                    # Pre-calculate denominator as float32
                    factor = np.float32(100)
                    tmp = np.exp(-2 * (dist_diff / factor)) / np.float32(4)
                    if threshold is not None:
                        cutoff = np.exp(-2 * (np.float32(threshold) / factor)) / np.float32(4)
                        tmp[tmp < cutoff] = 0
                elif method == 2:  # Santos et al.'s approach
                    factor = np.float32(200)
                    tmp = -1 * (dist_diff / factor) + np.float32(0.25)
                    if threshold is not None:
                        cutoff = -1 * (np.float32(threshold) / factor) + np.float32(0.25)
                    else:
                        cutoff = -1 * (np.float32(50) / factor) + np.float32(0.25)
                    tmp[tmp < cutoff] = 0
            elif mposunit in ("reco", "RECO"):
                if mpo[0] != 0:
                    raise ValueError(f"First value for reco rate on chr {chrm} isn't zero")
                # Use float32 for indices as well
                idx = np.arange(mpo.size, dtype=np.float32)
                dist_diff = np.abs(idx - idx[:, None])
                if method == 1:  # Bonk et al.'s approach
                    aaa = (1 - (2 * mpo)) / np.float32(4)
                    tmp = aaa[dist_diff.astype(np.int64)]  # indices must be integer
                    if threshold is not None:
                        cutoff = (1 - (2 * np.float32(threshold))) / np.float32(4)
                        tmp[tmp < cutoff] = 0
                elif method == 2:  # Santos et al.'s approach
                    aaa = (-1 * (mpo / np.float32(2))) + np.float32(0.25)
                    tmp = aaa[dist_diff.astype(np.int64)]
                    if threshold is not None:
                        cutoff = -1 * (np.float32(threshold) / np.float32(200)) + np.float32(0.25)
                    else:
                        cutoff = -1 * (np.float32(50) / np.float32(200)) + np.float32(0.25)
                    tmp[tmp < cutoff] = 0
            else:
                raise ValueError("Invalid mposunit value")
            mylist[ngp].append(tmp)
    # Return the LD matrices as a dictionary or list
    if no_grp > 1:
        mylist = dict(zip(probn, mylist))
    return mylist

def formatgen(gmat, progress):
    """
    Format genotype data and output major and minor alleles.
    - If gmat has one column, each cell is a string representing the genotype
      (e.g. "010101").
    - If gmat has two columns, the first is assumed to be an ID and the second
      is assumed to be a genotype string (e.g. "010101") for each row.
    - If gmat has more than two columns, the first column is assumed to be an ID
      and the remaining columns are assumed to contain numeric genotype data.
    The genotype data are converted to np.int8 to reduce memory usage.
    Args:
        gmat (pd.DataFrame or np.ndarray): Genotype data.
        progress (bool): Flag indicating whether to display progress.
    Returns:
        tuple: (genotype_array, allele_array)
            - genotype_array is a NumPy array of dtype np.int8.
            - allele_array is a NumPy array (dtype=np.int8) containing the major
              and minor alleles (and the missing allele if present).
    Raises:
        ValueError: If the marker is not biallelic (or includes more than one missing allele).
    """
    if progress:
        print("Formatting phased haplotypes")
    # Process input based on type and shape, while avoiding unnecessary copying.
    if isinstance(gmat, pd.DataFrame):
        n_cols = gmat.shape[1]
        if n_cols == 1:
            # Single-column: assume each cell is a genotype string.
            genotype_list = [
                np.array(list(val), dtype=np.int8) 
                for val in gmat.iloc[:, 0].values
            ]
            genotype_array = np.vstack(genotype_list)
        elif n_cols == 2:
            # Two columns: first is ID, second is a genotype string.
            genotype_list = [
                np.array(list(val), dtype=np.int8) 
                for val in gmat.iloc[:, 1].values
            ]
            genotype_array = np.vstack(genotype_list)
        else:
            # More than 2 columns: first is ID; the remaining columns are numeric genotype data.
            arr = gmat.iloc[:, 1:].to_numpy(copy=False)
            if arr.dtype != np.int8:
                genotype_array = arr.astype(np.int8, copy=False)
            else:
                genotype_array = arr
    else:
        # If gmat is already a NumPy array, follow the same logic by shape.
        if gmat.ndim == 2:
            n_cols = gmat.shape[1]
            if n_cols == 1:
                # Single column: genotype strings.
                genotype_list = [
                    np.array(list(x), dtype=np.int8) 
                    for x in gmat.ravel()
                ]
                genotype_array = np.vstack(genotype_list)
            elif n_cols == 2:
                # Two columns: first is ID, second is genotype string.
                genotype_list = [
                    np.array(list(x), dtype=np.int8) 
                    for x in gmat[:, 1]
                ]
                genotype_array = np.vstack(genotype_list)
            elif n_cols > 2:
                # Check if the first column is all single-digit integers
                col0 = gmat[:, 0]  # first column
                # A robust check to ensure they're non-negative integers less than 10.
                # We'll also check if dtype is integer-like (or can be interpreted as integer).
                try:
                    col0_int = col0.astype(int)  # attempt integer conversion
                    is_single_digit = np.all((col0_int >= 0) & (col0_int <= 9) & (col0_int == col0))
                except ValueError:
                    is_single_digit = False
            
                if is_single_digit:
                    # If the first column is indeed single-digit genotype data, parse entire gmat as numeric
                    arr = gmat
                else:
                    # Otherwise, drop the first column (treat it as an ID) and parse the remainder
                    arr = gmat[:, 1:]
                # Cast to np.int8 if needed
                if arr.dtype != np.int8:
                    genotype_array = arr.astype(np.int8, copy=False)
                else:
                    genotype_array = arr
        else:
            raise ValueError("gmat should be 2D (rows x columns).")
    print(f"phased genotype data has {genotype_array.shape[0]} rows and {genotype_array.shape[1]} columns")
    # Compute allele frequencies.
    if np.issubdtype(genotype_array.dtype, np.integer) and genotype_array.min() >= 0:
        flat = genotype_array.ravel()
        counts_all = np.bincount(flat)
        allele_vals = np.nonzero(counts_all)[0]
        counts = counts_all[allele_vals]
    else:
        allele_vals, counts = np.unique(genotype_array, return_counts=True)
    # Only supports biallelic (and optionally one missing allele).
    if len(allele_vals) not in [2, 3]:
        raise ValueError("Method only supports biallelic markers (with an optional missing allele)")
    # Determine major and minor alleles (using int8 for consistency).
    major_allele = allele_vals[np.argmax(counts)]
    minor_allele = allele_vals[np.argmin(counts)]
    if progress:
        for a, cnt in zip(allele_vals, counts):
            print(f"Allele {a}: {cnt}")
        print(f"Major allele: {major_allele}")
    # Build the allele array.
    if len(allele_vals) == 2:
        allele_array = np.array([major_allele, minor_allele], dtype=np.int8)
    return genotype_array, allele_array

@njit
def calculate_mspar_spec_me(gmat, alleles):
    """
    Recodes haplotypes for parent-specific marker effects using int16 for intermediate
    computations and returns an int8 array.
    
    Args:
        gmat (np.ndarray): Genetic marker matrix with an even number of rows.
                           Rows are assumed to be paired (first row: paternal haplotype,
                           second row: maternal haplotype, and so on).
        alleles (list): List of alleles (int8), where the first element is the major allele,
                        the second is the minor allele, and optionally the third is the missing allele.
                        
    Returns:
        np.ndarray: A 2D int8 array of recoded parent-specific marker effects with shape 
                    (n_individuals, n_markers), where n_individuals = gmat.shape[0] // 2.
                    
    Notes:
        - For biallelic markers:
            * If both haplotypes equal the major allele, the effect is 0.
            * If both equal the minor allele, the effect is 0.
            * If the haplotypes are (major, minor), the effect is 1.
            * If the haplotypes are (minor, major), the effect is -1.
        - When a third allele is provided, any genotype containing the missing allele is set to 0.
    """
    n_ind = gmat.shape[0] // 2
    n_markers = gmat.shape[1]
    
    # Allocate output array as int8.
    out = np.empty((n_ind, n_markers), dtype=np.int8)
    
    # Convert allele values to int16 to safely perform composite calculations.
    major = np.int16(alleles[0])
    minor = np.int16(alleles[1])
    comp_hom_major = major * 10 + major
    comp_hom_minor = minor * 10 + minor
    comp_het_dom   = major * 10 + minor
    comp_het_rec   = minor * 10 + major
    
    # Determine if a missing allele is provided.
    has_missing = False
    if len(alleles) == 3:
        missing = np.int16(alleles[2])
        has_missing = True
        
    # Loop over individuals and markers.
    for i in range(n_ind):
        for j in range(n_markers):
            # Convert the haplotype values to int16 for composite calculation.
            a = np.int16(gmat[2 * i, j])
            b = np.int16(gmat[2 * i + 1, j])
            composite = a * 10 + b
            
            if composite == comp_hom_major or composite == comp_hom_minor:
                out[i, j] = 0
            elif composite == comp_het_dom:
                out[i, j] = 1
            elif composite == comp_het_rec:
                out[i, j] = -1
            else:
                if has_missing and (a == missing or b == missing):
                    out[i, j] = 0
                else:
                    # Fallback: cast composite to int8 (may underflow if composite is large)
                    out[i, j] = np.int8(composite)
    return out

@njit(fastmath=True)
def trait_matrices(mspar_spec_me, meff):
    num_traits = meff.shape[1]
    mylist = []
    for nt in range(num_traits):
        # Explicitly reshape the 1D array to (1, n_markers)
        trait_matrix = mspar_spec_me * meff[:, nt].reshape(1, -1)
        mylist.append(trait_matrix)
    return mylist

def addimsenumba(gmat, meff, alleles, center):
    """
    Computes parent-specific additive marker effects for each trait.
    
    Args:
        gmat (np.ndarray): Haplotype matrix.
        meff (np.ndarray): Marker effect matrix.
        alleles (list): List of alleles.
        center (bool): Flag indicating whether to center the trait matrix.
        
    Returns:
        list: List of IMSE matrices for each trait.
    """
    # Compute the recoded haplotype matrix (using your efficient int8-based method)
    mspar_spec_me = calculate_mspar_spec_me(gmat, alleles)
    mspar_spec_me = mspar_spec_me.astype(np.float32)        # cast to float32
    # If centering is needed, convert to float32 (for both speed and memory efficiency)
    if center:
        # Compute the column means in float32
        col_means = mspar_spec_me.mean(axis=0, dtype=np.float32)
        # Subtract the means in-place to avoid extra memory allocation
        mspar_spec_me -= col_means
    return trait_matrices(mspar_spec_me, meff)

def progr(itern, total):
    """
    Print progress of a task.
    Args:
        itern (int): Current iteration.
        total (int): Total number of iterations.
    Returns:
        None.
    """
    pct = int(itern / total * 100)
    fill_len = int(50 * itern / total)
    bar = "#" * fill_len + "-" * (50 - fill_len)
    # Write a carriage return but don't add a newline
    sys.stdout.write(f"\rProgress: |{bar}| {pct}% Complete")
    sys.stdout.flush()
    if itern == total:
        # final newline
        sys.stdout.write("\n")

def namesdf(notr, trait_names):
    """
    Create names of dataframe columns for Mendelian co(var).
    Args:
        notr (int): Number of traits.
        trait_names (list): List of trait names.
    Returns:
        list: List of column names for the dataframe.
    """
    colnam = []
    for i in range(notr):
        for trt in range(notr):
            if i == trt:
                colnam.append(trait_names[trt])
            elif i >= trt:
                colnam.append(f"{trait_names[i]}_{trait_names[trt]}")
    return colnam

def mrmmult(temp, exp_ldmat):
    """
    Perform matrix multiplication (MRM' or m'Rm).
    Args:
        temp (np.ndarray): Matrix operand.
        exp_ldmat (np.ndarray): Matrix operand.
    Returns:
        np.ndarray: Result of the matrix product.
    """
    return np.matmul(np.matmul(temp, exp_ldmat), temp.T)

def subindcheck(group, sub_id):
    """
    Check if individuals provided in a pd.DataFrame (sub_id) are present in the group data.
    Args:
        group (pd.DataFrame): Group data.
        sub_id (pd.DataFrame): Individuals' IDs.
    Returns:
        np.ndarray: Array of indices corresponding to the matched individuals.
    Raises:
        ValueError: If individuals' IDs are not provided as a single-column data frame.
        ValueError: If some individual IDs could not be found in the group data.
    """
    # Subset of individual IDs
    sub_id = pd.DataFrame(sub_id)
    if sub_id.shape[1] != 1:
        raise ValueError("Individual IDs (sub_id) should be provided as a single-column data frame")
    group_ids = group.iloc[:, 0].astype(str).tolist()
    sub_ids = sub_id.squeeze().astype(str).tolist()
    # Ensure all individuals in subset are in group dataframe
    indices = [group_ids.index(id_) if id_ in group_ids else None for id_ in sub_ids]
    # Index of individual IDs
    indices = np.array(indices)
    if len(indices) != len(sub_ids):
        raise ValueError("Some individual IDs could not be found in group data")
    return indices

def combine_arrays(arr1, arr2):
    """
    Combine two arrays by alternating their elements.
    Args:
        arr1 (np.ndarray): First array.
        arr2 (np.ndarray): Second array.
    Returns:
        np.ndarray: Combined array with alternating elements from arr1 and arr2.
    """
    combined = np.empty(arr1.size + arr2.size, dtype=arr1.dtype)
    combined[0::2] = arr1
    combined[1::2] = arr2
    return combined

def split_integer(start, end, segments):
    """
    Split an integer range into segments.
    
    Args:
        start (int): Starting value of the range.
        end (int): Ending value of the range.
        segments (int): Number of segments to split the range into.
    Returns:
        list: List of split points.
    Raises:
        ValueError: If the number of segments is less than 2.
    """
    if segments < 2:
        raise ValueError("Number of segments must be at least 2.")
    total_range = end - start
    seg_size = total_range // (segments - 1)
    remainder = total_range % (segments - 1)
    # Each point is computed as: start + i * seg_size + min(i, remainder)
    return [start + i * seg_size + min(i, remainder) for i in range(segments)]

def msvarcov_g_loop(gmat, gmap, meff, exp_ldmat, group, indwt, sub_id, center, progress):
    """
    Compute Mendelian sampling co(variance) using genotype, marker effects, and map data.
    This version avoids float-indexing issues by casting chromosome values to integer
    and ensures exp_ldmat is not a function.
    Parameters:
        gmat : pandas.DataFrame
            Genotype matrix.
        gmap : pandas.DataFrame
            Genetic map information (assumes the first column is chromosome).
        meff : pandas.DataFrame
            Marker effect matrix.
        exp_ldmat : dict or list
            Linkage disequilibrium information.
        group : pandas.DataFrame
            Group information.
        indwt : numpy.ndarray or None
            Index weights for marker effects. (Ignored if only one trait.)
        sub_id : pandas.DataFrame or None
            Subset of individual IDs (single-column DataFrame). If None, all individuals are used.
        center : bool
            Whether to center the recoded marker matrix.
        progress : bool
            Whether to display progress.
    Returns:
        pandas.DataFrame: DataFrame with Mendelian sampling covariance information.
    """
    # ------------------------------------------------------------------------
    # 2) Format genotype data (splitting strings if needed, etc.)
    # ------------------------------------------------------------------------
    gmat1, alleles = formatgen(gmat, progress)
    
    # Convert marker effects to float32 and compute parent-specific marker effects.
    meff = meff.astype(np.float32)  # if meff is a pandas DataFrame/Series
    meff32 = np.asarray(meff.values, dtype=np.float32)
    mselist = addimsenumba(gmat1, meff32, np.array(alleles.tolist(), dtype=np.int8), center)
    # Subset individuals if sub_id is provided.
    if sub_id is not None:
        subset_indices = subindcheck(group, sub_id)
        idn = group.iloc[subset_indices, 0].astype(str).reset_index(drop=True)
        groupsex = group.iloc[subset_indices, 1].astype(str).reset_index(drop=True)
        # Subset each element of mselist in place.
        mselist = [array[subset_indices, :] for array in mselist]
    else:
        idn = group.iloc[:, 0].astype(str).reset_index(drop=True)
        groupsex = group.iloc[:, 1].astype(str).reset_index(drop=True)
    # Validate that the number of LD maps matches (if >1).
    no_grp = gmap.shape[1] - 3
    if no_grp > 1 and (no_grp != len(exp_ldmat)):
        raise ValueError("Number of maps does not match the number of groups in expected LD matrix")
    # Validate index weights.
    if meff.shape[1] > 1 and indwt is None:
        raise ValueError("Please provide index weights")
    if meff.shape[1] == 1 and indwt is not None:
        indwt = None
        print("There is only one marker effect; index weights are ignored")
    if indwt is not None:
        indwt = np.array(indwt, dtype=np.float32)
        if meff.shape[1] != len(indwt):
            raise ValueError("Length of index weights differs from number of marker effects columns")
    notr = meff.shape[1]
    # Check dimensions.
    if gmat1.shape[1] != meff.shape[0]:
        raise ValueError("Number of markers in genotype and marker effects data do not match")
    if gmap.shape[0] != meff.shape[0]:
        raise ValueError("Number of markers (effects) does not match the genetic map")
    # checks
    n_unique_groups = len(group.iloc[:, 1].astype(str).unique().tolist())
    if n_unique_groups > 1 and no_grp == 1:
        print(f"Assuming a group-average map: {n_unique_groups} groups but only one map")
    # Calculate the number of groups
    if no_grp > 1 and no_grp != n_unique_groups:
        raise ValueError(
            f"{no_grp} maps!= {n_unique_groups} groups; they must match for distinct group maps"
        )
    tnames = meff.columns
    no_ind = mselist[0].shape[0]
    prog_milestone = split_integer(0, no_ind - 1, 10)
    if progress:
        progr(0, no_ind)
    # ------------------------------------------------------------------------
    # 3) Convert chromosome values to integers to avoid float indexing
    # ------------------------------------------------------------------------
    chrom_vals = gmap.iloc[:, 0].to_numpy()
    # If they're float, cast to int32. If they're string, you may need a custom mapping.
    if np.issubdtype(chrom_vals.dtype, np.floating):
        chrom_vals = chrom_vals.astype(np.int32)
    elif np.issubdtype(chrom_vals.dtype, np.integer):
        pass  # already int
    else:
        # If they're strings like "chr1", you'd need a custom mapping here
        raise TypeError("Chromosome column is non-numeric. Modify code to handle string-based chromosomes.")
    # Reassign the integer chromosome values back so subsequent logic sees them
    gmap[gmap.columns[0]] = chrom_vals
    # Now extract unique chromosomes as integers
    unique_chromosomes = np.unique(chrom_vals)
    # Build index arrays for each chromosome
    snp_indices = [np.where(chrom_vals == chrm)[0] for chrm in unique_chromosomes]
    # ------------------------------------------------------------------------
    # 4) Allocate output array
    # ------------------------------------------------------------------------
    if notr == 1:
        msvmsc = np.empty((no_ind, 1), dtype=np.float32)
    else:
        # (notr+1)*(notr+2)//2 is the # of lower-tri elements, including diagonal
        mad = (notr + 1) * (notr + 2) // 2
        msvmsc = np.empty((no_ind, mad), dtype=np.float32)
    groupsex_arr = groupsex.to_numpy()
    # ------------------------------------------------------------------------
    # 5) Main loop: build the Mendelian sampling covariance for each individual
    # ------------------------------------------------------------------------
    for i in range(no_ind):
        # Initialize mscov either as a float or a matrix of zeros
        if notr == 1:
            mscov = 0.0  # float32 scalar
        else:
            mscov = np.zeros((notr + 1, notr + 1), dtype=np.float32)
        # Loop over each chromosome, using snp_indices
        for chrm, s_ind in zip(unique_chromosomes, snp_indices):
            # Construct the temp matrix: shape = (notr+1) x (# SNPs on this chromosome)
            if notr > 1:
                temp = np.empty((notr + 1, len(s_ind)), dtype=np.float32)
                for ntr in range(notr):
                    temp[ntr, :] = mselist[ntr][i, s_ind].astype(np.float32, copy=False)
                # Weighted index combination for the (notr+1)-th row
                temp[notr, :] = np.dot(indwt.T, temp[:notr, :])
            else:
                temp = np.empty((1, len(s_ind)), dtype=np.float32)
                temp[0, :] = mselist[0][i, s_ind].astype(np.float32, copy=False)
            # Select the appropriate LD matrix (float32). 
            # If only one group-dist column in gmap => single LD list
            if (gmap.shape[1] - 3) == 1:
                ld = np.array(exp_ldmat[0][chrm - 1], dtype=np.float32, copy=False)
            else:
                # multiple groups => dictionary keyed by groupsex_arr[i],
                # then index by (chrm-1)
                ld = np.array(exp_ldmat[groupsex_arr[i]][chrm - 1], dtype=np.float32, copy=False)
            # Add the (temp * LD) piece for this chromosome
            mscov += mrmmult(temp, ld)
        # Record the result in msvmsc
        if notr > 1:
            # Flatten only the lower triangle (including diagonal)
            msvmsc[i, :] = mscov[np.tril_indices(notr + 1)]
        else:
            msvmsc[i, 0] = mscov
        # Progress bar update
        if progress and i in prog_milestone:
            progr(i + 1, no_ind)
    # Convert final numpy array -> pd.DataFrame
    msvmsc = pd.DataFrame(msvmsc)
    if notr > 1:
        tnames = np.concatenate((tnames, ["AG"]))
        colnam = namesdf(notr + 1, tnames)
        msvmsc.columns = colnam
    else:
        msvmsc.columns = tnames
    # Insert ID and Group columns at the front
    msvmsc.insert(0, "ID", idn, True)
    msvmsc.insert(1, "Group", groupsex, True)
    return msvmsc

def msvarcov(gmat, gmap, meff, exp_ldmat, group, **kwargs):
    """
    Derive Mendelian sampling co(variance) and aggregate genotype.
    Parameters:
    -----------
    gmat : pandas.DataFrame
        Genotype matrix.
    gmap : pandas.DataFrame
        Genetic map information.
    meff : pandas.DataFrame
        Marker effect matrix.
    exp_ldmat : dict
        Linkage disequilibrium information.
    group : pandas.DataFrame
        Group information.
    **kwargs : Optional keyword arguments for customization.
        indwt : numpy.ndarray, optional
            Index weights. Default is None.
        sub_id : pandas.DataFrame or None, optional
            Subset IDs. Default is None.
        center : bool, optional
            Flag indicating whether to center the mspar_spec_me values. Default is False.
        progress : bool, optional
            Flag indicating whether to display the progress bar. Default is False.
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing Mendelian sampling co(variance) information.
    Notes:
    ------
    - If sub_id is None, Mendelian sampling co(variance) will be estimated for all individuals.
    - Otherwise, Mendelian sampling co(variance) will be estimated for the individuals in sub_id.
    """
    # Convert input arguments to DataFrames
    if not isinstance(gmap, pd.DataFrame):
        gmap = pd.DataFrame(gmap)
    if not isinstance(meff, pd.DataFrame):
        meff = pd.DataFrame(meff)
    if not isinstance(group, pd.DataFrame):
        group = pd.DataFrame(group)
    # Get the optional arguments
    indwt = kwargs.get("indwt", None)
    sub_id = kwargs.get("sub_id", None)
    center = kwargs.get("center", False)
    progress = kwargs.get("progress", False)
    # Call the msvarcov_g_loop function
    msvmsc = msvarcov_g_loop(gmat, gmap, meff, exp_ldmat, group, indwt, sub_id, center, progress)
    # other approaches to be developed in the future
    return msvmsc

def cov2corr(cov):
    """
    Convert covariance matrix to correlation matrix.
    
    Parameters:
        cov : numpy.ndarray
            Covariance matrix.
            
    Returns:
        numpy.ndarray
            Correlation matrix.
    """
    # Compute standard deviations from the diagonal.
    std = np.sqrt(np.diag(cov)).astype(np.float32, copy=False)
    # Compute the outer product of std.
    outer_std = np.outer(std, std)
    # Preallocate output array with the same shape and type as cov.
    corr = np.empty_like(cov, dtype=np.float32)
    # Perform element-wise division; avoid division by zero.
    np.divide(cov, outer_std, out=corr, where=outer_std != 0)
    # Set the diagonal to 1.0
    np.fill_diagonal(corr, 1.0)
    return corr

def array2sym(array):
    """
    Convert a 1D array representing the lower triangle (including the diagonal)
    of a symmetric matrix into a standardized correlation matrix, and then
    return the lower triangle (excluding the diagonal) as a 1D array.
    
    This optimized version avoids allocating the full n x n matrix.
    
    Parameters:
        array (numpy.ndarray): 1D input array to convert. Its size must equal n*(n+1)/2 for some n.
    
    Returns:
        numpy.ndarray: 1D array corresponding to the lower triangle (excluding the diagonal)
                       of the correlation matrix.
    """
    dfmsize = array.size
    # Solve n*(n+1)/2 = dfmsize for n
    n = int((-1 + math.sqrt(1 + 8 * dfmsize)) / 2)
    if n * (n + 1) // 2 != dfmsize:
        raise ValueError("Input array size does not correspond to a valid triangular matrix.")
    
    # Compute the indices of the diagonal elements in the flattened lower-triangle.
    # For row i, the diagonal element is at index: (i+1)*(i+2)//2 - 1.
    diag_indices = np.array([ (i+1)*(i+2)//2 - 1 for i in range(n) ], dtype=np.intp)
    variances = array[diag_indices].astype(np.float32, copy=False)
    
    # Get row and column indices for the lower triangle (excluding the diagonal)
    rows, cols = np.tril_indices(n, k=-1)
    # Compute the flat index for each element (i, j) in the lower triangle:
    # The index is given by: i*(i+1)//2 + j
    flat_indices = rows * (rows + 1) // 2 + cols
    
    # Compute correlation for each off-diagonal element: corr(i,j) = cov(i,j)/sqrt(var[i]*var[j])
    corrs = array[flat_indices].astype(np.float32) / np.sqrt(variances[rows] * variances[cols])
    return corrs

def msvarcov_corr(msvmsc):
    """
    Standardize Mendelian sampling co(variance) to correlation.
    
    Parameters:
        msvmsc (pandas.DataFrame): DataFrame containing the Mendelian sampling (co)variance
            and aggregate genotype information created using the `msvarcov` function.
            The DataFrame is assumed to have columns in the order: ID, Group, Trait1, Trait2, ...
    
    Returns:
        pandas.DataFrame: DataFrame containing the standardized correlation values.
    
    Raises:
        ValueError: If the input DataFrame has only three columns (i.e. only one trait).
    """
    if msvmsc.shape[1] == 3:
        raise ValueError("Correlation cannot be derived for a single trait")
    
    # Exclude ID and Group; assume that each row of the remaining data is a flattened lower triangle
    # of a covariance matrix.
    cov_flat = msvmsc.iloc[:, 2:].values.astype(np.float32, copy=False)
    # The column names for the correlation DataFrame: we assume (as before) that names with '_' indicate pair names.
    orig_names = msvmsc.columns.tolist()[2:]
    corr_col_names = [name for name in orig_names if "_" in name]
    
    # Apply the optimized conversion row‐by‐row
    corr_values = np.apply_along_axis(array2sym, 1, cov_flat)
    dfcorr = pd.DataFrame(corr_values, columns=corr_col_names)
    
    # Reattach the ID and Group columns
    result = pd.concat([msvmsc.iloc[:, 0:2].reset_index(drop=True), dfcorr], axis=1)
    return result

def chr_int(chr_of_interest):
    """
    Format the chromosome of interest parameter.
    Parameters:
    ----------
    chr_of_interest : str or list or None
        Chromosome of interest parameter. If "all", all chromosomes are selected.
        If a list of integers is provided, those specific chromosomes are selected.
        If None, no specific chromosome is selected.
    Returns:
    -------
    str or numpy.ndarray or None
        Formatted chromosome of interest parameter.
    Notes:
    ------
    The function takes the chromosome of interest parameter and returns a formatted version.
    If the parameter is already "all" or None, it is returned as is.
    If the parameter is a list of integers, it is converted to a NumPy array of integers.
    """
    if chr_of_interest is not None:
        if chr_of_interest == "all" or chr_of_interest == "All":
            return chr_of_interest
        return np.array([int(i) for i in chr_of_interest])
    return chr_of_interest

def compute_covariance(mselist, exp_ldmat, us_ind, s_ind, i, chrm_value, gmap_shape, pfnp):
    """
    Compute the covariance matrix for a specific chromosome and trait.
    Parameters:
    ----------
    mselist : list of numpy.ndarray
        List of parent-specific marker effects for different traits.
    exp_ldmat : numpy.ndarray
        Array of expected linkage disequilibrium matrices.
    us_ind : numpy.ndarray
        Array of indices for individuals used in the calculation.
    s_ind : int
        Index for the specific snps/markers.
    i : int
        Index for the specific trait.
    chrm_value : int
        Value representing the chromosome.
    gmap_shape : int
        number of maps.
    pfnp : int
        Index value for the expected linkage disequilibrium matrices (chromosome).
    Returns:
    -------
    numpy.ndarray
        Covariance matrix for the specific chromosome and trait.
    Notes:
    ------
    The function computes the covariance matrix by multiplying the genotype matrix (`amat`)
    with the expected linkage disequilibrium matrix (`ldmat`) for a specific chromosome and trait.
    The resulting covariance matrix is returned.
    If the size of the `amat` exceeds 300, matrix multiplication (`@`) is used; otherwise, 
    element-wise multiplication (`*`) is used. The resulting matrix is element-wise absolute
    (`np.abs`) to ensure positive values.
    """
    amat = mselist[i][:, s_ind] if gmap_shape == 1 else mselist[i][us_ind.reshape(-1, 1), s_ind]
    if gmap_shape == 1:
        ldmat = exp_ldmat[0][chrm_value - 1]
    else:
        ldmat = np.array(exp_ldmat[pfnp][chrm_value - 1])
    # if amat.shape[0] > 300:
    #     covtmpx = np.abs(dgmrm(amat, ldmat))
    # else:
    return np.abs(mrmmult(amat, ldmat)).astype(np.float32)


def write_matrices_to_file(chrinterest, stdsim, chrm, trtnam, covtmpx, save=False):
    """
    Write matrices to a file based on the chromosome and trait.
    Parameters:
    ----------
    chrinterest : list or None
        List of chromosomes of interest or None if all chromosomes are of interest.
    stdsim : bool
        Flag indicating whether to compute and save the standardized similarity matrix.
    chrm : int
        Chromosome number.
    trtnam : str
        Name of the trait.
    covtmpx : numpy.ndarray
        Covariance matrix for the given chromosome and trait.
    Returns:
    -------
    None

    Notes:
    ------
    The function writes the covariance matrix and, if specified, the standardized similarity matrix
    to a file based on the chromosome and trait.
    If the `chrinterest` parameter is not None, the function checks if the chromosome is in the 
    list of chromosomes of interest. If it is, the covariance matrix is saved to a file.
    If the `stdsim` parameter is True and the chromosome is in the list of chromosomes of interest
    or the list of chromosomes of interest is 'all', the standardized similarity matrix is 
    also saved to a file.
    The file names are formatted as "{current_working_directory}/Sim mat for {trait_name} 
    chrm {chromosome}.npy" and "{current_working_directory}/Stdsim mat for {trait_name} 
    chrm {chromosome}.npy" for the covariance matrix and standardized similarity matrix,respectively
    """
    # Do not write files if chrinterest is None
    if chrinterest is None:
        return  # Exit function early
    
    if save:
        chrm_str = str(chrm)
        if 'all' in chrinterest or 'All' in chrinterest or chrm_str in chrinterest:
            chrfile1 = f"{os.getcwd()}/Sim mat for {trtnam} chrm {chrm}.npy"
            np.save(chrfile1, covtmpx)
        if stdsim and ('all' in chrinterest or 'All' in chrinterest or chrm_str in chrinterest):
            chrfilec = f"{os.getcwd()}/Stdsim mat for {trtnam} chrm {chrm}.npy"
            np.save(chrfilec, cov2corr(covtmpx))

def subindcheckzyg(group, sub_idz):
    """
    Check the sex and if mate pairs provided in sub_idz are in the group data.
    Parameters:
    -----------
    group : pd.DataFrame
        Grouping information.
    sub_idz : pd.DataFrame
        Subset of individuals with mate pairs.
    Returns:
    --------
    mal1 : np.ndarray
        Indices of males in the group data.
    fem1 : np.ndarray
        Indices of females in the group data.
    idn : np.ndarray
        Subset of individuals with mate pairs.
    probn : list
        Group classes of males and females.
    Raises:
    -------
    ValueError
        If the subset ID does not have two columns or if the subset ID is not present in 
            the group data,
        if the group class in the subset ID is not unique to males or females, or 
            if multiple sexes are detected in the data.
    Notes:
    ------
    The function checks if the provided subset of individuals with mate pairs is present in the 
    group data.
    It verifies the sex of individuals and ensures that the group class in the subset ID is 
    unique to males and females.
    The subset ID should have two columns, one for each group (e.g., Male, Female). The function
    compares the IDs
    in the subset with the IDs in the group data to determine the indices of males and females.
    The function returns the indices of males and females, the subset of individuals with mate 
    pairs, and the group
    classes of males and females.
    Example:
    ---------
    group:
    +-----+-------+
    |  ID | Group |
    +-----+-------+
    |  M1 |   A   |
    |  M2 |   A   |
    |  F1 |   B   |
    |  F2 |   B   |
    +-----+-------+
    sub_idz:
    +-----+-----+
    |  M1 |  F1 |
    |  M2 |  F2 |
    +-----+-----+
    The function will return:
    mal1 = [0, 1]
    fem1 = [2, 3]
    idn = [['M1', 'F1'], ['M2', 'F2']]
    probn = ['A', 'B']
    """
    sub_idz = pd.DataFrame(sub_idz)
    if sub_idz.shape[1] != 2:
        raise ValueError("Please provide two columns, one for each group (e.g., Male, Female)")
    group_ids = group.iloc[:, 0].astype(str).tolist()
    id_to_index = {id_: idx for idx, id_ in enumerate(group_ids)}
    sub_idz = sub_idz.reset_index(drop=True).squeeze()
    mal_ids = sub_idz.iloc[:, 0].astype(str).tolist()
    fem_ids = sub_idz.iloc[:, 1].astype(str).tolist()
    mal1 = [id_to_index[x] for x in mal_ids if x in id_to_index]
    fem1 = [id_to_index[x] for x in fem_ids if x in id_to_index]
    if len(mal1) != len(mal_ids) or len(fem1) != len(fem_ids):
        raise ValueError("Subset ID not found in group data")
    mal1 = np.array(mal1)
    fem1 = np.array(fem1)
    male_group_classes = group.iloc[mal1, 1]
    female_group_classes = group.iloc[fem1, 1]
    if len(male_group_classes.unique()) != 1:
        raise ValueError("Group class in sub_idz is not unique to ID of males")
    if len(female_group_classes.unique()) != 1:
        raise ValueError("Group class in sub_idz is not unique to ID of females")
    idn = sub_idz.reset_index(drop=True).to_numpy()
    mgp = male_group_classes.unique().tolist()
    fgp = female_group_classes.unique().tolist()
    if len(mgp) > 1 or len(fgp) > 1:
        raise ValueError("Multiple sexes detected in data")
    probn = [mgp[0], fgp[0]]
    return mal1, fem1, idn, probn

def simmat(gmat, gmap, meff, group, exp_ldmat, **kwargs):
    """
    Compute similarity matrix for gametes or zygotes.
    Parameters:
    -----------
    gmat : ndarray or DataFrame
        Genotype data.
    gmap : ndarray or DataFrame
        Genetic map information.
    meff : ndarray or DataFrame
        Marker effects.
    group : ndarray or DataFrame
        Grouping information.
    exp_ldmat : list or dict
        Expected within-family LD matrix.
    **kwargs : Optional keyword arguments for customization.
        sub_id : ndarray or DataFrame, optional
            Subset of individuals for subsetting. Default is None.
        indwt : ndarray or None, optional
            Index weights for marker effects. If None, index weights are ignored. Default is None.
        chrinterest : ndarray or list or int, optional
            Chromosome(s) of interest. Default is None.
        save : bool, optional
            Flag indicating whether to save the similarity matrix to a file. Default is False.
        stdsim : bool, optional
            Flag indicating whether to standardize the similarity matrix. Default is False.
        center : bool, optional
            Flag indicating whether to center the similarity matrix. Default is False.
        progress : bool, optional
            Flag indicating whether to enable progress reporting. Default is False.
    Returns:
    --------
    multgrpcov : list or dict
        List or dictionary of similarity matrices computed based on groups and traits.
    Raises:
    -------
    ValueError
        If the expected within-family LD matrix is not a list or dict.
    Notes:
    ------
    The function computes a similarity matrix for gametes or zygotes based on the provided
    genotype data,
    genetic map information, marker effects, and grouping information.
    The expected within-family LD matrix should be precomputed using the function 'expldmat'.
    Optional keyword arguments can be used for customization, such as subsetting individuals,
    providing index weights,
    specifying chromosome of interest, saving the similarity matrix to a file, standardizing 
    the similarity matrix,
    centering the similarity matrix, and enabling progress reporting.
    If the 'sub_id' parameter is not provided or has one column, the function computes the 
    similarity matrix for gametes.
    If the 'sub_id' parameter has two columns, the function computes  similarity matrix for zygotes.
    Example:
    ---------
    gmat: ndarray or DataFrame
    gmap: ndarray or DataFrame
    meff: ndarray or DataFrame
    group: ndarray or DataFrame
    exp_ldmat: list or dict
    simmat(gmat, gmap, meff, group, exp_ldmat, sub_id=sub_id, indwt=indwt, chrinterest=chrinterest,
           save=True, stdsim=True)
    """
    # Convert input data to DataFrame if necessary
    if not isinstance(gmap, pd.DataFrame):
        gmap = pd.DataFrame(gmap, copy=False)
    if not isinstance(meff, pd.DataFrame):
        meff = pd.DataFrame(meff, copy=False)
    if not isinstance(group, pd.DataFrame):
        group = pd.DataFrame(group, copy=False)
    if not isinstance(exp_ldmat, (list, dict)):
        raise ValueError("Expected within-family LD matrix should be a list or dict. Derive it using expldmat")
    sub_id = kwargs.get("sub_id", None)
    indwt = kwargs.get("indwt", None)
    chrinterest = kwargs.get("chrinterest", None)
    save = kwargs.get("save", False)
    stdsim = kwargs.get("stdsim", False)
    center = kwargs.get("center", False)
    progress = kwargs.get("progress", False)
    if sub_id is not None:
        sub_id = pd.DataFrame(sub_id)
    if sub_id is None or sub_id.shape[1] == 1:
        return simmat_g(gmat, gmap, meff, group, exp_ldmat, sub_id, indwt,
                        chrinterest, save, stdsim, center, progress)
    elif sub_id.shape[1] == 2:
        return simmat_z(gmat, gmap, meff, group, exp_ldmat, sub_id, indwt,
                        chrinterest, save, stdsim, center, progress)
    else:
        raise ValueError("sub_id must have either 1 or 2 columns.")

def simmat_g(gmat, gmap, meff, group, exp_ldmat, sub_id, indwt, chrinterest, save, stdsim,
             center, progress):
    """
    Compute similarity matrices using the gametic approach.
    Parameters:
        gmat (ndarray or DataFrame): Genotype data.
        gmap (ndarray or DataFrame): Genetic map information.
        meff (ndarray or DataFrame): Marker effects.
        group (ndarray or DataFrame): Grouping information.
        exp_ldmat (list or dict): Expected within-family LD matrix.
        sub_id (ndarray or DataFrame or None): Subset of individuals.
        indwt (ndarray or None): Index weights.
        chrinterest (list or None): List of chromosomes of interest.
        save (bool): Flag indicating whether to save matrices to files.
        stdsim (bool): Flag indicating whether to compute standardized similarity matrices.
        center (bool): Flag indicating whether to center marker effects.
        progress (bool): Flag indicating whether to display progress information.
    Returns:
        multgrpcov (list or dict): List or dictionary of similarity matrices computed based
        on groups and traits.
    """
    # Subset individuals if sub_id is provided
    if sub_id is not None:
        # Subset using sub_id
        subset_indices = subindcheck(group, sub_id)
        subset_indices2 = subset_indices * 2
        subset_indices2 = combine_arrays(subset_indices2, subset_indices2 + 1)
        if isinstance(gmat, pd.DataFrame):
            matsub = gmat.iloc[subset_indices2, :]
        else:
            matsub = gmat[subset_indices2, :]
        numbers = group.iloc[subset_indices, 1].astype(str).reset_index(drop=True)
        probn = pd.unique(numbers).tolist()
        alt_no = np.arange(0, len(probn), 1)
        noli = numbers.tolist()
        numbers = np.array([dict(zip(probn, alt_no)).get(x, x) for x in noli])
    else:
        # Use all individuals in group
        sub_id = pd.DataFrame(group.iloc[:, 0])
        subset_indices = subindcheck(group, sub_id)
        subset_indices2 = subset_indices * 2
        subset_indices2 = combine_arrays(subset_indices2, subset_indices2 + 1)
        if isinstance(gmat, pd.DataFrame):
            matsub = gmat.iloc[subset_indices2, :]
        else:
            matsub = gmat[subset_indices2, :]
        numbers = group.iloc[subset_indices, 1].astype(str).reset_index(drop=True)
        probn = pd.unique(numbers).tolist()
        alt_no = np.arange(0, len(probn), 1)
        noli = numbers.tolist()
        numbers = np.array([dict(zip(probn, alt_no)).get(x, x) for x in noli])
    # Check if the same map will be used for all groups
    if gmap.shape[1] - 3 == 1 and len(probn) > 1:
        print("The same map will be used for all groups")
    if gmap.shape[1] - 3 > 1 and (gmap.shape[1] - 3) != len(probn) > 1:
        raise ValueError(
            f"{gmap.shape[1] - 3} maps!= {len(probn)} groups; they must match for distinct group maps"
            )
    if gmap.shape[1] - 3 != len(exp_ldmat) > 1:
        raise ValueError("Number of maps does not match the number of groups in expected LD matrix")
    if meff.shape[1] > 1 and indwt is None:
        raise ValueError("Please provide index weights")
    if meff.shape[1] == 1 and indwt is not None:
        indwt = None
        print("There is only marker effect for one trait, index weight(s) is ignored")
    if indwt is not None:
        indwt = np.array(indwt)
        if meff.shape[1] != len(indwt):
            raise ValueError("Length of index weights differs from marker effects columns")
    chrinterest = chr_int(chrinterest)
    gmat1, alleles = formatgen(matsub, progress)
    mselist = addimsenumba(gmat1, np.array(meff.values, dtype=np.float32), np.array(alleles.tolist()), center)
    grp = gmap.shape[1] - 3
    trtnam = meff.columns
    multgrpcov = []
    unique_chromosomes = np.unique(gmap.iloc[:, 0].to_numpy())
    snp_indices = [np.where(gmap.iloc[:, 0].to_numpy() == chrm)[0] for chrm in unique_chromosomes]
    for gnp in range(grp):
        if grp > 1:
            tng = numbers == gnp
            us_ind = np.where(tng)[0]
            if progress:
                print("Processing group", probn[gnp])
        else:
            us_ind = np.arange(len(numbers))
        if progress:
            progr(0, meff.columns.size)
        for i in range(meff.shape[1]):
            cov_indxx = np.zeros((len(us_ind), len(us_ind)), dtype=np.float32)
            for chrm, s_ind in zip(unique_chromosomes, snp_indices):
                chrm_value = int(chrm)
                s_ind = np.array(s_ind)
                covtmpx = compute_covariance(mselist, exp_ldmat, us_ind, s_ind, i, chrm_value,
                                             gmap.shape[1] - 3, probn[gnp])
                cov_indxx += covtmpx
                write_matrices_to_file(chrinterest, stdsim, chrm_value, trtnam[i], covtmpx, save)
            if stdsim and save:
                covxfile = os.path.join(os.getcwd(), f"Stdsim mat for {meff.columns[i]}{' grp ' + str(probn[gnp]) if grp > 1 else ''}.npy")
                np.save(covxfile, cov2corr(cov_indxx))
            elif save:
                covxfile = os.path.join(os.getcwd(), f"Sim mat for {meff.columns[i]}{' grp ' + str(probn[gnp]) if grp > 1 else ''}.npy")
                np.save(covxfile, cov_indxx)
            if progress:
                progr(i + 1, meff.columns.size)
        snpindexxxx = np.arange(0, gmap.shape[0], 1)
        if meff.shape[1] == 1 and not stdsim:
            mat = cov_indxx
        elif meff.shape[1] == 1 and stdsim:
            mat = cov2corr(cov_indxx)
        elif meff.shape[1] > 1:
            print('Creating similarity matrix based on aggregate genotype')
            if progress:
                progr(0, max(pd.unique(gmap.iloc[:, 0])))
            tmpmt1 = np.empty((len(us_ind), meff.shape[0]), dtype=np.float32)
            for j, iii in enumerate(us_ind):
                tempagg = np.array([mselist[trt][iii, :] for trt in range(indwt.size)], dtype=np.float32)
                tmpmt1[j, :] = np.dot(indwt.T, tempagg)
            mat = np.zeros((len(us_ind), len(us_ind)), dtype=np.float32)
            for chrm in unique_chromosomes:
                chrm_value = int(chrm)
                s_ind = snpindexxxx[gmap.iloc[:, 0] == chrm_value]
                amat = tmpmt1[:, s_ind]
                ldmat = exp_ldmat[0][chrm_value - 1] if gmap.shape[1] - 3 == 1 else exp_ldmat[probn[gnp]][chrm_value - 1]
                temp1111 = np.dot(amat, ldmat)
                covtmpx = np.abs(np.dot(temp1111, amat.T))
                mat += covtmpx
                if progress:
                    progr(chrm_value, len(unique_chromosomes))
            if stdsim:
                mat = cov2corr(mat)
        multgrpcov.append(mat)
        if len(probn) == 1:
            break
    if grp > 1 and len(probn):
        multgrpcov = dict(zip(probn, multgrpcov))
    return multgrpcov

def simmat_z(gmat, gmap, meff, group, exp_ldmat, sub_id, indwt, chrinterest, save, stdsim,
             center, progress):
    """
    Compute similarity matrices using the zygotic approach.
    Parameters:
        gmat (ndarray or DataFrame): Genotype data.
        gmap (ndarray or DataFrame): Genetic map information.
        meff (ndarray or DataFrame): Marker effects.
        group (ndarray or DataFrame): Grouping information.
        exp_ldmat (list or dict): Expected within-family LD matrix.
        sub_id (ndarray or DataFrame or None): Subset of individuals.
        indwt (ndarray or None): Index weights.
        chrinterest (list or None): List of chromosomes of interest.
        save (bool): Flag indicating whether to save matrices to files.
        stdsim (bool): Flag indicating whether to compute standardized similarity matrices.
        center (bool): Flag indicating whether to center marker effects.
        progress (bool): Flag indicating whether to display progress information.
    Returns:
        mat (ndarray): Similarity matrix computed based on the zygotic approach.
    """
    # --- 1) Preliminary checks ---
    if gmap.shape[1] - 3 == 1 and len(pd.unique(group.iloc[:, 1])) > 1:
        print("The same map will be used for all groups")
    if gmap.shape[1] - 3 != len(exp_ldmat):
        raise ValueError("Number of maps does not match the number of groups in expected LD matrix")
    if meff.shape[1] > 1 and indwt is None:
        raise ValueError("Please provide index weights")
    if meff.shape[1] == 1 and indwt is not None:
        indwt = None
        print("There are only marker effects for one trait, index weight(s) is ignored")
    if indwt is not None:
        indwt = np.array(indwt)
        if meff.shape[1] != len(indwt):
            raise ValueError("Length of index weights differs from marker effects columns")
    # --- 2) Convert group labels in probn to integers (if float) ---
    mal1, fem1, _, probn = subindcheckzyg(group, sub_id)
    # If probn has float labels like 10009.0, 10131.0, convert them to int
    for i, val in enumerate(probn):
        if isinstance(val, float):
            probn[i] = int(val)
    chrinterest = chr_int(chrinterest)
    
    # --- 3) Build male/female indices ---
    male_indices = mal1 * 2
    male_indices = combine_arrays(male_indices, male_indices + 1)
    female_indices = fem1 * 2
    female_indices = combine_arrays(female_indices, female_indices + 1)
    print(f"Processing {probn[0]}")
    if isinstance(gmat, pd.DataFrame):
        gmat_male, alleles_male = formatgen(gmat.iloc[male_indices, :], progress)
    else:
        gmat_male, alleles_male = formatgen(gmat[male_indices, :], progress)
    print(f"Processing {probn[1]}")
    if isinstance(gmat, pd.DataFrame):
        gmat_female, alleles_female = formatgen(gmat.iloc[female_indices, :], progress)
    else:
        gmat_female, alleles_female = formatgen(gmat[female_indices, :], progress)
    mefff64 = np.array(meff.values, dtype=np.float32)
    mselist_male = addimsenumba(gmat_male, mefff64, np.array(alleles_male.tolist()), center)
    mselist_female = addimsenumba(gmat_female, mefff64, np.array(alleles_female.tolist()), center)
    print("Processing similarity matrices for zygotes")
    # --- 4) Loop over traits to build the similarity matrix ---
    for i in range(meff.shape[1]):
        mat = np.zeros((len(mal1), len(mal1)), dtype=np.float32)
        # Convert the chromosome column to int to avoid float indexing
        chrom_col = gmap.iloc[:, 0].astype(int).to_numpy()
        # Loop over unique chromosomes as int
        for chrm_int in np.unique(chrom_col):
            s_ind = np.where(chrom_col == chrm_int)[0]
            # If there's only 1 map-dist column => use exp_ldmat[0][...]
            if gmap.shape[1] - 3 == 1:
                # Make sure chrm_int - 1 is valid
                covtmpx = (np.abs(mrmmult(mselist_male[i][:, s_ind], exp_ldmat[0][chrm_int - 1])) +
                           np.abs(mrmmult(mselist_female[i][:, s_ind], exp_ldmat[0][chrm_int - 1])))
            else:
                # multi-group => use exp_ldmat[probn[0] or probn[1]]
                # Ensure probn[0], probn[1] are also int or valid keys
                covtmpx = (
                    np.abs(mrmmult(mselist_male[i][:, s_ind], exp_ldmat[probn[0]][chrm_int - 1])) +
                    np.abs(mrmmult(mselist_female[i][:, s_ind], exp_ldmat[probn[1]][chrm_int - 1]))
                )
            mat += covtmpx
            # Pass chrm_int (an int) to the file-writing function
            write_matrices_to_file(chrinterest, stdsim, chrm_int, meff.columns[i], covtmpx, save)
        # Save the full matrix if needed
        if stdsim:
            if save:
                covxfile = os.path.join(os.getcwd(), f"Sim mat zygote {meff.columns[i]}.npy")
                np.save(covxfile, cov2corr(mat))
        else:
            if save:
                covxfile = os.path.join(os.getcwd(), f"Sim mat zygotes {meff.columns[i]}.npy")
                np.save(covxfile, mat)
        if progress:
            progr(i + 1, meff.shape[1])
    # -----------------------------------------------------------------------
    # 5) If more than one trait, build an aggregate similarity matrix
    # -----------------------------------------------------------------------
    snpindexxx = np.arange(0, gmap.shape[0], 1)
    if meff.shape[1] == 1:
        if stdsim:
            mat = cov2corr(mat)
    else:
        # multi-trait
        if progress:
            print('Creating similarity matrix based on aggregate genotype')
            # must cast to int for the progress bar
            max_chr = int(np.max(chrom_col))
            progr(0, max_chr)
        tmpmm = np.empty((sub_id.shape[0], gmap.shape[0]), dtype=np.float32)
        tmpmfm = np.empty((sub_id.shape[0], gmap.shape[0]), dtype=np.float32)
        for i in range(sub_id.shape[0]):
            tmpmt1 = np.zeros((indwt.size, gmap.shape[0]), dtype=np.float32)
            tmpmt2 = np.zeros((indwt.size, gmap.shape[0]), dtype=np.float32)
            for trt in range(indwt.size):
                tmpmt1[trt, :] = mselist_male[trt][i, :]
                tmpmt2[trt, :] = mselist_female[trt][i, :]
            tmpmm[i, :] = np.matmul(indwt.T, tmpmt1)
            tmpmfm[i, :] = np.matmul(indwt.T, tmpmt2)
        mat = np.zeros((sub_id.shape[0], sub_id.shape[0]), dtype=np.float32)
        for chrm_int in np.unique(chrom_col):
            s_ind = snpindexxx[chrom_col == chrm_int]
            # If single map-dist col => use exp_ldmat[0][chrm_int - 1]
            amm = exp_ldmat[0][chrm_int - 1] if (gmap.shape[1] - 3 == 1) else exp_ldmat[probn[0]][chrm_int - 1]
            amf = exp_ldmat[0][chrm_int - 1] if (gmap.shape[1] - 3 == 1) else exp_ldmat[probn[1]][chrm_int - 1]
            amatm = tmpmm[:, s_ind]
            tmp1111 = np.matmul(amatm, amm)
            mat += np.abs(np.matmul(tmp1111, amatm.T))
            amatf = tmpmfm[:, s_ind]
            tmp1111 = np.matmul(amatf, amf)
            mat += np.abs(np.matmul(tmp1111, amatf.T))
            if progress:
                progr(chrm_int, max_chr)
        if stdsim:
            mat = cov2corr(mat)
    return mat

@njit
def calculate_par_spec_me(gmat, alleles, haplotype):
    """
    Recodes haplotypes for parent-specific marker effects.
    Args:
        gmat (np.ndarray): Genetic marker matrix.
        alleles (list): List of alleles.
        haplotype (bool): Flag indicating whether to gmat is haplotype data.
    Returns:
        np.ndarray: recoded haplotpes matrix.
    """
    # Sum up rows and columns
    if haplotype:
        par_spec_me = (gmat[::2] * 10) + gmat[1::2]
        # Determine the number of alleles and ref allele
        no_alleles = len(alleles)
        ref_allele = alleles[0]
        if no_alleles == 2:
            # Two alleles: ref allele and alt allele
            alt_allele = alleles[1]
            # Calculate the genotype combinations
            homo_dom = ref_allele * 10 + ref_allele
            homo_rec = alt_allele * 10 + alt_allele
            het_dom = ref_allele * 10 + alt_allele
            het_rec = alt_allele * 10 + ref_allele
            # Assign values to the par_spec_me array based on genotype combinations
            par_spec_me = np.where(par_spec_me == homo_rec, -1, par_spec_me)
            par_spec_me = np.where((par_spec_me == het_dom) |
                                   (par_spec_me == het_rec), 0, par_spec_me)
            par_spec_me = np.where(par_spec_me == homo_dom, 1, par_spec_me)
        else:
            # Three alleles: ref allele, alt allele, and missing allele
            alt_allele = alleles[1]
            missing_allele = alleles[2]
            # Calculate the genotype combinations
            homo_dom = ref_allele * 10 + ref_allele
            homo_rec = alt_allele * 10 + alt_allele
            het_dom = ref_allele * 10 + alt_allele
            het_rec = alt_allele * 10 + ref_allele
            mis_maj = missing_allele * 10 + ref_allele
            mis_min = missing_allele * 10 + alt_allele
            maj_mis = ref_allele * 10 + missing_allele
            min_mis = alt_allele * 10 + missing_allele
            # Assign values to the par_spec_me array based on genotype combinations
            for i in range(par_spec_me.shape[0]):
                for j in range(par_spec_me.shape[1]):
                    if par_spec_me[i, j] == homo_rec:
                        par_spec_me[i, j] = -1
                    elif par_spec_me[i, j] == het_dom or par_spec_me[i, j] == het_rec:
                        par_spec_me[i, j] = 0
                    elif par_spec_me[i, j] == homo_dom:
                        par_spec_me[i, j] = 1
                    elif par_spec_me[i, j] == mis_maj or par_spec_me[i, j] == mis_min or \
                            par_spec_me[i, j] == maj_mis or par_spec_me[i, j] == min_mis:
                        par_spec_me[i, j] = 0
    else:
        par_spec_me = gmat - 1
    
    return par_spec_me

@njit(fastmath=True)
def derivebv(par_spec_me, meff):
    """derives breeding values"""
    return par_spec_me @ meff

def calcbvnumba(gmat, meff, alleles, haplotype, center):
    """
    Computes parent-specific additive marker effects for each trait.
    Args:
        gmat (np.ndarray): Haplotype matrix.
        meff (np.ndarray): Marker effect matrix.
        alleles (list): List of alleles.
        haplotype (bool): Flag indicating whether to gmat is haplotype data.
        center (bool): Flag indicating whether to center the trait matrix.
    Returns:
        list: List of IMSE matrices for each trait.
    """
    par_spec_me = calculate_par_spec_me(gmat, alleles, haplotype)
    
    # column centering
    if center:
        mean_trait = np.mean(par_spec_me, axis=0)
        par_spec_me = par_spec_me.astype(np.float32)  # Convert to float32
        par_spec_me -= mean_trait
    
    if not center:
        par_spec_me = par_spec_me.astype(np.float32)
    
    breedval = derivebv(par_spec_me, meff)
    return breedval

def calcgbv(gmat, meff, group, indwt, sub_id, haplotype, center, progress):
    """
    Calculate breeding values for each trait.
    Parameters:
    -----------
    gmat : pandas.DataFrame
        Genotype matrix.
    meff : numpy.ndarray or pandas.DataFrame
        Marker effect matrix.
    group : pandas.DataFrame
        Group information.
    indwt : numpy.ndarray or None
        Index weights for marker effects. If None, index weights are ignored.
    sub_id : pandas.DataFrame or None, optional
        Subset IDs of individuals. If None, all individuals are considered. Default is None.
    haplotype : bool
        Whether gmat is haplotypes. Defaults to true
    center : bool
        Flag indicating whether to center the breeding values. 
    progress : bool
        Flag indicating whether to display the progress.
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing breeding values for each trait.
    Notes:
    ------
    If sub_id is not None, breeding values will be calculated only for the individuals 
    specified in sub_id.
    The calculated breeding values are stored in a DataFrame with columns for each trait 
    and additional columns
    for ID, Group, and ABV (if multiple traits exist).
    """
    meff = pd.DataFrame(meff)
    if sub_id is not None:
        sub_id = pd.DataFrame(sub_id)
        aaa = subindcheck(group, sub_id)
        idn = group.iloc[aaa, 0].astype(str).reset_index(drop=True)  # ID
        groupsex = group.iloc[aaa, 1].astype(str).reset_index(drop=True)
        aaa = aaa * 2
        aaa = combine_arrays(aaa, aaa + 1)
        matsub = gmat.iloc[aaa, :]
    else:
        idn = group.iloc[:, 0].astype(str).reset_index(drop=True)  # ID
        groupsex = group.iloc[:, 1].astype(str).reset_index(drop=True)
        matsub = gmat
    # Checks to ensure there is index weight
    if meff.shape[1] > 1 and indwt is None:
        raise ValueError("Please provide index weights")
    # ignore index weight of single trait
    if meff.shape[1] == 1 and indwt is not None:
        indwt = None
        print("There is only marker effect for one trait, index weights are ignored")
    if indwt is not None:
        indwt = np.array(indwt)
        if meff.shape[1] != len(indwt):
            raise ValueError("Please provide index weights")
    trait_names = meff.columns  # traits names
    gmat1, alleles = formatgen(matsub, progress)
    gbv = calcbvnumba(gmat1, np.array(meff.values, dtype=np.float32), np.array(alleles.tolist()),
                      haplotype, center)
    gbv = pd.DataFrame(gbv)
    gbv.columns = trait_names
    if len(trait_names) > 1:
        abv = gbv @ indwt
        gbv.insert(len(trait_names), "ABV", abv, True)
    gbv.insert(0, "ID", idn, True)  # insert ID
    gbv.insert(1, "Group", groupsex, True)  # insert group
    return gbv

def selint(propsel):
    """
    Compute the selection intensity.
    Parameters
    ----------
    propsel : float
        Proportion value.
    Returns
    -------
    sel_int : float
        Selection intensity.
    """
    return norm.pdf(norm.ppf(1 - propsel)) / propsel

# def stdtruncp(propsel):
#     """
#     Calculate the standard truncated proportion.
#     Parameters
#     ----------
#     propsel : float
#         Proportion value.
#     Returns
#     -------
#     std_trunc_p : float
#         Standard truncated proportion.
#     """
#     if propsel == 1:
#         propsel = 0.999999999
#     elif propsel == 0:
#         propsel = 1 - 0.999999999
#     prob = propsel
#     appp = 1 - prob if prob > 0.5 else prob
#     ttt = math.sqrt(math.log(1 / (appp * appp)))
#     xnorx = ttt - (2.530517 + ttt * (0.802853 + ttt * 0.010328)) / (
#             1 + ttt * (1.432788 + ttt * (0.189269 + ttt * 0.001308)))
#     if appp != prob:
#         xnorx = -xnorx
#     return round(xnorx, 4)

def stdtruncp(propsel):
    """
    Calculate the selection intensity for a standard normal distribution,
    given the proportion propsel of individuals selected from the top tail.
    The formula is i = dnorm(qnorm(1 - propsel)) / propsel.
    Parameters
    ----------
    propsel : float
        Proportion selected (0 < propsel < 1).

    Returns
    -------
    float
        The selection intensity, rounded to 4 decimal places.
    """
    # Handle edge cases to avoid inf/nan
    if propsel <= 0:
        propsel = 1e-9
    elif propsel >= 1:
        propsel = 1 - 1e-9
    # z = inverse CDF (qnorm in R) for the upper tail of size propsel
    z = norm.ppf(1 - propsel)
    # PDF at z
    pdfz = norm.pdf(z)
    # selection intensity i = PDF / p
    i_val = pdfz / propsel
    return round(i_val, 4)

def calcindex(gmat, meff, group, indwt, sub_id, msvmsc, criterion, prop_sel, aggregate,
    haplotype, center, progress):
    """
    Calculate selection criteria (index or UC) for matepairs using a gametic approach.
    Parameters:
    -----------
    gmat : pd.DataFrame
        Genotype matrix.
    meff : pd.DataFrame
        Marker effects.
    group : pd.DataFrame
        Group information.
    indwt : ndarray or None
        Index weights for marker effects. Required for multiple marker effects.
    sub_id : ndarray or None
        Sub-individual IDs.
    msvmsc : pd.DataFrame
        Multi-trait selection variance-covariance matrix.
    criterion : str
        Selection strategy ('index' or 'UC').
    prop_sel : float
        Proportion of selection.
    aggregate : bool, optional
        Whether to use aggregate genotype or weighted (index) sum of trait-specific msvs.
        Defaults to True.
    haplotype : bool
        Whether gmat is haplotypes. Defaults to true
    center : bool
        Whether the to center the genotype matrix.
    progress : bool
        Whether to show the progress bar.
    Returns:
    --------
    pd.DataFrame
        Selection criteria values for matepairs.
    """
    meff = pd.DataFrame(meff)
    if sub_id is not None:
        aaa = subindcheck(msvmsc.iloc[:, 0:2], sub_id)
    else:
        aaa = subindcheck(msvmsc.iloc[:, 0:2], msvmsc.iloc[:, 0])
        sub_id = msvmsc.iloc[:, 0]
    trait_names = meff.columns  # traits names
    gbv = calcgbv(gmat, meff, group, indwt, sub_id, haplotype, center, progress)
    no_individuals = gbv.shape[0]  # Number of individuals
    notr = trait_names.size
    indexdf = np.zeros((no_individuals, notr + 1))
    if notr == 1:
        if criterion.lower() == "index":
            indexdf[:, 0] = gbv.iloc[:, 2] + (np.sqrt(2) * stdtruncp(prop_sel) * np.sqrt(
                msvmsc.iloc[aaa, 2])).reset_index(drop=True)
        elif criterion.lower() == "uc":
            indexdf[:, 0] = gbv.iloc[:, 2] + (selint(prop_sel) * np.sqrt(
                msvmsc.iloc[aaa, 2])).reset_index(drop=True)
    else:
        sub_msvmsc = msvmsc[trait_names]
        for i in range(notr):
            if criterion.lower() == "index":
                indexdf[:, i] = gbv.iloc[:, (i + 2)] + (np.sqrt(2) * stdtruncp(prop_sel) * np.sqrt(
                    sub_msvmsc.iloc[aaa, i])).reset_index(drop=True)
            elif criterion.lower() == "uc":
                indexdf[:, i] = gbv.iloc[:, (i + 2)] + (selint(prop_sel) * np.sqrt(
                    sub_msvmsc.iloc[aaa, i])).reset_index(drop=True)
        if criterion.lower() == "index":
            if aggregate:
                agg = msvmsc["AG"]
                indexdf[:, notr] = gbv.iloc[:, (notr + 2)] + (np.sqrt(2) * stdtruncp(
                    prop_sel) * np.sqrt(agg[aaa])).reset_index(drop=True)
            else:
                indexdf[:, notr] = gbv.iloc[:, (notr + 2)] + (np.sqrt(2) * stdtruncp(
                    prop_sel) * np.sqrt(sub_msvmsc.iloc[aaa, :].sum(axis=1))).reset_index(drop=True)
        elif criterion.lower() == "uc":
            if aggregate:
                agg = msvmsc["AG"]
                indexdf[:, notr] = gbv.iloc[:, (notr + 2)] + (selint(prop_sel) * np.sqrt(
                    agg[aaa])).reset_index(drop=True)
            else:
                indexdf[:, notr] = gbv.iloc[:, (notr + 2)] + (selint(prop_sel) * np.sqrt(
                    sub_msvmsc.iloc[aaa, :].sum(axis=1))).reset_index(drop=True)
    indexdf = pd.DataFrame(indexdf)
    colnames = np.concatenate((trait_names, ["ABV"]), axis=None)
    indexdf.columns = colnames
    indexdf = pd.concat([gbv.iloc[:, 0:2], indexdf], axis=1)
    return indexdf

def calcgbvzyg(gmat, meff, group, indwt, sub_id, haplotype, center, progress):
    """
    Calculate breeding values for matepairs using a gametic approach.
    Parameters:
    -----------
    gmat : pd.DataFrame
        Genotype matrix.
    meff : pd.DataFrame
        Marker effects.
    group : pd.DataFrame
        Group information.
    indwt : ndarray or None
        Index weights for marker effects.
    sub_id : ndarray or None
        Sub-individual IDs.
   haplotype : bool
       Whether gmat is haplotypes. Defaults to true
    center : bool
        Whether to center the genotype matrix.
    progress : bool
        Whether to show the progress bar.
    Returns:
    --------
    pd.DataFrame
        Breeding values for matepairs.
    """
    mal1, fem1, idn, gid = subindcheckzyg(group, sub_id)
    sub_idxxxx = None
    data = calcgbv(gmat, meff, group, indwt, sub_idxxxx, haplotype, center, progress)
    gbv1 = data.iloc[mal1, 2:].reset_index(drop=True)
    gbv2 = data.iloc[fem1, 2:].reset_index(drop=True)
    result = (gbv1 + gbv2) / 2
    result.insert(0, gid[1] + "_ID", idn[:, 1], True)  # insert FemaleID column
    result.insert(0, gid[0] + "_ID", idn[:, 0], True)  # insert MaleID column
    return result

def calcindexzyg(gmat, meff, group, indwt, sub_id, msvmsc, criterion, prop_sel, aggregate,
                 haplotype, center, progress):
    """
    Calculate selection criteria (index or UC) for matepairs using a gametic approach.
    Parameters:
    -----------
    gmat : pd.DataFrame
        Genotype matrix.
    meff : pd.DataFrame
        Marker effects.
    group : pd.DataFrame
        Group information.
    indwt : ndarray or None
        Index weights for marker effects. Required for multiple marker effects.
    sub_id : ndarray or None
        Sub-individual IDs.
    msvmsc : pd.DataFrame
        Multi-trait selection variance-covariance matrix.
    criterion : str
        Selection strategy ('index' or 'UC').
    prop_sel : float
        Proportion of selection.
    aggregate : bool, optional
        Whether to use aggregate genotype or weighted (index) sum of trait-specific msvs. 
        Defaults to True.
    haplotype : bool
        Whether gmat is haplotypes. Defaults to true
    center : bool
        Whether to center the genotype matrix.
    progress : bool
        Whether to show the progress bar.
    Returns:
    --------
    pd.DataFrame
        Selection criteria values for matepairs.
    """
    mal1, fem1, idn, gid = subindcheckzyg(msvmsc.iloc[:, 0:2], sub_id)
    sub_idxxxx = None
    trait_names = meff.columns  # traits names
    gbv = calcgbv(gmat, meff, group, indwt, sub_idxxxx, haplotype, center, progress)
    gbv1 = gbv.iloc[mal1, 2:].reset_index(drop=True)
    gbv2 = gbv.iloc[fem1, 2:].reset_index(drop=True)
    result = (gbv1 + gbv2) / 2
    msv1 = msvmsc.iloc[mal1, 2:].reset_index(drop=True)
    msv2 = msvmsc.iloc[fem1, 2:].reset_index(drop=True)
    msvdf = msv1 + msv2
    no_individuals = result.shape[0]  # Number of individuals
    notr = trait_names.size
    if notr == 1:
        indexdf = np.zeros((no_individuals, notr))
        if criterion in {"index", "INDEX", "Index"}:
            indexdf[:, 0] = result.iloc[:, 0] + (np.sqrt(2) * stdtruncp(prop_sel) * np.sqrt(
                msvdf.iloc[:, 0]))
        elif criterion in {"UC", "uc", "Uc"}:
            indexdf[:, 0] = result.iloc[:, 0] + (selint(prop_sel) * np.sqrt(msvdf.iloc[:, 0]))
        indexdf = pd.DataFrame(indexdf)
        indexdf.columns = trait_names
    elif notr > 1:
        indexdf = np.zeros((no_individuals, notr + 1))
        sub_msvmsc = msvdf[trait_names]
        for i in range(notr):
            if criterion in {"index", "INDEX", "Index"}:
                indexdf[:, i] = result.iloc[:, i] + (
                        np.sqrt(2) * stdtruncp(prop_sel) * np.sqrt(sub_msvmsc.iloc[:, i]))
            elif criterion in {"UC", "uc", "Uc"}:
                indexdf[:, i] = result.iloc[:, i] + (selint(prop_sel) * np.sqrt(
                    sub_msvmsc.iloc[:, i]))
        if criterion in {"index", "INDEX", "Index"}:
            if aggregate:
                indexdf[:, notr] = result.iloc[:, notr] + (
                        np.sqrt(2) * stdtruncp(prop_sel) * np.sqrt(msvdf["AG"]))
            else:
                indexdf[:, notr] = result.iloc[:, notr] + (
                        np.sqrt(2) * stdtruncp(prop_sel) * np.sqrt(sub_msvmsc.sum(axis=1)))
        elif criterion in {"UC", "uc", "Uc"}:
            if aggregate:
                indexdf[:, notr] = result.iloc[:, notr] + (selint(prop_sel) * np.sqrt(msvdf["AG"]))
            else:
                indexdf[:, notr] = result.iloc[:, notr] + (selint(prop_sel) * np.sqrt(
                    sub_msvmsc.sum(axis=1)))
        indexdf = pd.DataFrame(indexdf)
        colnames = np.concatenate((trait_names, "ABV"), axis=None)
        indexdf.columns = colnames
    indexdf.insert(0, gid[1] + "_ID", idn[:, 1], True)  # insert FemaleID column
    indexdf.insert(0, gid[0] + "_ID", idn[:, 0], True)  # insert MaleID column
    return indexdf

def selstrat(gmat, meff, group, **kwargs):
    """
    Calculate selection criteria (GEBV, UC, or index) for gametes or zygotes.
    Parameters:
    -----------
    gmat : pd.DataFrame
        Genotype matrix.
    meff : pd.DataFrame
        Marker effects.
    group : pd.DataFrame
        Group information.
    **kwargs : Optional keyword arguments for customization.
        indwt : ndarray or None
            Index weights for marker effects. Required for multiple marker effects. 
            If None, index weights are ignored.
        sub_id : ndarray or None
            Subset IDs of individuals. If None, all individuals are considered.
        msvmsc : pd.DataFrame or None
            Multi-trait selection variance-covariance matrix. Required for UC or index calculations.
        criterion : str, optional
            Selection strategy ('gebv', 'uc', or 'index'). Defaults to 'gebv'.
        prop_sel : float, optional
            Proportion of selection. Required for UC or index calculations. 
            Should be in the range [0, 1].
        aggregate : bool, optional
            Whether to use aggregate genotype or weighted (index) sum of trait-specific msvs. 
            Defaults to True.
        haplotype : bool
            Whether gmat is haplotypes. Defaults to True
        center : bool, optional
            Whether to center the genotype matrix. Defaults to True.
        progress : bool, optional
            Whether to show the progress bar. Defaults to False.
    Returns:
    --------
    pd.DataFrame
        Selection criteria values.
    Raises:
    -------
    ValueError
        If the number of index weights does not match the number of traits (marker effects).
        If the selection strategy is invalid.
        If msvmsc is not provided or prop_sel is not in the range [0, 1] for UC or index calcs.
    """
    # Convert input data to DataFrame if necessary
    meff = pd.DataFrame(meff)
    group = pd.DataFrame(group)
    # Get optional arguments
    indwt = kwargs.get("indwt", None)
    sub_id = kwargs.get("sub_id", None)
    msvmsc = kwargs.get("msvmsc", None)
    criterion = kwargs.get("criterion", "gebv").lower()
    prop_sel = kwargs.get("prop_sel")
    aggregate = kwargs.get("aggregate", True)
    haplotype = kwargs.get("haplotype", True)
    center = kwargs.get("center", True)
    progress = kwargs.get("progress", False)
    # Check if index weights are provided when multiple marker effects are present
    if meff.shape[1] > 1 and indwt is None:
        raise ValueError("Please provide index weights")
    # Ignore index weights if there is only one marker effect
    if meff.shape[1] == 1 and indwt is not None:
        indwt = None
        print("There is only marker effect for one trait, index weights are ignored")
    # Check if the number of index weights matches the number of marker effects
    if indwt is not None:
        indwt = np.array(indwt)
        if meff.shape[1] != len(indwt):
            raise ValueError("length of index weights differs from marker effects columns")
    # Check valid selection strategy
    valid_strategies = ['gebv', 'uc', 'index']
    if criterion not in valid_strategies:
        raise ValueError(f"Invalid selection strategy. Expected one of {valid_strategies}.")
    # Check msvmsc and prop_sel for 'uc' or 'index' strategy
    if criterion in ['uc', 'index'] and (msvmsc is None or not 0 <= prop_sel <= 1):
        raise ValueError("For UC or index calculation, 'msvmsc' dataframe should be provided,"
                         " and 'prop_sel' should be in the range [0, 1].")
    if sub_id is None or pd.DataFrame(sub_id).shape[1] == 1:
        if criterion == 'gebv':
            data = calcgbv(gmat, meff, group, indwt, sub_id, haplotype, center, progress)
        elif criterion in ['uc', 'index']:
            data = calcindex(gmat, meff, group, indwt, sub_id, msvmsc, criterion, prop_sel,
                             aggregate, haplotype, center, progress)
    elif pd.DataFrame(sub_id).shape[1] == 2:
        if criterion == 'gebv':
            data = calcgbvzyg(gmat, meff, group, indwt, sub_id, haplotype, center, progress)
        elif criterion in ['uc', 'index']:
            data = calcindexzyg(gmat, meff, group, indwt, sub_id, msvmsc, criterion, prop_sel,
                                aggregate, haplotype, center, progress)
    return data
