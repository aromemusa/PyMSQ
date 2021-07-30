"""Genetic evaluation of individuals."""
import os
import sys
# import time
from collections import Counter
from itertools import compress
from numba import njit
import pkg_resources
import numpy as np
import pandas as pd
import scipy.linalg
import scipy.stats


def example_data():
    """Provide data to the package."""
    cwd = os.getcwd()
    stream = pkg_resources.resource_stream(__name__, 'data/chr.txt')
    chrmosomedata = pd.read_table(stream, sep=" ")
    stream = pkg_resources.resource_stream(__name__, 'data/group.txt')
    groupdata = pd.read_table(stream, sep=" ")
    stream = pkg_resources.resource_stream(__name__, 'data/effects.txt')
    markereffdata = pd.read_table(stream, sep=" ")
    stream = pkg_resources.resource_stream(__name__, 'data/phase.txt')
    genodata = pd.read_table(stream, header=None, sep=" ")
    stream = pkg_resources.resource_stream(__name__, 'data/ped.txt')
    ped = pd.read_table(stream, header=None, sep=" ")
    os.chdir(cwd)
    return chrmosomedata, markereffdata, genodata, groupdata, ped


if __name__ == "__main__":
    example_data()


@njit
def fnrep2(gen, aaxx, aaxx1):
    """Code phased genotypes into 1, 2, 3 and 4."""
    qqq = np.empty((int(gen.shape[0]/2), gen.shape[1]), np.int_)
    for i in range(qqq.shape[0]):
        for j in range(qqq.shape[1]):
            if gen[2*i, j] == aaxx and gen[2*i+1, j] == aaxx:
                qqq[i, j] = 1
            elif gen[2*i, j] == aaxx1 and gen[2*i+1, j] == aaxx1:
                qqq[i, j] = 2
            elif gen[2*i, j] == aaxx and gen[2*i+1, j] == aaxx1:
                qqq[i, j] = 3
            else:
                qqq[i, j] = 4
    return qqq


def haptogen(gen, progress=False):
    """Convert haplotypes to coded genotypes."""
    if progress:
        print("Converting phased haplotypes to genotypes")
    if gen.shape[1] == 2:
        gen = np.array(gen.iloc[:, 1])  # del col containing ID
        # convert string to 2D array of integers
        gen = [list(gen[i].rstrip()) for i in range(gen.shape[0])]
        gen = np.array(gen, int)
        # derives the frequency of alleles to determine the major allele
        allele = np.asarray(np.unique(gen, return_counts=True)).T.astype(int)
        if len(allele[:, 0]) != 2:
            sys.exit("method only supports biallelic markers")
        aaxx = allele[:, 0][np.argmax(allele[:, 1])]  # major allele
        aaasns = np.isin(allele[:, 0], aaxx, invert=True)
        aaxx1 = int(allele[:, 0][aaasns])  # minor allele
        gen = np.array(gen, int)
        gen = fnrep2(gen, aaxx, aaxx1)
    elif gen.shape[1] > 2:
        gen = gen.iloc[:, 1:gen.shape[1]]  # del col containing ID
        # derives the frequency of alleles to determine the major allele
        allele = np.asarray(np.unique(gen, return_counts=True)).T.astype(int)
        if len(allele[:, 0]) != 2:
            sys.exit("method only supports biallelic markers")
        aaxx = allele[:, 0][np.argmax(allele[:, 1])]  # major allele
        aaasns = np.isin(allele[:, 0], aaxx, invert=True)
        aaxx1 = int(allele[:, 0][aaasns])  # minor allele
        gen = np.array(gen, int)
        gen = fnrep2(gen, aaxx, aaxx1)
    return gen


class Datacheck:
    """Check the input data for errors and store relevant info as an object."""

    def __init__(self, gmap, meff, gmat, group, indwt, progress=False):
        """
        Check input data for errors and store relevant info as class object.

        Parameters
        ----------
        gmap : pandas.DataFrame
            Index: RangeIndex
            Columns:
            Name: CHR, dtype: int64; chromosome number
            Name: SNPName, dtype: object; marker name
            Name: Position: dtype: int64; marker position in bp
            Name: group: dtype: float64; marker distance (cM) or reco rates
        meff : pandas.DataFrame
            Index: RangeIndex
            Columns:
            Name: trait names: float64; no. of columns = no of traits
        gmat : pandas.DataFrame
            Index: RangeIndex
            Columns:
            Name: ID, dtype: int64 or str; identification of individuals
            Name: haplotypes, dtype: object; must be biallelic
        group : pandas.DataFrame
            Index: RangeIndex
            Columns:
            Name: group, dtype: object; group code of individuals, e.g., M, F
            Name: ID, dtype: int64 or str; identification of individuals
        indwt : list of index weights for each trait
        progress : bool, optional; print progress of the function if True
        Returns stored input files
        -------
        """
        # check: ensures number of traits match size of index weights
        indwt = np.array(indwt)
        if (meff.shape[1]-1) != indwt.size:
            sys.exit('no. of index weights do not match marker effects cols')
        # check: ensure individuals' genotypes match group and ID info
        id_indgrp = pd.Series(group.iloc[:, 1]).astype(str)  # no of inds
        if not pd.Series(
                pd.unique(gmat.iloc[:, 0])).astype(str).equals(id_indgrp):
            sys.exit("ID of individuals in group & genotypic data don't match")
        # check: ensure marker names in marker map and effects match
        if not (gmap.iloc[:, 1].astype(str)).equals(meff.iloc[:, 0].astype(str)):
            print("Discrepancies between marker names")
            sys.exit("Check genetic map and marker effects")
        # check: ensure marker or allele sub effect are all numeric
        meff = meff.iloc[:, 1:meff.shape[1]]
        test = meff.apply(
            lambda s: pd.to_numeric(s, errors='coerce').notnull().all())
        if not test.all():
            sys.exit("Marker or allele sub effects contain non-numeric values")
        # check: ensure unique maps match no of groups if map more than 1
        grpg = pd.unique(group.iloc[:, 0])  # groups of individuals
        grp_chrom = gmap.shape[1]-3  # no of unique maps
        gmat = haptogen(gmat, progress)
        if grp_chrom > 1 and grp_chrom != grpg.size:
            sys.exit("no. of unique maps does not match no. of groups")
        # check no of markers in genotype and map and marker effects match
        no_markers = gmap.shape[0]  # no of markers
        if no_markers != gmat.shape[1] or no_markers != meff.shape[0]:
            sys.exit("markers nos in gen, chrm or marker effects don't match")
        # check: ordered marker distance or recombination rates
        for grn in range(grp_chrom):
            for chrm in pd.unique(gmap.iloc[:, 0]):
                mpx = np.array(gmap.iloc[:, 3+grn][gmap.iloc[:, 0] == chrm])
                if not (mpx == np.sort(sorted(mpx))).any():
                    sys.exit(
                        f"Faulty marker map on chr {chrm} for grp {grpg[grn]}")
        if progress:
            print('Data passed the test!')
            print("Number of individuals:  ", len(id_indgrp))
            print("Number of groups:       ", len(grpg), ": ", grpg)
            print("Number of specific maps:", grp_chrom)
            print("Number of chromosomes:  ", len(pd.unique(gmap.iloc[:, 0])))
            print("Total no. markers:      ", no_markers)
            print("Number of trait(s):     ", meff.columns.size)
            print("Trait name(s) and Index weight(s)")
            if meff.columns.size == 1:
                print(meff.columns[0], ": ", indwt[0])
            elif meff.columns.size > 1:
                for i in range(meff.columns.size):
                    print(meff.columns[i], ": ", indwt[i])
        self.gmap = gmap
        self.meff = meff
        self.gmat = gmat
        self.group = group
        self.indwt = indwt


def elem_cor(mylist, mprc, ngp, mposunit, method, chrm):
    """Derive pop cov matrix."""
    if method == 1:  # Bonk et al's approach
        if mposunit in ("cM", "cm", "CM", "Cm"):
            tmp = np.exp(-2*(np.abs(mprc - mprc[:, None])/100))/4
        elif mposunit in ("reco", "RECO"):
            if mprc[0] != 0:
                sys.exit(f"First value for reco rate on chr {chrm} isn't zero")
            aaa = (1-(2*mprc))/4
            ida = np.arange(aaa.size)
            tmp = aaa[np.abs(ida - ida[:, None])]
    elif method == 2:  # Santos et al's approach
        if mposunit in ("cM", "cm", "CM", "Cm"):
            tmp = (-1*(np.abs(mprc - mprc[:, None])/200))+0.25
            cutoff = (-1*(50/200))+0.25
            tmp = np.where(tmp < cutoff, 0, tmp)
        elif mposunit in ("reco", "RECO"):
            if mprc[0] != 0:
                sys.exit(f"First value for reco rate on chr {chrm} isn't zero")
            aaa = (-1*(mprc/2))+0.25
            ida = np.arange(aaa.size)
            tmp = aaa[np.abs(ida - ida[:, None])]
            cutoff = (-1*(0.5/2))+0.25
            tmp = np.where(tmp < cutoff, 0, tmp)
    # append chromosome-specific covariance matrix to list
    mylist[int(ngp)].append(tmp)
    return mylist


def popcovmat(info, mposunit, method):
    """
    Derive population-specific covariance matrices.

    Parameters
    ----------
    info : class object
        A class object created using the function "datacheck"
    mposunit : string
        A sting with containing "cM" or "reco".
    method : int
        An integer with a value of 1 for Bonk et al.'s approach or
        2 for Santos et al's approach'
    Returns
    -------
    mylist : list
        A list containing group-specific pop covariance matrices for each chr.
    """
    if mposunit not in ("cM", "cm", "CM", "Cm", "reco", "RECO"):
        sys.exit("marker unit should be either cM or reco")
    # unique group name for naming the list if map is more than 1
    probn = pd.unique(info.group.iloc[:, 0].astype(str)).tolist()
    chromos = pd.unique(info.gmap.iloc[:, 0])  # chromosomes
    no_grp = info.gmap.shape[1]-3  # no of maps
    mylist = []  # list stores chromosome-wise covariance matrix
    for ngp in range(no_grp):
        mylist.append([])
        # marker position in cM or recombination rates
        grouprecodist = info.gmap.iloc[:, 3+ngp]
        for chrm in chromos:
            mpo = np.array(grouprecodist[info.gmap.iloc[:, 0] == (chrm)])
            elem_cor(mylist, mpo, ngp, mposunit, method, chrm)
    if no_grp > 1:
        # if map is more than one, name list using group names
        mylist = dict(zip(probn, mylist))
    return mylist


@njit
def makemems(gmat, meff):
    """Set up family-specific marker effects (Mendelian sampling)."""
    qqq = np.zeros((gmat.shape))
    for i in range(gmat.shape[0]):
        for j in range(gmat.shape[1]):
            if gmat[i, j] == 4:
                qqq[i, j] = meff[j]*-1
            elif gmat[i, j] == 3:
                qqq[i, j] = meff[j]
            else:
                qqq[i, j] = 0
    return qqq


@njit
def makemebv(gmat, meff):
    """Set up family-specific marker effects (GEBV)."""
    qqq = np.zeros((gmat.shape))
    for i in range(gmat.shape[0]):
        for j in range(gmat.shape[1]):
            if gmat[i, j] == 2:
                qqq[i, j] = meff[j]*-1
            elif gmat[i, j] == 1:
                qqq[i, j] = meff[j]
            else:
                qqq[i, j] = 0
    return qqq


def traitspecmatrices(gmat, meff):
    """Store trait-specific matrices in a list."""
    notr = meff.shape[1]  # number of traits
    slist = []  # list stores trait-specific matrices
    slist.append([])
    for i in range(notr):
        # specify data type for numba
        mefff = np.array(meff.iloc[:, i], float)
        matrix_ms = makemems(gmat, mefff)
        slist[0].append(matrix_ms)
    return slist


def namesdf(notr, trait_names):
    """Create names of dataframe columns for Mendelian co(var)."""
    tnn = np.zeros((notr, notr), 'U20')
    tnn = np.chararray(tnn.shape, itemsize=30)
    for i in range(notr):
        for trt in range(notr):
            if i == trt:
                tnn[i, trt] = str(trait_names[i])
            elif i != trt:
                tnn[i, trt] = "{}_{}".format(trait_names[i], trait_names[trt])
    colnam = tnn[np.tril_indices(notr)]
    return colnam


def mrmmult(temp, covmat):
    """Matrix multiplication (MRM' or m'Rm)."""
    return temp @ covmat @ temp.T


def dgmrm(temp, covmat):
    """Matrix multiplication (MRM') for bigger matrices."""
    temp1111 = scipy.linalg.blas.dgemm(alpha=1.0, a=temp, b=covmat)
    return scipy.linalg.blas.dgemm(alpha=1.0, a=temp1111, b=temp.T)


def progr(itern, total):
    """Print progress of a task."""
    fill, printend, prefix, suffix = '█', "\r", 'Progress:', 'Complete'
    deci, length = 0, 50
    percent = ("{0:." + str(deci) + "f}").format(100 * (itern / float(total)))
    filledlen = int(length * itern // total)
    bars = fill * filledlen + '-' * (length - filledlen)
    print(f'\r{prefix} |{bars}| {percent}% {suffix}', end=printend)
    if itern == total:
        print()


def subindcheck(info, sub_id):
    """Check if inds provided in pd.DataFrame (sub_id) are in group data."""
    sub_id = pd.DataFrame(sub_id).reset_index(drop=True)
    if sub_id.shape[1] != 1:
        sys.exit("Individuals' IDs (sub_id) should be provided in one column")
    numbers = info.group.iloc[:, 1].astype(str).tolist()
    sub_id = sub_id.squeeze().astype(str).tolist()
    aaa = [numbers.index(x) if x in numbers else None for x in sub_id]
    aaa = np.array(aaa)
    if len(aaa) != len(sub_id):
        sys.exit("Some individual ID could not be found in group data")
    return aaa


def msvarcov_g_st(info, covmat, sub_id, progress=False):
    """Derive Mendelian sampling co(variance) for single trait."""
    if sub_id is not None:
        aaa = subindcheck(info, sub_id)
        idn = info.group.iloc[aaa, 1].reset_index(drop=True).astype(str)  # ID
        groupsex = info.group.iloc[aaa, 0].reset_index(drop=True).astype(str)
        matsub = info.gmat[aaa, :]
    else:
        idn = info.group.iloc[:, 1].reset_index(drop=True).astype(str)  # ID
        groupsex = info.group.iloc[:, 0].reset_index(drop=True).astype(str)
        matsub = info.gmat
    if (info.gmap.shape[1]-3 == 1 and len(pd.unique(groupsex)) > 1):
        print("The same map will be used for all groups")
    if progress:
        progr(0, matsub.shape[0])  # print progress bar
    snpindexxx = np.arange(start=0, stop=info.gmap.shape[0], step=1)
    notr = info.meff.columns.size
    slist = traitspecmatrices(matsub, info.meff)
    # dataframe to save Mendelian sampling (co)variance and aggregate breeding
    msvmsc = np.empty((matsub.shape[0], 1))
    for i in range(matsub.shape[0]):  # loop over no of individuals
        mscov = np.zeros((notr, notr))  # Mendelian co(var) mat for ind i
        for chrm in pd.unique(info.gmap.iloc[:, 0]):
            # snp index for chromosome chrm
            s_ind = np.array(snpindexxx[info.gmap.iloc[:, 0] == (chrm)])
            # family-specific marker effects for ind i
            temp = np.zeros((notr, len(s_ind)))
            for trt in range(notr):
                temp[trt, :] = slist[0][trt][i, s_ind]
            if info.gmap.shape[1]-3 == 1:
                mscov = mscov + mrmmult(temp, covmat[0][chrm-1])
            else:
                mscov = mscov + mrmmult(temp, covmat[groupsex[i]][chrm-1])
        msvmsc[i, 0] = mscov
        if progress:
            progr(i + 1, matsub.shape[0])  # print progress bar
    msvmsc = pd.DataFrame(msvmsc)
    msvmsc.columns = info.meff.columns
    msvmsc.insert(0, "ID", idn, True)
    msvmsc.insert(1, "Group", groupsex, True)  # insert group
    return msvmsc


def msvarcov_g_mt(info, covmat, sub_id, progress=False):
    """Derive Mendelian sampling co(variance) for multiple traits."""
    if sub_id is not None:
        aaa = subindcheck(info, sub_id)
        idn = info.group.iloc[aaa, 1].reset_index(drop=True).astype(str)  # ID
        groupsex = info.group.iloc[aaa, 0].reset_index(drop=True).astype(str)
        matsub = info.gmat[aaa, :]
    else:
        idn = info.group.iloc[:, 1].reset_index(drop=True).astype(str)  # ID
        groupsex = info.group.iloc[:, 0].reset_index(drop=True).astype(str)
        matsub = info.gmat
    if (info.gmap.shape[1]-3 == 1 and len(pd.unique(groupsex)) > 1):
        print("The same map will be used for all groups")
    if progress:
        progr(0, matsub.shape[0])  # print progress bar
    snpindexxx = np.arange(start=0, stop=info.gmap.shape[0], step=1)
    notr = info.meff.columns.size
    slist = traitspecmatrices(matsub, info.meff)
    # dataframe to save Mendelian sampling (co)variance and aggregate breeding
    mad = len(np.zeros((notr+1, notr+1))[np.tril_indices(notr+1)])
    msvmsc = np.empty((matsub.shape[0], mad))
    for i in range(matsub.shape[0]):  # loop over no of individuals
        mscov = np.zeros((notr+1, notr+1))  # Mendelian co(var) mat for ind i
        for chrm in pd.unique(info.gmap.iloc[:, 0]):
            # snp index for chromosome chrm
            s_ind = np.array(snpindexxx[info.gmap.iloc[:, 0] == (chrm)])
            # family-specific marker effects for ind i
            temp = np.zeros((notr+1, len(s_ind)))
            for trt in range(notr):
                temp[trt, :] = slist[0][trt][i, s_ind]
                temp[notr, :] = np.matmul(info.indwt.T, temp[0:notr, :])
            if info.gmap.shape[1]-3 == 1:
                mscov = mscov + mrmmult(temp, covmat[0][chrm-1])
            else:
                mscov = mscov + mrmmult(temp, covmat[groupsex[i]][chrm-1])
        msvmsc[i, :] = mscov[np.tril_indices(notr+1)]
        if progress:
            progr(i + 1, matsub.shape[0])  # print progress bar
    msvmsc = pd.DataFrame(msvmsc)
    tnames = np.concatenate((info.meff.columns, "AG"), axis=None)
    colnam = namesdf(notr+1, tnames).decode('utf-8')
    msvmsc.columns = colnam
    msvmsc.insert(0, "ID", idn, True)
    msvmsc.insert(1, "Group", groupsex, True)  # insert group
    return msvmsc


def msvarcov_g(info, covmat, sub_id, progress=False):
    """
    Derive Mendelian sampling co(variance) and aggregate genotype.

    Parameters
    ----------
    info : class object
        A class object created using the function "datacheck"
    covmat : A list of pop cov matrices created using "popcovmat" function
    sub_id : pandas.DataFrame with one column
        Index: RangeIndex (minimum of 2 rows)
        Containing ID numbers of specific individuals to be evaluated
    progress : bool, optional; print progress of the function if True
    Returns
    -------
    msvmsc : pandas.DataFrame
        containing the Mendelian sampling (co)variance and aggregate genotype
    Note: If sub_id is None, Mendelian (co-)variance will be estimated for
    all individuals. Otherwise, Mendelian (co-)variance will be estimated for
    the individuals in sub_id
    """
    notr = info.meff.columns.size
    if notr == 1:
        msvmsc = msvarcov_g_st(info, covmat, sub_id, progress)
    elif notr > 1:
        msvmsc = msvarcov_g_mt(info, covmat, sub_id, progress)
    return msvmsc


def array2sym(array):
    """Convert array to stdized symm mat, and back to array without diags."""
    dfmsize = array.size
    for notr in range(1, 10000):
        if dfmsize == len(np.zeros((notr, notr))[np.tril_indices(notr)]):
            break
    iii, jjj = np.tril_indices(notr)
    mat = np.empty((notr, notr), float)
    mat[iii, jjj], mat[jjj, iii] = array, array
    mat = np.array(mat)
    mat1 = cov2corr(mat)
    return np.array(mat1[np.tril_indices(notr, k=-1)])


def msvarcov_gcorr(msvmsc):
    """
    Standardize Mendelian sampling co(variance) and aggregate genotype.

    Parameters
    ----------
    msvmsc : pandas.DataFrame
        containing the Mendelian sampling (co)variance and aggregate genotype
        created using msvarcov_g function
    Returns
    -------
    dfcor : pandas.DataFrame
        containing standardized Mendelian sampling (co)variance
    """
    if msvmsc.columns.size == 3:
        sys.exit("Correlation cannot be derived for a single trait")
    dfm = msvmsc.iloc[:, 2:msvmsc.shape[1]]  # exclude ID and group
    dfmsize = dfm.shape[1]
    # derive number of traits
    for notr in range(1, 10000):
        if dfmsize == len(np.zeros((notr, notr))[np.tril_indices(notr)]):
            break
    # standardize covariance between traits
    dfcor = dfm.apply(array2sym, axis=1)
    # extract column names
    listnames = dfm.columns.tolist()
    cnames = [x for x in listnames if "_" in x]
    # convert pd.series of list to data frame
    dfcor = pd.DataFrame.from_dict(dict(zip(dfcor.index, dfcor.values))).T
    dfcor.columns = cnames
    # insert ID and group info
    dfcor = [pd.DataFrame(msvmsc.iloc[:, 0:2]), dfcor]  # add ID and GRP
    dfcor = pd.concat(dfcor, axis=1)
    return dfcor


def calcgbv(info, sub_id):
    """Calculate breeding values for each trait."""
    if sub_id is not None:
        aaa = subindcheck(info, sub_id)
        idn = info.group.iloc[aaa, 1].reset_index(drop=True).astype(str)  # ID
        groupsex = info.group.iloc[aaa, 0].reset_index(drop=True).astype(str)
        matsub = info.gmat[aaa, :]
    else:
        idn = info.group.iloc[:, 1].reset_index(drop=True).astype(str)  # ID
        groupsex = info.group.iloc[:, 0].reset_index(drop=True).astype(str)
        matsub = info.gmat
    no_individuals = matsub.shape[0]  # Number of individuals
    trait_names = info.meff.columns  # traits names
    notr = trait_names.size  # number of traits
    if notr == 1:
        gbv = np.zeros((no_individuals, notr))
        mefff = np.array(info.meff.iloc[:, 0], float)  # type spec for numba
        matrix_me = makemebv(matsub, mefff)  # fam-spec marker effects BV
        gbv[:, 0] = matrix_me.sum(axis=1)  # sum all effects
        gbv = pd.DataFrame(gbv)
        gbv.columns = trait_names
    elif notr > 1:
        gbv = np.zeros((no_individuals, notr+1))
        for i in range(notr):
            mefff = np.array(info.meff.iloc[:, i], float)  # type spec 4 numba
            matrix_me = makemebv(matsub, mefff)  # fam-spec marker effects BV
            gbv[:, i] = matrix_me.sum(axis=1)  # sum all effects for each trait
            gbv[:, notr] = gbv[:, notr] + info.indwt[i]*gbv[:, i]  # Agg gen
        gbv = pd.DataFrame(gbv)
        colnames = np.concatenate((trait_names, "ABV"), axis=None)
        gbv.columns = colnames
    gbv.insert(0, "ID", idn, True)  # insert ID
    gbv.insert(1, "Group", groupsex, True)  # insert group
    return gbv


def calcprob(info, msvmsc, thresh):
    """Calculate the probability of breeding top individuals."""
    aaa = subindcheck(info, pd.DataFrame(msvmsc.iloc[:, 0]))
    gbvall = calcgbv(info, None)  # calc GEBV for all inds used by thresh
    gbv = gbvall.iloc[aaa, :].reset_index(drop=True)  # GEBV matching msvmsc
    no_individuals = gbv.shape[0]  # Number of individuals
    trait_names = info.meff.columns  # traits names
    notr = trait_names.size  # number of traits
    if notr == 1:
        probdf = np.zeros((no_individuals, notr))
        ttt = np.quantile(gbvall.iloc[:, (0+2)], q=1-thresh)  # threshold
        probdf[:, 0] = 1 - scipy.stats.norm.cdf(
            ttt, loc=gbv.iloc[:, (0+2)], scale=np.sqrt(msvmsc.iloc[:, 0+2]))
        probdf = pd.DataFrame(probdf)
        probdf.columns = trait_names
    elif notr > 1:
        colnam = np.concatenate((info.meff.columns, "AG"), axis=None)
        colnam = namesdf(notr+1, colnam).decode('utf-8')
        ttt = np.quantile(gbvall.iloc[:, (notr+2)], q=1-thresh)  # threshold
        probdf = np.zeros((no_individuals, notr+1))
        t_ind = np.arange(colnam.shape[0])[np.in1d(colnam, trait_names)]
        for i in range(notr):
            ttt = np.quantile(gbvall.iloc[:, (i+2)], q=1-thresh)  # threshold
            probdf[:, i] = scipy.stats.norm.cdf(
                ttt, loc=gbv.iloc[:, (i+2)], scale=np.sqrt(
                    msvmsc.iloc[:, (t_ind[i])+2]))
            probdf[:, i] = np.nan_to_num(probdf[:, i])  # convert Inf to zero
            probdf[:, i] = 1 - probdf[:, i]
        ttt = np.quantile(gbvall.iloc[:, (notr+2)], q=1-thresh)
        probdf[:, notr] = scipy.stats.norm.cdf(
            ttt, loc=gbv.iloc[:, (notr+2)], scale=np.sqrt(
                msvmsc["AG"]))
        probdf[:, notr] = np.nan_to_num(probdf[:, notr])  # Agg
        probdf[:, notr] = 1 - probdf[:, notr]
        probdf = pd.DataFrame(probdf)  # convert matrix to dataframe
        colnames = np.concatenate((trait_names, "ABV"), axis=None)
        probdf.columns = colnames
    probdf = [pd.DataFrame(gbv.iloc[:, 0:2]), probdf]  # add ID and GRP
    probdf = pd.concat(probdf, axis=1)
    return probdf


def calcindex(info, msvmsc, const):
    """Calculate the index if constant is known."""
    sub_id = pd.DataFrame(msvmsc.iloc[:, 0])
    gbv = calcgbv(info, sub_id)  # calc GEBV
    no_individuals = gbv.shape[0]  # Number of individuals
    trait_names = info.meff.columns  # traits names
    notr = trait_names.size
    if notr == 1:
        indexdf = np.zeros((no_individuals, notr))
        indexdf[:, 0] = (gbv.iloc[:, (0+2)]/2) + np.sqrt(
            msvmsc.iloc[:, 0+2])*const
        indexdf = pd.DataFrame(indexdf)
        indexdf.columns = trait_names
    elif notr > 1:
        colnam = np.concatenate((info.meff.columns, "AG"), axis=None)
        colnam = namesdf(notr+1, colnam).decode('utf-8')
        indexdf = np.zeros((no_individuals, notr+1))
        t_ind = np.arange(colnam.shape[0])[np.in1d(colnam, trait_names)]
        for i in range(notr):
            indexdf[:, i] = (gbv.iloc[:, (i+2)]/2) + np.sqrt(
                msvmsc.iloc[:, (t_ind[i]+2)])*const
        indexdf[:, notr] = (gbv.iloc[:, (notr+2)]/2) + np.sqrt(
            msvmsc["AG"])*const
        indexdf = pd.DataFrame(indexdf)
        colnames = np.concatenate((trait_names, "ABV"), axis=None)
        indexdf.columns = colnames
    indexdf = [pd.DataFrame(gbv.iloc[:, 0:2]), indexdf]  # add ID and GRP
    indexdf = pd.concat(indexdf, axis=1)
    return indexdf


def selstrat_g(selstrat, info, sub_id, msvmsc, throrconst):
    """
    Calc selection criteria (GEBV, PBTI, or index using gametic approach.

    Parameters
    ----------
    selstrat : str
        A str containing any of GEBV, PBTI or index
    info : class object
        A class object created using the function "datacheck"
    sub_id : pandas.DataFrame with one column
        Index: RangeIndex (minimum of 2 rows)
        Containing ID numbers of specific individuals to be evaluated
    msvmsc : pandas.DataFrame
        DF created using the function "msvarcov_g"
    throrconst : float
        If selstrat is PBTI, a throrconst of value 0.05 sets threshold at
        top 5% of GEBV. If selstrat is index, throrconst is a constant.
        If selstrat is GEBV, throrconst can be any random value.

    Returns
    -------
    data : pandas.DataFrame
        Index: RangeIndex
        Columns:
        ID, Group, trait names and Aggregate Breeding Value (ABV)
    Note: If selstrat is GEBV, None may be used for throrconst and msvmsc.
    If sub_id is None and selstrat is GEBV, GEBVs will be estimated for all
    individuals. However, if selstrat is not GEBV, the chosen selection
    criterion will be estimated for all individuals in msvmsc data frame.
    """
    if selstrat in ("PBTI", "pbti", "index", "INDEX") and msvmsc is None:
        sys.exit("Provide Mendelian (co-)variance dataframe: 'msvmsc'")
    if selstrat in ("PBTI", "pbti", "index", "INDEX") and throrconst is None:
        sys.exit("Provide value for throrconst parameter")
    if selstrat not in ('GEBV', 'gebv', 'PBTI', 'pbti', 'index', 'INDEX'):
        sys.exit("selection strategy should be one of GEBV, PBTI or INDEX")
    if selstrat in ('GEBV', 'gebv'):
        data = calcgbv(info, sub_id)
    elif selstrat in ('PBTI', 'pbti'):
        if throrconst > 1 or throrconst < 0:
            sys.exit("value must be in the range of 0 and 1")
        data = calcprob(info, msvmsc, throrconst)
    elif selstrat in ('index', 'INDEX'):
        data = calcindex(info, msvmsc, throrconst)
    return data


def cov2corr(cov):
    """Convert covariance to correlation matrix."""
    cov = np.asanyarray(cov)
    std_ = np.sqrt(np.diag(cov))
    with np.errstate(invalid='ignore'):
        corr = cov / np.outer(std_, std_)
    return corr


def aggen(us_ind, no_markers, slst, indwt):
    """Set up additive effects matrix of aggregate genotype."""
    mmfinal = np.empty((len(us_ind), no_markers))
    xxx = 0
    for iii in us_ind:
        tmpmt1 = np.array([slst[0][trt][iii, :] for trt in range(indwt.size)])
        mmfinal[xxx, :] = np.matmul(indwt.transpose(), tmpmt1)
        xxx = xxx + 1
    return mmfinal


def chr_int(xxxxx):
    """Format chromomosome of interest parameter."""
    if 'all' in xxxxx:
        xxxxx = 'all'
    elif 'none' in xxxxx:
        xxxxx = 'none'
    else:
        xxxxx = np.array([int(i) for i in xxxxx])
    return xxxxx


def writechr(covtmpx, chrinterest, chrm, trtnam, probx, stdsim):
    """Write matrices to file."""
    if isinstance(chrinterest, str):
        if chrinterest == 'all':
            chrfile1 = "{}/Sim mat for {} chrm {} grp {}.npy".format(
                os.getcwd(), trtnam, chrm, probx)
            np.save(chrfile1, covtmpx)
    elif chrm in chrinterest:
        chrfile1 = "{}/Sim mat for {} chrm {} grp {}.npy".format(
            os.getcwd(), trtnam, chrm, probx)  # output file
        np.save(chrfile1, covtmpx)
    if stdsim:
        if isinstance(chrinterest, str):
            if chrinterest == 'all':
                chrfilec = "{}/Stdsim mat for {} chrm {} grp {}.npy".format(
                    os.getcwd(), trtnam, chrm, probx)  # output file
                np.save(chrfilec, cov2corr(covtmpx))
        elif chrm in chrinterest:
            chrfilec = "{}/Stdsim mat for {} chrm {} grp {}.npy".format(
                os.getcwd(), trtnam, chrm, probx)  # output file
            np.save(chrfilec, cov2corr(covtmpx))


def writechrunspec(covtmpx, chrinterest, chrm, trtnam, stdsim):
    """Write matrices to file."""
    if isinstance(chrinterest, str):
        if chrinterest == 'all':
            chrfile1 = "{}/Sim mat for {} chrm {}.npy".format(
                os.getcwd(), trtnam, chrm)
            np.save(chrfile1, covtmpx)
    elif chrm in chrinterest:
        chrfile1 = "{}/Sim mat for {} chrm {}.npy".format(
            os.getcwd(), trtnam, chrm)  # output file
        np.save(chrfile1, covtmpx)
    if stdsim:
        if isinstance(chrinterest, str):
            if chrinterest == 'all':
                chrfilec = "{}/Stdsim mat for {} chrm {}.npy".format(
                    os.getcwd(), trtnam, chrm)  # output file
                np.save(chrfilec, cov2corr(covtmpx))
        elif chrm in chrinterest:
            chrfilec = "{}/Stdsim mat for {} chrm {}.npy".format(
                os.getcwd(), trtnam, chrm)  # output file
            np.save(chrfilec, cov2corr(covtmpx))


def grtonum(numnx):
    """Map chracters to numeric (0-no of groups)."""
    numnx = numnx.reset_index(drop=True)
    probn = pd.unique(numnx).tolist()
    alt_no = np.arange(0, len(probn), 1)
    noli = numnx.tolist()
    numnx = np.array(list(map(dict(zip(probn, alt_no)).get, noli, noli)))
    return numnx, probn


def datret(info, rw_nms, pfnp, us_ind, slist, covmat, cov_indxx, stdsim,
           progress):
    """Return sim mat based on aggregate genotypes."""
    snpindexxxx = np.arange(start=0, stop=info.gmap.shape[0], step=1)
    if info.meff.shape[1] == 1 and not stdsim:
        mat = cov_indxx
    elif info.meff.shape[1] == 1 and stdsim:
        mat = cov2corr(cov_indxx)
    elif info.meff.shape[1] > 1:
        if info.gmap.shape[1]-3 > 1:
            rw_nms = pd.DataFrame(rw_nms)
            rw_nms.to_csv(f"order of inds in mat grp {pfnp}.csv", index=False)
        if progress:
            print('Creating similarity matrix based on aggregate genotype')
            progr(0, max(pd.unique(info.gmap.iloc[:, 0])))
        tmpmt1 = aggen(us_ind, info.gmap.shape[0], slist, info.indwt)
        # stores ABV covariance btw inds
        mat = np.zeros((len(us_ind), len(us_ind)))
        # loop over chromososomes
        for chrm in pd.unique(info.gmap.iloc[:, 0]):
            s_ind = np.array(snpindexxxx[info.gmap.iloc[:, 0] == (chrm)])
            if info.gmap.shape[1]-3 == 1:
                covtmpx = abs(dgmrm(tmpmt1[:, s_ind], covmat[0][chrm-1]))
            else:
                covtmpx = abs(dgmrm(tmpmt1[:, s_ind], covmat[pfnp][chrm-1]))
            mat = mat + covtmpx
            if progress:
                progr(chrm, max(pd.unique(info.gmap.iloc[:, 0])))
        if stdsim:
            mat = cov2corr(mat)
    return mat


def mrmcals(info, us_ind, stdsim, slist, covmat, probn, chrinterest, save,
            progress):
    """Compute similarity matrix for each chromosome."""
    if progress:
        progr(0, info.meff.columns.size)
    for i in range(info.meff.columns.size):
        cov_indxx = np.zeros((len(us_ind), len(us_ind)))
        for chrm in pd.unique(info.gmap.iloc[:, 0]):
            s_ind = np.array(np.arange(0, info.gmap.shape[0], 1
                                       )[info.gmap.iloc[:, 0] == (chrm)])
            if info.gmap.shape[1]-3 == 1:  # map is 1
                covtmpx = abs(dgmrm(slist[0][i][:, s_ind], covmat[0][chrm-1]))
            else:  # if map is more than 1
                covtmpx = abs(dgmrm(slist[0][i][us_ind[:, None], s_ind],
                                    covmat[probn][chrm-1]))
            cov_indxx = cov_indxx + covtmpx  # sums up chrm-specific sims
            if len(pd.unique(info.group.iloc[:, 0].astype(str))) == 1:
                writechrunspec(covtmpx, chrinterest, chrm,
                               info.meff.columns[i], stdsim)
            else:
                writechr(covtmpx, chrinterest, chrm, info.meff.columns[i],
                         probn, stdsim)  # write sim to file
        if stdsim:
            if save is True:
                if info.gmap.shape[1]-3 == 1:
                    covxfile = "{}/Stdsim mat for {}.npy".format(
                        os.getcwd(), info.meff.columns[i])
                else:
                    covxfile = "{}/Stdsim mat for {} grp {}.npy".format(
                        os.getcwd(), info.meff.columns[i], probn)
                np.save(covxfile, cov2corr(cov_indxx))  # write std sim mats
        else:
            if save is True:
                if info.gmap.shape[1]-3 == 1:
                    covxfile = "{}/Sim mat for {}.npy".format(
                        os.getcwd(), info.meff.columns[i])
                else:
                    covxfile = "{}/Sim mat for {} grp {}.npy".format(
                        os.getcwd(), info.meff.columns[i], probn)
                np.save(covxfile, cov_indxx)  # write sim matrices
        if progress:
            progr(i + 1, info.meff.columns.size)
    return cov_indxx


def simmat_g(info, covmat, sub_id, chrinterest, save=False, stdsim=False,
             progress=False):
    """
    Compute similarity matrices using gametic approach.

    Parameters
    ----------
    info : class object
        A class object created using the function "datacheck"
    covmat : A list of pop cov matrices created using "popcovmat" function
    sub_id : pandas.DataFrame with one column
        Index: RangeIndex (minimum of 2 rows)
        Containing ID numbers of specific individuals to be evaluated
    chrinterest : str or list of int
        list of chromosome numbers of interest or str with "all" or "none"
    save : bool, optional; write trait-specific sim mats to file if true
    stdsim : bool, optional; print write std sim mats to file if true
    progress : bool, optional; print progress of the task if true
    Returns
    -------
    multgrpcov : list containing simimlarity matrices for each group
    """
    if sub_id is None:
        inda = np.arange(0, info.gmat.shape[0], 1)
        sub_id = pd.DataFrame(info.group.iloc[inda, 1])
        aaa = subindcheck(info, sub_id)
    else:
        aaa = subindcheck(info, sub_id)
    chrinterest = chr_int(chrinterest)
    slist = traitspecmatrices(info.gmat[aaa, :], info.meff)  # trt-spec mat
    grp = info.gmap.shape[1]-3
    if (grp == 1 and len(pd.unique(info.group.iloc[:, 0].astype(str))) > 1):
        print("The same map will be used for all groups")
    numbers, probn = grtonum(info.group.iloc[aaa, 0].astype(str))
    multgrpcov = []
    for gnp in range(grp):
        multgrpcov.append([])
        if grp == 1:
            us_ind = np.arange(start=0, stop=info.gmat[aaa, :].shape[0],
                               step=1)
        else:
            tng = numbers == gnp
            us_ind = np.array(list(compress(np.arange(0, len(tng), 1),
                                            tng))).T
            print("Processing group ", probn[gnp])
        rw_nms = info.group.iloc[aaa, 1].reset_index(drop=True).astype(
            str)[us_ind]
        cov_indxx = mrmcals(info, us_ind, stdsim, slist, covmat, probn[gnp],
                            chrinterest, save, progress)
        multgrpcov[int(gnp)].append(
            datret(info, rw_nms, probn[gnp], us_ind, slist, covmat,
                   cov_indxx, stdsim, progress))
        if len(probn) == 1:
            break
    if grp > 1 and len(probn):
        multgrpcov = dict(zip(probn, multgrpcov))
    return multgrpcov


def submsvmsc(msvmsc, sub_idz):
    """Extract index in msvmsc data frame."""
    sub_idz = pd.DataFrame(sub_idz)
    numbs = msvmsc.iloc[:, 0].astype(str).tolist()
    sub_idz = sub_idz.reset_index(drop=True).squeeze()
    mal = sub_idz.iloc[:, 0].astype(str).tolist()
    fem = sub_idz.iloc[:, 1].astype(str).tolist()
    if sub_idz is not None:
        for i in mal:
            if i not in numbs:
                sys.exit("Individuals are not in msvmsc parameter")
        for i in fem:
            if i not in numbs:
                sys.exit("Individuals are not in msvmsc parameter")
    mal1 = [numbs.index(x) if x in numbs else None for x in mal]
    fem1 = [numbs.index(x) if x in numbs else None for x in fem]
    return mal1, fem1


def pot_parents(info, data, selmale, selfm):
    """Subset individuals of interest."""
    trait_names = info.meff.columns
    if trait_names.size == 1:
        datamale = data[data.iloc[:, 1] == selmale[0]]
        pos = subindcheck(info, pd.DataFrame(datamale.iloc[:, 0]))
        datamale.insert(0, "pos", pos, True)
        no_sire = int(datamale.shape[0] * selmale[1])
        datamale = datamale.sort_values(
            by=[trait_names[0]], ascending=False).iloc[0:no_sire, :]
        datafemale = data[data.iloc[:, 1] == selfm[0]]
        pos = subindcheck(info, pd.DataFrame(datafemale.iloc[:, 0]))
        datafemale.insert(0, "pos", pos, True)
        no_dam = int(datafemale.shape[0] * selfm[1])
        datafemale = datafemale.sort_values(
            by=[trait_names[0]], ascending=False).iloc[0:no_dam, :]
    elif trait_names.size > 1:
        datamale = data[data.iloc[:, 1] == selmale[0]]
        pos = subindcheck(info, pd.DataFrame(datamale.iloc[:, 0]))
        datamale.insert(0, "pos", pos, True)
        no_sire = int(datamale.shape[0] * selmale[1])
        datamale = datamale.sort_values(
            by=['ABV'], ascending=False).iloc[0:no_sire, :]
        datafemale = data[data.iloc[:, 1] == selfm[0]]
        pos = subindcheck(info, pd.DataFrame(datafemale.iloc[:, 0]))
        datafemale.insert(0, "pos", pos, True)
        no_dam = int(datafemale.shape[0] * selfm[1])
        datafemale = datafemale.sort_values(
            by=['ABV'], ascending=False).iloc[0:no_dam, :]
    matlist = np.array(np.meshgrid(
        datamale.iloc[:, 0], datafemale.iloc[:, 0])).T.reshape(-1, 2)
    ids = np.array(np.meshgrid(
        datamale.iloc[:, 1], datafemale.iloc[:, 1])).T.reshape(-1, 2)
    if trait_names.size == 1:
        matndat = pd.DataFrame(index=range(matlist.shape[0]), columns=range(
            4+trait_names.size))
    else:
        matndat = pd.DataFrame(
            index=range(matlist.shape[0]), columns=range(5+trait_names.size))
    matndat.iloc[:, [0, 1]] = ids
    matndat.iloc[:, [2, 3]] = matlist
    return matndat


def selsgebv(notr, matndat, gbv, maxmale):
    """Calculate breeding values for each trait (zygote)."""
    mal = matndat.iloc[:, 2].tolist()
    fem = matndat.iloc[:, 3].tolist()
    if notr == 1:
        matndat.iloc[:, 4] = (np.array(gbv.iloc[mal, (0+2)]) + np.array(
            gbv.iloc[fem, (0+2)]))/2
    elif notr > 1:
        matndat.iloc[:, 4:(5+notr)] = (np.array(
            gbv.iloc[mal, 2:(notr+3)]) + np.array(gbv.iloc[fem, 2:(notr+3)]))/2
    idfxxx = np.unique(matndat.iloc[:, 3])
    mmat = pd.DataFrame(index=range(len(idfxxx)),
                        columns=range(matndat.shape[1]))
    for mmm in np.arange(0, len(idfxxx), 1):
        axx = matndat.loc[matndat.iloc[:, 3] == idfxxx[mmm]]
        tsire = np.array(axx.iloc[:, 2])
        mmat.iloc[mmm, :] = axx.iloc[np.argmax(
            axx.iloc[:, axx.columns.size-1]), :]
        norepssire = Counter(mmat.iloc[:, 2])
        lents = len(tsire)
        for nrs in range(lents):
            if norepssire[tsire[nrs]] <= maxmale-1:
                mmat.iloc[mmm, :] = np.array(axx[axx.iloc[:, 2] == tsire[nrs]])
                break
    matndat = mmat
    if notr == 1:
        matndat.columns = np.concatenate((
            ['MaleID', 'FemaleID', 'MaleIndex', 'FemaleIndex'],
            gbv.columns[gbv.columns.size-1]), axis=None)
    else:
        matndat.columns = np.concatenate((
            ['MaleID', 'FemaleID', 'MaleIndex', 'FemaleIndex'],
            gbv.columns[2:gbv.columns.size].tolist()), axis=None)
    return matndat


def selspbtizyg(notr, gbv, matndat, msvmsc, throrconst, maxmale):
    """Calculate prob of breeding top inds (zygote)."""
    mal1, fem1 = submsvmsc(msvmsc, pd.DataFrame(matndat.iloc[:, 0:2]))
    mal = matndat.iloc[:, 2].tolist()
    fem = matndat.iloc[:, 3].tolist()
    if notr == 1:
        matndat.iloc[:, 4] = (np.array(gbv.iloc[mal, (0+2)]) + np.array(
                gbv.iloc[fem, (0+2)]))/2
        ttt = np.quantile(gbv.iloc[:, 0+2], q=1-throrconst)
        msvtemp = np.array(msvmsc.iloc[mal1, 0+2]) + np.array(
            msvmsc.iloc[fem1, 0+2])
        matndat.iloc[:, 4] = 1 - scipy.stats.norm.cdf(
            ttt, loc=matndat.iloc[:, 4], scale=np.sqrt(
                msvtemp))
    elif notr > 1:
        trait_names = gbv.columns[2:2+notr]
        colnam = np.concatenate((trait_names, "AG"), axis=None)
        colnam = namesdf(notr+1, colnam).decode('utf-8')
        t_ind = np.arange(colnam.shape[0])[np.in1d(colnam, trait_names)]
        for i in range(notr):
            matndat.iloc[:, 4+i] = (
                np.array(gbv.iloc[mal, (i+2)]) + np.array(
                    gbv.iloc[fem, (i+2)]))/2
            ttt = np.quantile(gbv.iloc[:, 2+i], q=1-throrconst)
            msvtemp = np.array(msvmsc.iloc[mal1, t_ind[i]+2]) + np.array(
                msvmsc.iloc[fem1, t_ind[i]+2])
            matndat.iloc[:, 4+i] = 1 - scipy.stats.norm.cdf(
                ttt, loc=matndat.iloc[:, 4+i], scale=np.sqrt(msvtemp))
        matndat.iloc[:, 4+notr] = (
            np.array(gbv.iloc[mal, (notr+2)]) + np.array(
                    gbv.iloc[fem, (notr+2)]))/2
        ttt = np.quantile(gbv.iloc[:, 2+notr], q=1-throrconst)
        msvtemp = np.array(msvmsc.loc[mal1, ["AG"]]) + np.array(
            msvmsc.loc[fem1, ["AG"]])
        matndat.iloc[:, 4+notr] = 1 - scipy.stats.norm.cdf(
            ttt, loc=matndat.iloc[:, 4+notr], scale=np.sqrt(msvtemp.ravel()))
    idfxxx = np.unique(matndat.iloc[:, 3])
    mmat = pd.DataFrame(index=range(len(idfxxx)),
                        columns=range(matndat.shape[1]))
    for mmm in np.arange(0, len(idfxxx), 1):
        axx = matndat.loc[matndat.iloc[:, 3] == idfxxx[mmm]]
        tsire = np.array(axx.iloc[:, 2])
        mmat.iloc[mmm, :] = axx.iloc[np.argmax(
            axx.iloc[:, axx.columns.size-1]), :]
        norepssire = Counter(mmat.iloc[:, 2])
        lents = len(tsire)
        for nrs in range(lents):
            if norepssire[tsire[nrs]] <= maxmale-1:
                mmat.iloc[mmm, :] = np.array(axx[axx.iloc[:, 2] == tsire[nrs]])
                break
    matndat = mmat
    if notr == 1:
        matndat.columns = np.concatenate((
            ['MaleID', 'FemaleID', 'MaleIndex', 'FemaleIndex'],
            gbv.columns[gbv.columns.size-1]), axis=None)
    else:
        matndat.columns = np.concatenate((
            ['MaleID', 'FemaleID', 'MaleIndex', 'FemaleIndex'],
            gbv.columns[2:gbv.columns.size].tolist()), axis=None)
    return matndat


def selsindex(notr, gbv, matndat, msvmsc, throrconst, maxmale):
    """Calculate the index if constant is known (zygote)."""
    mal1, fem1 = submsvmsc(msvmsc, pd.DataFrame(matndat.iloc[:, 0:2]))
    mal = matndat.iloc[:, 2].tolist()
    fem = matndat.iloc[:, 3].tolist()
    if notr == 1:
        matndat.iloc[:, 4] = (np.array(gbv.iloc[mal, (0+2)]) + np.array(
                gbv.iloc[fem, (0+2)]))/2
        msvtemp = np.array(msvmsc.iloc[mal1, 0+2]) + np.array(
                msvmsc.iloc[fem1, 0+2])
        matndat.iloc[:, 4] = matndat.iloc[:, 4] + np.sqrt(msvtemp)*throrconst
    elif notr > 1:
        trait_names = gbv.columns[2:2+notr]
        colnam = np.concatenate((trait_names, "AG"), axis=None)
        colnam = namesdf(notr+1, colnam).decode('utf-8')
        t_ind = np.arange(colnam.shape[0])[np.in1d(colnam, trait_names)]
        for i in range(notr):
            matndat.iloc[:, 4+i] = (
                np.array(gbv.iloc[mal, (i+2)]) + np.array(
                    gbv.iloc[fem, (i+2)]))/2
            msvtemp = np.array(msvmsc.iloc[mal1, t_ind[i]+2]) + np.array(
                msvmsc.iloc[fem1, t_ind[i]+2])
            matndat.iloc[:, 4+i] = matndat.iloc[:, 4+i] + np.sqrt(
                msvtemp)*throrconst
        matndat.iloc[:, 4+notr] = (
             np.array(gbv.iloc[mal, (notr+2)]) + np.array(
                 gbv.iloc[fem, (notr+2)]))/2
        msvtemp = np.array(msvmsc.loc[mal1, ["AG"]]) + np.array(
            msvmsc.loc[fem1, ["AG"]])
        matndat.iloc[:, 4+notr] = matndat.iloc[:, 4+notr] + (
            np.sqrt(msvtemp)*throrconst).ravel()
    idfxxx = np.unique(matndat.iloc[:, 3])
    mmat = pd.DataFrame(index=range(len(idfxxx)),
                        columns=range(matndat.shape[1]))
    for mmm in np.arange(0, len(idfxxx), 1):
        axx = matndat.loc[matndat.iloc[:, 3] == idfxxx[mmm]]
        tsire = np.array(axx.iloc[:, 2])
        mmat.iloc[mmm, :] = axx.iloc[np.argmax(
            axx.iloc[:, axx.columns.size-1]), :]
        norepssire = Counter(mmat.iloc[:, 2])
        lents = len(tsire)
        for nrs in range(lents):
            if norepssire[tsire[nrs]] <= maxmale-1:
                mmat.iloc[mmm, :] = np.array(axx[axx.iloc[:, 2] == tsire[nrs]])
                break
    matndat = pd.DataFrame(mmat)
    if notr == 1:
        matndat.columns = np.concatenate((
            ['MaleID', 'FemaleID', 'MaleIndex', 'FemaleIndex'],
            gbv.columns[gbv.columns.size-1]), axis=None)
    else:
        matndat.columns = np.concatenate((
            ['MaleID', 'FemaleID', 'MaleIndex', 'FemaleIndex'],
            gbv.columns[2:gbv.columns.size].tolist()), axis=None)
    return matndat


def subindcheckzyg(info, sub_idz):
    """Check sex and if matepairs provided in sub_idz are in group data."""
    numbs = info.group.iloc[:, 1].astype(str).tolist()
    sub_idz = pd.DataFrame(sub_idz).reset_index(drop=True).squeeze()
    mal = sub_idz.iloc[:, 0].astype(str).tolist()
    fem = sub_idz.iloc[:, 1].astype(str).tolist()
    mal1 = [numbs.index(x) if x in numbs else None for x in mal]
    fem1 = [numbs.index(x) if x in numbs else None for x in fem]
    if len(pd.unique(info.group.iloc[mal1, 0])) != 1:
        sys.exit("Group class in sub_idz is not unique to ID of males")
    if len(pd.unique(info.group.iloc[fem1, 0])) != 1:
        sys.exit("Group class in sub_idz is not unique to ID of females")
    idn = sub_idz.reset_index(drop=True)
    mgp = list(set(info.group.iloc[mal1, 0]))
    fgp = list(set(info.group.iloc[fem1, 0]))
    if len(mgp) > 1 or len(fgp) > 1:
        sys.exit("multiple sexes detected in data")
    probn = [mgp[0], fgp[0]]
    return mal1, fem1, idn, probn


def calcgbvzygsub(info, sub_idz):
    """Calc breeding values for matepairs."""
    mal1, fem1, idn, _ = subindcheckzyg(info, sub_idz)
    no_individuals, trait_names = idn.shape[0], info.meff.columns
    notr = trait_names.size
    if notr == 1:
        gbv = np.zeros((no_individuals, notr))
        mefff = np.array(info.meff.iloc[:, 0], float)
        matrix_me1 = makemebv(info.gmat[mal1, :], mefff)
        matrix_me2 = makemebv(info.gmat[fem1, :], mefff)
        gbv[:, 0] = (matrix_me1.sum(axis=1) + matrix_me2.sum(axis=1))/2
        gbv = pd.DataFrame(gbv)
        gbv.columns = trait_names
    elif notr > 1:
        gbv = np.zeros((no_individuals, notr+1))
        for i in range(notr):
            mefff = np.array(info.meff.iloc[:, i], float)
            matrix_me1 = makemebv(info.gmat[mal1, :], mefff)
            matrix_me2 = makemebv(info.gmat[fem1, :], mefff)
            gbv[:, i] = (matrix_me1.sum(axis=1) + matrix_me2.sum(axis=1))/2
            gbv[:, notr] = gbv[:, notr] + info.indwt[i]*gbv[:, i]
        gbv = pd.DataFrame(gbv)
        colnames = np.concatenate((trait_names, "ABV"), axis=None)
        gbv.columns = colnames
    gbv.insert(0, "FemaleIndex", fem1, True)  # insert ID
    gbv.insert(0, "MaleIndex", mal1, True)  # insert ID
    gbv.insert(0, "FemaleID", idn.iloc[:, 1], True)  # insert ID
    gbv.insert(0, "MaleID", idn.iloc[:, 0], True)  # insert ID
    return gbv


def calcprobzygsub(info, msvmsc, thresh, sub_idz):
    """Calculate the probability of breeding top individuals."""
    subindcheckzyg(info, sub_idz)
    mal1, fem1 = submsvmsc(msvmsc, sub_idz)
    gbv = calcgbvzygsub(info, sub_idz)
    trait_names = info.meff.columns  # traits names
    notr = trait_names.size
    gbvall = calcgbv(info, None)
    if notr == 1:
        probdf = np.zeros((gbv.shape[0], notr))
        ttt = np.quantile(gbvall.iloc[:, (0+2)], q=1-thresh)
        msvmsc111 = np.array(msvmsc.iloc[mal1, (0+2)]) + np.array(
            msvmsc.iloc[fem1, (0+2)])
        probdf[:, 0] = 1 - scipy.stats.norm.cdf(
            ttt, loc=gbv.iloc[:, (0+4)], scale=np.sqrt(msvmsc111))
        probdf = pd.DataFrame(probdf)
        probdf.columns = trait_names
    elif notr > 1:
        colnam = np.concatenate((trait_names, "AG"), axis=None)
        colnam = namesdf(notr+1, colnam).decode('utf-8')
        probdf = np.zeros((gbv.shape[0], notr+1))
        t_ind = np.arange(colnam.shape[0])[np.in1d(colnam, trait_names)]
        for i in range(notr):
            ttt = np.quantile(gbvall.iloc[:, (i+2)], q=1-thresh)
            msvmsc111 = np.array(msvmsc.iloc[mal1, (t_ind[i])+2]) + np.array(
                msvmsc.iloc[fem1, (t_ind[i])+2])
            probdf[:, i] = scipy.stats.norm.cdf(
                ttt, loc=gbv.iloc[:, (i+4)], scale=np.sqrt(msvmsc111))
            probdf[:, i] = np.nan_to_num(probdf[:, i])
            probdf[:, i] = 1 - probdf[:, i]
        ttt = np.quantile(gbvall.iloc[:, (notr+2)], q=1-thresh)
        msvmsc111 = np.array(msvmsc.loc[mal1, ["AG"]]) + np.array(
            msvmsc.loc[fem1, ["AG"]])
        probdf[:, notr] = scipy.stats.norm.cdf(
            ttt, loc=gbv.iloc[:, (notr+4)], scale=np.sqrt(msvmsc111.ravel()))
        probdf[:, notr] = np.nan_to_num(probdf[:, notr])
        probdf[:, notr] = 1 - probdf[:, notr]
        probdf = pd.DataFrame(probdf)  # convert matrix to dataframe
        colnames = np.concatenate((trait_names, "ABV"), axis=None)
        probdf.columns = colnames
    probdf = pd.concat([gbv.iloc[:, 0:4], probdf], axis=1)
    return probdf


def calcindexzygsub(info, msvmsc, const, sub_idz):
    """Calc index matepairs if constant is known."""
    subindcheckzyg(info, sub_idz)
    mal1, fem1 = submsvmsc(msvmsc, sub_idz)
    gbv = calcgbvzygsub(info, sub_idz)
    trait_names = info.meff.columns  # traits names
    notr = trait_names.size
    if notr == 1:
        indexdf = np.zeros((gbv.shape[0], notr))
        msvmsc111 = np.array(msvmsc.iloc[mal1, (0+2)]) + np.array(
            msvmsc.iloc[fem1, (0+2)])
        indexdf[:, 0] = gbv.iloc[:, (0+4)] + np.sqrt(msvmsc111)*const
        indexdf = pd.DataFrame(indexdf)
        indexdf.columns = trait_names
    elif notr > 1:
        colnam = np.concatenate((trait_names, "AG"), axis=None)
        colnam = namesdf(notr+1, colnam).decode('utf-8')
        indexdf = np.zeros((gbv.shape[0], notr+1))
        t_ind = np.arange(colnam.shape[0])[np.in1d(colnam, trait_names)]
        for i in range(notr):
            msvmsc111 = np.array(msvmsc.iloc[mal1, (t_ind[i])+2]) + np.array(
                msvmsc.iloc[fem1, (t_ind[i])+2])
            indexdf[:, i] = gbv.iloc[:, (i+4)] + np.sqrt(msvmsc111)*const
        msvmsc111 = np.array(msvmsc.loc[mal1, ["AG"]]) + np.array(
            msvmsc.loc[fem1, ["AG"]])
        indexdf[:, notr] = gbv.iloc[:, (notr+4)] + (
            np.sqrt(msvmsc111)*const).ravel()
        indexdf = pd.DataFrame(indexdf)
        colnames = np.concatenate((trait_names, "ABV"), axis=None)
        indexdf.columns = colnames
    indexdf = pd.concat([gbv.iloc[:, 0:4], indexdf], axis=1)  # grp
    return indexdf


def selstrat_z(selstrat, info, sub_idz, msvmsc, throrconst, selmale, selfm,
               maxmale):
    """
    Calculate selection criteria (GEBV, PBTI, or index) for zygotes.

    Parameters
    ----------
    selstrat : str
        A str containing any of GEBV, PBTI or index
    info : class object
        A class object created using the function "datacheck"
    sub_idz : pandas.DataFrame
        Index: RangeIndex (minimum of 2 rows)
        Containing ID numbers of specific individuals to be evaluated.
        The 1st and 2nd columns must be IDS of males and females, respectively.
    msvmsc : pandas.DataFrame
        DF created using the function "msvarcov_g"
    throrconst : float
        If selstrat is PBTI, a throrconst of value 0.05 sets threshold at
        top 5% of GEBV of the population. If selstrat is index, throrconst
        a constant.
    selmale : list
        list of two items. 1st item is the str coding for males as in group
        dataframe. The 2nd item is a float representing x% of males to be used
    selfm : list
        list of two items. 1st item is the str coding for females as in group
        dataframe.The 2nd item is a float representing x% of females to be used
    maxmale : integer
        maximum number of allocations for males
    Returns
    -------
    matndat : pandas.DataFrame
        Index: RangeIndex
        Columns:
        MaleID, FemaleID, MaleIndex, FemaleIndex, trait names and ABV
    Note: If selstrat is GEBV, None may be used for throrconst and msvmsc.
    If sub_idz is None and selstrat is GEBV, GEBVs will be estimated for all
    individuals. However, if selstrat is not GEBV, the chosen selection
    criterion will be estimated for all individuals in msvmsc data frame.
    """
    if len(pd.unique(info.group.iloc[:, 0])) == 1:
        sys.exit("Inds are the same group. Use 'selstrat_g' function")
    if selstrat not in ('gebv', 'GEBV', 'pbti', 'PBTI', 'index', 'INDEX'):
        sys.exit("Options must be one of 'GEBV', PBTI', or 'INDEX'")
    if selstrat in ("PBTI", "pbti", "index", "INDEX") and msvmsc is None:
        sys.exit("Provide Mendelian (co-)variance dataframe: 'msvmsc'")
    if selstrat in ("PBTI", "pbti", "index", "INDEX") and throrconst is None:
        sys.exit("Provide value for throrconst parameter")
    if sub_idz is None:
        if maxmale is None:
            sys.exit("Provide maximum allocation for males 'maxmale'")
        elif selmale is None:
            sys.exit("Provide value for propoertio of males to be selected")
        elif selfm is None:
            sys.exit("Provide value for propoertio of females to be selected")
    if sub_idz is not None:
        if selstrat in ('gebv', 'GEBV'):
            matndat = calcgbvzygsub(info, sub_idz)
        elif selstrat in ('pbti', 'PBTI'):
            matndat = calcprobzygsub(info, msvmsc, throrconst, sub_idz)
        elif selstrat in ('index', 'INDEX'):
            matndat = calcindexzygsub(info, msvmsc, throrconst, sub_idz)
    else:
        for item in [selmale[0], selfm[0]]:
            if item not in pd.unique(info.group.iloc[:, 0]).tolist():
                sys.exit("Sex name does not match group names")
        if selstrat in ('gebv', 'GEBV'):
            matndat = pot_parents(info, calcgbv(info, None), selmale, selfm)
            matndat = selsgebv(info.meff.columns.size, matndat,
                               calcgbv(info, None), maxmale)
        elif selstrat in ('pbti', 'PBTI'):
            probdf = calcprob(info, msvmsc, throrconst)
            matndat = pot_parents(info, probdf, selmale, selfm)
            matndat = selspbtizyg(info.meff.columns.size, calcgbv(info, None),
                                  matndat, msvmsc, throrconst, maxmale)
        elif selstrat in ('index', 'INDEX'):
            indexdf = calcindex(info, msvmsc, throrconst)
            matndat = pot_parents(info, indexdf, selmale, selfm)
            matndat = selsindex(info.meff.columns.size, calcgbv(info, None),
                                matndat, msvmsc, throrconst, maxmale)
    return matndat


def agggenmsezygsub(no_markers, no_individuals, slist, slist1, indwt):
    """Set up add effects mat of agg gen (zygote) subset."""
    mmfinal = np.empty((no_individuals, no_markers))
    mmfinal1 = np.empty((no_individuals, no_markers))
    for i in range(no_individuals):
        tmpmt1 = np.zeros((indwt.size, no_markers))
        tmpmt2 = np.zeros((indwt.size, no_markers))
        for trt in range(indwt.size):
            tmpmt1[trt, :] = slist[0][trt][i, :]
            tmpmt2[trt, :] = slist1[0][trt][i, :]
        mmfinal[i, :] = np.matmul(indwt.transpose(), tmpmt1)
        mmfinal1[i, :] = np.matmul(indwt.transpose(), tmpmt2)
    return mmfinal, mmfinal1


def writechrzyg(covtmpx, chrinterest, chrm, trtnam, stdsim):
    """Write matrices to file (zygote)."""
    if isinstance(chrinterest, str):
        if chrinterest == 'all':
            chrfile1 = "{}/Sim mat zygotes {} chrm {}.npy".format(
                os.getcwd(), trtnam, chrm)  # output
            np.save(chrfile1, covtmpx)
    elif chrm in chrinterest:
        chrfile1 = "{}/Sim mat zygotes {} chrm {}.npy".format(
            os.getcwd(), trtnam, chrm)  # output file
        np.save(chrfile1, covtmpx)
    if stdsim:
        if isinstance(chrinterest, str):
            if chrinterest == 'all':
                chrfilec = "{}/Stdsim zygotes {} chrm {}.npy".format(
                    os.getcwd(), trtnam, chrm)
                np.save(chrfilec, cov2corr(covtmpx))
        elif chrm in chrinterest:
            chrfilec = "{}/Stdsim mat zygotes {} chrm {}.npy".format(
                os.getcwd(), trtnam, chrm)
            np.save(chrfilec, cov2corr(covtmpx))


def mrmcalszygsub(info, sub_idz, stdsim, covmat, chrinterest, save, progress):
    """Compute similarity matrix for each chromosome (zygote) for matepairs."""
    if progress:
        progr(0, info.meff.columns.size)
    mal1, fem1, _, probn = subindcheckzyg(info, sub_idz)
    chrinterest = chr_int(chrinterest)
    slist = traitspecmatrices(info.gmat[mal1, :], info.meff)
    slist1 = traitspecmatrices(info.gmat[fem1, :], info.meff)
    for i in range(info.meff.columns.size):
        mat = np.zeros((len(mal1), len(mal1)))
        for chrm in pd.unique(info.gmap.iloc[:, 0]):
            s_ind = np.array(np.arange(0, info.gmap.shape[0], 1
                                       )[info.gmap.iloc[:, 0] == (chrm)])
            if info.gmap.shape[1]-3 == 1:
                covtmpx = abs(
                    dgmrm(slist[0][i][:, s_ind], covmat[0][chrm-1])) + abs(
                        dgmrm(slist1[0][i][:, s_ind], covmat[0][chrm-1]))
            else:
                covtmpx = abs(
                    dgmrm(slist[0][i][:, s_ind], covmat[probn[0]][chrm-1])
                    ) + abs(dgmrm(slist1[0][i][:, s_ind],
                                  covmat[probn[1]][chrm-1]))
            mat = mat + covtmpx  # sums up chrm specific covariances
            writechrzyg(covtmpx, chrinterest, chrm, info.meff.columns[i],
                        stdsim)
        if stdsim:
            if save is True:
                covxfile = "{}/Sim mat zygote {}.npy".format(
                    os.getcwd(), info.meff.columns[i])  # output file
                np.save(covxfile, cov2corr(mat))
        else:
            if save is True:
                covxfile = "{}/Sim mat zygotes {}.npy".format(
                    os.getcwd(), info.meff.columns[i])
                np.save(covxfile, mat)
        if progress:
            progr(i + 1, info.meff.columns.size)
    return mat, probn, np.arange(0, info.gmap.shape[0], 1), slist, slist1


def simmat_z(info, covmat, sub_idz, chrinterest, save=False, stdsim=False,
             progress=False):
    """
    Compute similarity matrices using zygotic approach for specific matepairs.

    Parameters
    ----------
    info : class object
        A class object created using the function "datacheck"
    sub_idz : pandas.DataFrame
        Index: RangeIndex
        Containing ID numbers of specific individuals to be evaluated.
        The 1st and 2nd columns must be IDS of males and females, respectively.
    covmat : A list of pop cov matrices created using "popcovmat" function
    chrinterest : str or list of int
        list of chromosome numbers of interest or str with "all" or "none"
    save : bool, optional; write trait-specific sim mats to file if True
    stdsim : bool, optional; print write std sim mats to file if True
    progress : bool, optional; print progress of the task if True
    Returns
    -------
    multgrpcov : list containing simimlarity matrices for each group
    """
    if len(pd.unique(info.group.iloc[:, 0])) == 1:
        sys.exit("Inds are the same group. Try 'simmat_gsub' function")
    _, _, _, probn = subindcheckzyg(info, sub_idz)
    if (info.gmap.shape[1]-3 == 1 and len(probn) > 1):
        print("The same map will be used for both sexes")
    mat, probn, snpindexxx, slist, slist1 = mrmcalszygsub(
        info, sub_idz, stdsim, covmat, chr_int(chrinterest), save, progress)
    if info.meff.columns.size == 1:
        if stdsim:
            mat = cov2corr(mat)
    elif info.meff.columns.size > 1:
        if progress:
            print('Creating similarity matrix based on aggregate genotype')
            progr(0, max(pd.unique(info.gmap.iloc[:, 0])))
        tmpmm, tmpmfm = agggenmsezygsub(
            info.gmap.shape[0], sub_idz.shape[0], slist, slist1, info.indwt)
        mat = np.zeros((sub_idz.shape[0], sub_idz.shape[0]))
        # loop over chromososomes
        for chrm in pd.unique(info.gmap.iloc[:, 0]):
            s_ind = np.array(snpindexxx[info.gmap.iloc[:, 0] == (chrm)])
            if info.gmap.shape[1]-3 == 1:
                covtmpx = abs(dgmrm(tmpmm[:, s_ind], covmat[0][chrm-1])) + abs(
                        dgmrm(tmpmfm[:, s_ind], covmat[0][chrm-1]))
            else:
                covtmpx = abs(dgmrm(
                    tmpmm[:, s_ind], covmat[probn[0]][chrm-1])) + abs(
                        dgmrm(tmpmfm[:, s_ind], covmat[probn[1]][chrm-1]))
            mat = mat + covtmpx
            if progress:
                progr(chrm, max(pd.unique(info.gmap.iloc[:, 0])))
        if stdsim:
            mat = cov2corr(mat)
    return mat
