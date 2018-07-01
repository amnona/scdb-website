import numpy as np
from scipy import stats


# new fdr method
def dsfdr(data, labels, transform_type='rankdata', method='meandiff',
          alpha=0.1, numperm=1000, fdr_method='dsfdr'):
    '''
    calculate the Discrete FDR for the data
    input:
    data : N x S numpy array
        each column is a sample (S total), each row an OTU (N total)
    labels : a 1d numpy array (length S)
        the labels of each sample (same order as data) with the group
        (0/1 if binary, 0-G-1 if G groups, or numeric values for correlation)
    transform_type : str or None
        transformation to apply to the data before caluculating
        the test statistic
        'rankdata' : rank transfrom each OTU reads
        'log2data' : calculate log2 for each OTU using minimal cutoff of 2
        'normdata' : normalize the data to constant sum per samples
        'binarydata' : convert to binary absence/presence
        'clrdata' : clr transformation of data (after replacing 0 with 1)
         None : no transformation to perform
    method : str or function
        the method to use for calculating test statistics:
        'meandiff' : mean(A)-mean(B) (binary)
        'mannwhitney' : mann-whitney u-test (binary)
        'kruwallis' : kruskal-wallis test (multiple groups)
        'stdmeandiff' : (mean(A)-mean(B))/(std(A)+std(B)) (binary)
        'spearman' : spearman correlation (numeric)
        'pearson' : pearson correlation (numeric)
        'nonzerospearman' : spearman correlation only non-zero entries
                            (numeric)
        'nonzeropearson' : pearson correlation only non-zero entries (numeric)
        function : use this function to calculate the test statistic
        (input is data,labels, output is array of float)
    alpha : float
        the desired FDR control level
    numperm : int
        number of permutations to perform
    fdr_method : str
        the FDR procedure to determine significant bacteria
        'dsfdr' : discrete FDR method
        'bhfdr' : Benjamini-Hochberg FDR method
        'byfdr' : Benjamini-Yekutielli FDR method
        'filterBH' : Benjamini-Hochberg FDR method with filtering
    output:
    reject : np array of bool (length N)
        True for OTUs where the null hypothesis is rejected
    tstat : np array of float (length N)
        the test statistic value for each OTU (for effect size)
    pvals : np array of float (length N)
        the p-value for each OTU
    '''

    data = data.copy()

    # transform the data
    if transform_type == 'rankdata':
        data = rankdata_transform(data)
    elif transform_type == 'log2data':
        data = log2data(data)
    elif transform_type == 'binarydata':
        data = binarydata(data)
    elif transform_type is None:
            pass
    else:
        raise ValueError('transform type %s not supported' % transform_type)

    numbact = np.shape(data)[0]

    labels = labels.copy()

    numbact = np.shape(data)[0]
    labels = labels.copy()

    if method == 'meandiff':
        # fast matrix multiplication based calculation
        method = meandiff
        tstat = method(data, labels)
        t = np.abs(tstat)
        numsamples = np.shape(data)[1]
        p = np.zeros([numsamples, numperm])
        k1 = 1 / np.sum(labels == 0)
        k2 = 1 / np.sum(labels == 1)
        for cperm in range(numperm):
            np.random.shuffle(labels)
            p[labels == 0, cperm] = k1
        p2 = np.ones(p.shape) * k2
        p2[p > 0] = 0
        mean1 = np.dot(data, p)
        mean2 = np.dot(data, p2)
        u = np.abs(mean1 - mean2)

    # fix floating point errors (important for permutation values!)
    # https://github.com/numpy/numpy/issues/8116
    for crow in range(numbact):
        closepos = np.isclose(t[crow], u[crow, :])
        u[crow, closepos] = t[crow]

    # calculate permutation p-vals
    pvals = np.zeros([numbact])  # p-value for original test statistic t
    pvals_u = np.zeros([numbact, numperm])
    # pseudo p-values for permutated test statistic u
    for crow in range(numbact):
        allstat = np.hstack([t[crow], u[crow, :]])
        stat_rank = rankdata(allstat, method='min')
        allstat = 1 - ((stat_rank - 1) / len(allstat))
        # assign ranks to t from biggest as 1
        pvals[crow] = allstat[0]
        pvals_u[crow, :] = allstat[1:]

    # calculate FDR
    if fdr_method == 'dsfdr':
        # sort unique p-values for original test statistics biggest to smallest
        pvals_unique = np.unique(pvals)
        sortp = pvals_unique[np.argsort(-pvals_unique)]

        # find a data-dependent threshold for the p-value
        foundit = False
        allfdr = []
        allt = []
        for cp in sortp:
            realnum = np.sum(pvals <= cp)
            fdr = (realnum + np.count_nonzero(
                pvals_u <= cp)) / (realnum * (numperm + 1))
            allfdr.append(fdr)
            allt.append(cp)
            if fdr <= alpha:
                realcp = cp
                foundit = True
                break

        if not foundit:
            # no good threshold was found
            reject = np.repeat([False], numbact)
            return reject, tstat, pvals

        # fill the reject null hypothesis
        reject = np.zeros(numbact, dtype=int)
        reject = (pvals <= realcp)
    else:
        raise ValueError('fdr method %s not supported' % fdr_method)

    return reject, tstat, pvals


def rankdata_transform(data):
    rdata = np.zeros(np.shape(data))
    for crow in range(np.shape(data)[0]):
        rdata[crow, :] = rankdata(data[crow, :])
    return rdata


def log2data(data):
    data[data < 2] = 2
    data = np.log2(data)
    return data


def binarydata(data):
    data[data != 0] = 1
    return data


def meandiff(data, labels):
    mean0 = np.mean(data[:, labels == 0], axis=1)
    mean1 = np.mean(data[:, labels == 1], axis=1)
    tstat = mean1 - mean0
    return tstat


def rankdata(a, method='average'):
    """
    rankdata(a, method='average')
    Assign ranks to data, dealing with ties appropriately.
    Ranks begin at 1.  The `method` argument controls how ranks are assigned
    to equal values.  See [1]_ for further discussion of ranking methods.
    Parameters
    ----------
    a : array_like
        The array of values to be ranked.  The array is first flattened.
    method : str, optional
        The method used to assign ranks to tied elements.
        The options are 'average', 'min', 'max', 'dense' and 'ordinal'.
        'average':
            The average of the ranks that would have been assigned to
            all the tied values is assigned to each value.
        'min':
            The minimum of the ranks that would have been assigned to all
            the tied values is assigned to each value.  (This is also
            referred to as "competition" ranking.)
        'max':
            The maximum of the ranks that would have been assigned to all
            the tied values is assigned to each value.
        'dense':
            Like 'min', but the rank of the next highest element is assigned
            the rank immediately after those assigned to the tied elements.
        'ordinal':
            All values are given a distinct rank, corresponding to the order
            that the values occur in `a`.
        The default is 'average'.
    Returns
    -------
    ranks : ndarray
         An array of length equal to the size of `a`, containing rank
         scores.
    References
    ----------
    .. [1] "Ranking", http://en.wikipedia.org/wiki/Ranking
    Examples
    --------
    >>> from scipy.stats import rankdata
    >>> rankdata([0, 2, 3, 2])
    array([ 1. ,  2.5,  4. ,  2.5])
    >>> rankdata([0, 2, 3, 2], method='min')
    array([ 1,  2,  4,  2])
    >>> rankdata([0, 2, 3, 2], method='max')
    array([ 1,  3,  4,  3])
    >>> rankdata([0, 2, 3, 2], method='dense')
    array([ 1,  2,  3,  2])
    >>> rankdata([0, 2, 3, 2], method='ordinal')
    array([ 1,  2,  4,  3])
    """
    if method not in ('average', 'min', 'max', 'dense', 'ordinal'):
        raise ValueError('unknown method "{0}"'.format(method))

    arr = np.ravel(np.asarray(a))
    algo = 'mergesort' if method == 'ordinal' else 'quicksort'
    sorter = np.argsort(arr, kind=algo)

    inv = np.empty(sorter.size, dtype=np.intp)
    inv[sorter] = np.arange(sorter.size, dtype=np.intp)

    if method == 'ordinal':
        return inv + 1

    arr = arr[sorter]
    obs = np.r_[True, arr[1:] != arr[:-1]]
    dense = obs.cumsum()[inv]

    if method == 'dense':
        return dense

    # cumulative counts of each unique value
    count = np.r_[np.nonzero(obs)[0], len(obs)]

    if method == 'max':
        return count[dense]

    if method == 'min':
        return count[dense - 1] + 1

    # average method
    return .5 * (count[dense] + count[dense - 1] + 1)


def rank_data(data):
    '''Replacing the scipy rank_data
    '''
    pass