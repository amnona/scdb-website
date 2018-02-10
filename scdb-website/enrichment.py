import requests

import numpy as np
from mini_dsfdr import dsfdr
from utils import debug
from collections import defaultdict


def getannotationstrings2(cann):
    """
    get a nice string summary of a curation

    input:
    cann : dict from /sequences/get_annotations (one from the list)
    output:
    cdesc : str
        a short summary of each annotation
    """
    cdesc = ''
    if cann['description']:
        cdesc += cann['description'] + ' ('
    if cann['annotationtype'] == 'diffexp':
        chigh = []
        clow = []
        call = []
        for cdet in cann['details']:
            if cdet[0] == 'all':
                call.append(cdet[1])
                continue
            if cdet[0] == 'low':
                clow.append(cdet[1])
                continue
            if cdet[0] == 'high':
                chigh.append(cdet[1])
                continue
        cdesc += ' high in '
        for cval in chigh:
            cdesc += cval + ' '
        cdesc += ' compared to '
        for cval in clow:
            cdesc += cval + ' '
        cdesc += ' in '
        for cval in call:
            cdesc += cval + ' '
    elif cann['annotationtype'] == 'isa':
        cdesc += ' is a '
        for cdet in cann['details']:
            cdesc += 'cdet,'
    elif cann['annotationtype'] == 'contamination':
        cdesc += 'contamination'
    else:
        cdesc += cann['annotationtype'] + ' '
        for cdet in cann['details']:
            cdesc = cdesc + ' ' + cdet[1] + ','

    if len(cdesc) >= 1 and cdesc[-1] == ',':
        cdesc = cdesc[:-1]
    return cdesc


def get_seq_annotations_fast(sequences):
    import os
    rdata = {}
    rdata['sequences'] = sequences
    if 'OPENU_FLAG' in os.environ:
        res = requests.get('http://0.0.0.0:5001/sequences/get_fast_annotations', json=rdata)
    else:
        res = requests.get('http://dbbact.org/REST-API/sequences/get_fast_annotations', json=rdata)
    if res.status_code != 200:
        debug(5, 'error getting fast annotations for sequence list')
        return None, None, None
    res = res.json()

    sequence_terms = {}
    sequence_annotations = {}
    for cseq in sequences:
        sequence_terms[cseq] = []
        sequence_annotations[cseq] = []
    for cseqannotation in res['seqannotations']:
        cpos = cseqannotation[0]
        # need str since json dict is always string
        cseq = sequences[cpos]
        sequence_annotations[cseq].extend(cseqannotation[1])
        for cannotation in cseqannotation[1]:
            for k, v in res['annotations'][str(cannotation)]['parents'].items():
                if k == 'high' or k == 'all':
                    for cterm in v:
                        sequence_terms[cseq].append(cterm)
                elif k == 'low':
                    for cterm in v:
                        sequence_terms[cseq].append('-' + cterm)

    annotations = res['annotations']
    # replace the string in the key with an int (since in json key is always str)
    keys = list(annotations.keys())
    for cid in keys:
        annotations[int(cid)] = annotations.pop(cid)

    total_annotations = 0
    for cseq_annotations in sequence_annotations.values():
        total_annotations += len(cseq_annotations)
    debug(2, 'Got %d annotations' % total_annotations)
    return sequence_terms, sequence_annotations, res['annotations']


def _get_term_features(features, feature_terms):
    '''Get dict of number of appearances in each sequence keyed by term

    Parameters
    ----------
    features : list of str
        A list of DNA sequences
    feature_terms : dict of {feature: list of tuples of (term, amount)}
        The terms associated with each feature in exp
        feature (key) : str the feature (out of exp) to which the terms relate
        feature_terms (value) : list of tuples of (str or int the terms associated with this feature, count)

    Returns
    -------
    numpy array of T (terms) * F (features)
        total counts of each term (row) in each feature (column)
    list of str
        list of the terms corresponding to the numpy array rows
    '''
    # get all terms
    terms = {}
    cpos = 0
    for cfeature, ctermlist in feature_terms.items():
        for cterm, ccount in ctermlist:
            if cterm not in terms:
                terms[cterm] = cpos
                cpos += 1

    tot_features_inflated = 0
    feature_pos = {}
    for cfeature in features:
        ctermlist = feature_terms[cfeature]
        feature_pos[cfeature] = tot_features_inflated
        tot_features_inflated += len(ctermlist)

    res = np.zeros([len(terms), tot_features_inflated])

    for cfeature in features:
        for cterm, ctermcount in feature_terms[cfeature]:
            res[terms[cterm], feature_pos[cfeature]] += ctermcount
    term_list = sorted(terms, key=terms.get)
    return res, term_list


def _get_all_annotation_string_counts(features, sequence_annotations, annotations):
    feature_annotations = {}
    for cseq, annotations_list in sequence_annotations.items():
        if cseq not in features:
            continue
        newdesc = []
        for cannotation in annotations_list:
            cdesc = getannotationstrings2(annotations[cannotation])
            newdesc.append((cdesc, 1))
        feature_annotations[cseq] = newdesc
    return feature_annotations


def _get_all_term_counts(features, feature_annotations, annotations):
    feature_terms = {}
    for cfeature in features:
        annotation_list = [annotations[x] for x in feature_annotations[cfeature]]
        feature_terms[cfeature] = get_annotation_term_counts(annotation_list)
    return feature_terms


def get_annotation_term_counts(annotations):
    '''Get the annotation type corrected count for all terms in annotations

    Parameters
    ----------
    annotations : list of dict
        list of annotations

    Returns
    -------
        list of tuples (term, count)
    '''
    term_count = defaultdict(int)
    for cannotation in annotations:
        if cannotation['annotationtype'] == 'common':
            for cdesc in cannotation['details']:
                term_count[cdesc[1]] += 1
            continue
        if cannotation['annotationtype'] == 'highfreq':
            for cdesc in cannotation['details']:
                term_count[cdesc[1]] += 2
            continue
        if cannotation['annotationtype'] == 'other':
            for cdesc in cannotation['details']:
                term_count[cdesc[1]] += 0.5
            continue
        if cannotation['annotationtype'] == 'contamination':
            term_count['contamination'] += 1
            continue
        if cannotation['annotationtype'] == 'diffexp':
            for cdesc in cannotation['details']:
                if cdesc[0] == 'all':
                    term_count[cdesc[1]] += 1
                    continue
                if cdesc[0] == 'high':
                    term_count[cdesc[1]] += 2
                    continue
                if cdesc[0] == 'low':
                    term_count[cdesc[1]] -= 2
                    continue
                debug(4, 'unknown detail type %s encountered' % cdesc[0])
            continue
        if cannotation['annotationtype'] == 'other':
            continue
        debug(4, 'unknown annotation type %s encountered' % cannotation['annotationtype'])
    res = []
    for k, v in term_count.items():
        # flip and add '-' to term if negative
        if v < 0:
            k = '-' + k
            v = -v
        res.append((k, v))
    return res


def enrichment(seqs1, seqs2, term_type="term"):
    '''
    Do dbbact term and annotation enrichment analysis for 2 fasta files

    Parameters
    ----------
    fasta1 : str
        name of first fasta file
    fasta2 : str
        name of second fasta file (can also contain sequences from fasta1 - it is the background file)
    term_type : str (optional)
        type of the term to analyze for enrichment. can be:
        "term" : analyze the terms per annotation (not including parent terms)
        "annotation" : analyze the annotations associated with each sequence


    Returns
    -------
    err : str
        empty if ok, otherwise the error encountered
    term_list : list of str
        the terms which are enriched
    pvals : list of float
        the p-value for each term
    odif : list of float
        the effect size for each term
    '''
    # set the same seed (since we use a random permutation test)
    np.random.seed(2018)

    all_seqs = set(seqs1).union(set(seqs2))
    seqs2 = list(all_seqs - set(seqs1))
    if len(seqs2) == 0:
        return 'No sequences remaining in background fasta after removing the sequences of interest', None, None, None
    all_seqs = list(all_seqs)

    # get the annotations for the sequences
    info = {}
    info['sequence_terms'], info['sequence_annotations'], info['annotations'] = get_seq_annotations_fast(all_seqs)

    if term_type == 'term':
        feature_terms = _get_all_term_counts(all_seqs, info['sequence_annotations'], info['annotations'])
    elif term_type == 'annotation':
        feature_terms = _get_all_annotation_string_counts(all_seqs, info['sequence_annotations'], info['annotations'])
    else:
        debug(8, 'strange term_type encountered: %s' % term_type)
    feature_array, term_list = _get_term_features(seqs1, feature_terms)
    bg_array, term_list = _get_term_features(seqs2, feature_terms)

    all_feature_array = np.hstack([feature_array, bg_array])

    labels = np.zeros(all_feature_array.shape[1])
    labels[:feature_array.shape[1]] = 1

    keep, odif, pvals = dsfdr(all_feature_array, labels, method='meandiff', transform_type=None, alpha=0.1, numperm=1000, fdr_method='dsfdr')
    keep = np.where(keep)[0]
    if len(keep) == 0:
        debug(2, 'no enriched terms found')
    term_list = np.array(term_list)[keep]
    odif = odif[keep]
    pvals = pvals[keep]
    si = np.argsort(odif)
    odif = odif[si]
    pvals = pvals[si]
    term_list = term_list[si]
    return '', term_list, pvals, odif
