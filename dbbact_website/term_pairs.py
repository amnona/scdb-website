from collections import defaultdict
import requests
import operator

from .utils import get_db_address, debug


def get_enriched_term_pairs(annotations, seqannotations=None, min_exp=3, get_pairs=True, get_singles=True):
	'''Get enrichment score terms including term pairs, for a given set of sequences and annotations

	Parameters
	----------
	annotations
	seqannotations: list of dict of {seqid: list of annotationids}, optional
		the annotations each sequence appears in (for the sequences to be tested for term enrichment)
	min_exp: int, optional
		the minimal number of experiments a term appears in (terms appearing in less will not be returned)
	get_pairs: bool, optional
		get enrichment for term pairs in each annotation
	get_singles: bool, optional
		get enrichment for single terms

	'''
	term_pair_score, term_pair_exps, term_total_counts = get_term_pairs_score(annotations, seqannotations=seqannotations, min_exp=min_exp, get_pairs=get_pairs, get_singles=get_singles)
	term_pairs = list(term_pair_score.keys())
	res = requests.get(get_db_address() + '/ontology/get_term_pair_count', json={'term_pairs': term_pairs})
	term_count = res.json()['term_count']

	new_score = {}
	for cterm_pair, cterm_score in term_pair_score.items():
		if cterm_pair in term_count:
			cscore = cterm_score / term_count[cterm_pair]
		else:
			cscore = 0
		new_score[cterm_pair] = cscore

	# print('----------')
	# print(term_count)
	# print('*********')
	# print(term_pair_score)

	fscores = {}
	for cterm, cv in new_score:
		if cterm not in term_total_counts:
			debug(6, 'term %s is in new_score but not in term_total_counts' % cterm)
			continue
		fscores[cterm] = cv * term_total_counts[cterm] / (cv + term_total_counts[cterm])
	return fscores

	# keep only the best term pair from each set of experiments
	used_exps = set()
	filtered_score = {}
	sorted_score = sorted(new_score.items(), key=operator.itemgetter(1), reverse=True)
	for cpair, cscore in sorted_score:
		cexps = tuple(sorted(list(term_pair_exps[cpair])))
		if cexps not in used_exps:
			used_exps.add(cexps)
			filtered_score[cpair] = new_score[cpair]
			print('keeping %s: %f (exps: %s)' % (cpair, cscore, cexps))
		# else:
		# 	print('already used: %s: %f (exps: %s)' % (cpair, cscore, cexps))
	# return new_score
	return filtered_score


def get_term_pairs_score(annotations, seqannotations=None, min_exp=3, get_pairs=True, get_singles=True):
	'''Get the term-pairs (i.e. homo spaiens+feces) score based on the annotations

	Parameters
	----------
	annotations: list of dict
		list of the dbbact annotations
	seqannotations: list of dict of {seqid: list of annotationids}, optional
		the annotations each sequence appears in (for the sequences to be tested for term enrichment)
	min_exp: int, optional
		the minimal number of experiments for the term-pair to appear in order to use it

	Returns
	-------
	dict of {str: float}
		key is the term-pair string ("homo sapiens+feces")
		value is the score (sum over all experiments of the fraction of annotations (in the experiment) containing this term pair where it appears)
	dict of {term pair(str): exp_list(set)}
		the experiments in which each term pair appears in
	term_total_counts: dict of {term: fraction of sequence annotations that contain this term}
	'''
	# group the annotations by experiment (a dict of annotationid:annotation for each experiment)
	experiments = set()
	for cann in annotations:
		experiments.add(cann['expid'])

	# get all the experiment annotations and count the number of appearances of each term pair in each annotation
	exp_term_pairs = {}
	term_pair_exps = defaultdict(set)
	for cexp in experiments:
		cexp_term_pairs = defaultdict(float)
		res = requests.get(get_db_address() + '/experiments/get_annotations', json={'expId': cexp})
		cannotations = res.json()['annotations']
		for ccann in cannotations:
			cterm_pairs = get_annotation_term_pairs(ccann, get_pairs=get_pairs, get_singles=get_singles)
			for ccterm_pair in cterm_pairs:
				# add one count to the number of annotations in the experiment where the term appears
				cexp_term_pairs[ccterm_pair] += 1
				# and add this experiment to the experiment list for the term pair
				term_pair_exps[ccterm_pair].add(cexp)
		exp_term_pairs[cexp] = cexp_term_pairs

	# count in how many experiments each term/term pair appears
	term_pair_score = defaultdict(float)
	for cann in annotations:
		cexp = cann['expid']
		cterm_pairs = get_annotation_term_pairs(cann)
		for ccterm_pair in cterm_pairs:
			# if only from one experiment - ignore
			if len(term_pair_exps[ccterm_pair]) < 1:
				continue
			# if len(term_pair_exps[ccterm_pair]) >= min_exp:
			# 	term_pair_score[ccterm_pair] += 1 / exp_term_pairs[cexp][ccterm_pair]
			term_pair_score[ccterm_pair] += 1 / (exp_term_pairs[cexp][ccterm_pair] + min_exp)

	# calculate the precision for each term (i.e. how many of the annotations for the sequences are for this term)
	term_total_counts = defaultdict(float)
	total_annotations = 0
	for cseq, cseqannotations in seqannotations:
		for cannotation in cseqannotations:
			cannotation = str(cannotation)
			cterm_pairs = get_annotation_term_pairs(annotations[cannotation], get_pairs=get_pairs, get_singles=get_singles)
			total_annotations += 1
			for cterm in cterm_pairs:
				term_total_counts[cterm] += 1
	for ck, cv in term_total_counts.items():
		term_total_counts[ck] = cv / total_annotations

	return term_pair_score, term_pair_exps, term_total_counts


def get_annotation_term_pairs(cann, max_terms=20, get_pairs=True, get_singles=True):
	'''Get the pairs of terms in the annotation and their type

	Parameters
	----------
	cann : dict
		items of the output of get_seq_annotations()

	Returns
	-------
	list of str of term1 + "+" + term2 (sorted alphabetically term1<term2)
	if term is "lower in", it will be preceeded by "-"
	'''
	term_pairs = []
	details = cann['details']
	# add single terms
	if get_singles:
		for cdetail in details:
			cterm = cdetail[1]
			ctype = cdetail[0]
			if ctype == 'low':
				cterm = '-' + cterm
			term_pairs.append(cterm)

	# add term pairs
	if get_pairs:
		if len(details) <= max_terms:
			for p1 in range(len(details)):
				# print('now detail term idx %d' % p1)
				for p2 in range(p1 + 1, len(details)):
					det1 = details[p1]
					det2 = details[p2]
					term1 = det1[1]
					term2 = det2[1]
					type1 = det1[0]
					type2 = det2[0]
					if type1 == 'low':
						term1 = '-' + term1
					if type2 == 'low':
						term2 = '-' + term2
					cnew_type = 'all'
					if type1 == type2:
						cnew_type == type1
					cnew_term = sorted([term1, term2])
					cnew_term = "+".join(cnew_term)
					# cnew_term = '%s+%s' % (term1, term2)
					term_pairs.append(cnew_term)
			# print('new details: %d' % len(details))
	return term_pairs
