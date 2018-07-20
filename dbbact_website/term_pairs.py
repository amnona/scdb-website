from collections import defaultdict
import requests

from .utils import get_db_address, debug


def get_enrichment_score(annotations, seqannotations, ignore_exp=[], term_info=None):
	'''Get f score, recall and precision for set of annotations

	Parameters
	----------
	annotations: dict of {annotationid (str): annotation(dict)}
	seqannotations: list of (seqid, [annotation ids])
	ingore_exp: list of str, optional
		list of experiment ids to ignore in the analysis
	term_info: dict of {term (str): details {"total_annotations": float, "total_sequences": float}} (see dbbact rest-api /ontology/get_term_stats) or None, optional
		The statistics about each term. if None, the function will contact dbbact to get the term_info

	Returns
	-------
	fscore: dict of {term(str): fscore(float)}
	recall: dict of {term(str): recall(float)}
	precision: dict of {term(str): precision(float)}
	term_count: dict of {term(str): total experiments(float)}
		the number of experiments where annotations for each term appear
	'''
	debug(2, 'getting enrichment scores from %d sequences' % len(seqannotations))
	debug(1, 'getting recall')
	recall = get_recall(annotations, seqannotations, ignore_exp=ignore_exp)
	debug(1, 'getting precision')
	precision = get_precision(annotations, seqannotations, ignore_exp=ignore_exp)
	debug(1, 'getting term count from get_enrichent_score()')
	term_count = get_term_total_counts(annotations, seqannotations, ignore_exp=ignore_exp)
	debug(1, 'calculating the enrichment scores')

	fscore = {}
	for cterm, crecall in recall.items():
		cprecision = precision[cterm]
		fscore[cterm] = 2 * (crecall * cprecision) / (crecall + cprecision)

	for cterm in fscore.keys():
		debug(1, 'term %s, recall %f, precision %f, fscore %f' % (cterm, recall[cterm], precision[cterm], fscore[cterm]))

	zz = sorted(fscore.items(), key=lambda x: x[1])
	debug(2, 'fscore')
	debug(2, zz[-10:])

	zz = sorted(recall.items(), key=lambda x: x[1])
	debug(2, 'recall')
	debug(2, zz[-10:])

	zz = sorted(precision.items(), key=lambda x: x[1])
	debug(2, 'precision')
	debug(2, zz[-10:])
	return fscore, recall, precision, term_count


def get_term_total_counts(annotations, seqannotations, ignore_exp=[]):
	'''Get the number of experiments containing each term from our annotations

	Used to calculate the color for the wordcloud

	Parameters
	----------
	annotations: dict of {annotationid (str): annotation(dict)}
	seqannotations: list of (seqid, [annotation ids])
	ingore_exp: list of str, optional
		list of experiment ids to ignore in the analysis

	Returns
	-------
	dict of {term(str): number of experiments(int)}
	'''
	term_exps = defaultdict(set)
	for cseqid, cseq_annotations in seqannotations:
		for cannotationid in cseq_annotations:
			cannotation = annotations[str(cannotationid)]
			cexpid = cannotation['expid']
			if cexpid in ignore_exp:
				continue
			for cterm in get_terms(cannotation):
				term_exps[cterm].add(cannotation['expid'])

	term_exp_count = {}
	for cterm, cexps in term_exps.items():
		term_exp_count[cterm] = len(cexps)
	return term_exp_count


def get_recall(annotations, seqannotations, method='exp-mean', ignore_exp=[], term_info=None):
	'''Calculate the recall (what fraction of the database enteries for this term are covered by the query)

	Parameters
	----------
	annotations: dict of {annotationid (str): annotation(dict)}
	seqannotations: list of (seqid, [annotation ids])
	method: str, optional
		the method to calculate the recall. options are:
		'exp-mean': calculate the per-experiment mean for the term
	term_info: dict of {term (str): details {"total_annotations": float, "total_sequences": float}} (see dbbact rest-api /ontology/get_term_stats) or None, optional
		The statistics about each term. if None, the function will contact dbbact to get the term_info

	Returns
	-------
	dict of {term (str): recall(float)}
	'''
	# get the term counts for all the terns
	debug(1, 'calculating recall')
	recall = defaultdict(float)
	all_terms = set()
	for cannotation in annotations.values():
		all_terms = all_terms.union(set(get_terms(cannotation)))
	all_terms_positive = [x[1:] if x[0] == '-' else x for x in all_terms]
	debug(1, 'total terms in all annotations: %d' % len(all_terms))

	if term_info is None:
		debug(2, 'term_info was None, getting from dbbact')
		term_info = get_term_info(all_terms_positive)

	num_sequences = len(seqannotations)
	debug(1, 'total sequences: %d' % num_sequences)
	for cseq, cseq_annotations in seqannotations:
		debug(1, 'processing seq %s' % cseq)
		# for each term (in all the annotations), get all the experiments where it appears
		cseq_term_exps = defaultdict(set)
		for cannotationid in cseq_annotations:
			cannotation = annotations[str(cannotationid)]
			terms = get_terms(cannotation)
			cexp = cannotation['expid']
			if cexp in ignore_exp:
				continue
			for cterm in terms:
				cseq_term_exps[cterm].add(cexp)
		# if 'whale blow' not in cseq_term_exps:
		# 	print('*whale blow not found for sequence')
		# 	print(cseq_term_exps)
		# and add the normalized count
		debug(1, 'going over exp list')
		for cterm, cexplist in cseq_term_exps.items():
			debug(1, 'processing term %s' % cterm)
			if cterm[0] != '-':
				if cterm not in term_info:
					debug(4, 'term %s not in term_info' % cterm)
					continue
				crecall = len(cexplist) / term_info[cterm]['total_experiments']
			else:
				if cterm[1:] not in term_info:
					debug(4, 'term %s not in term_info' % cterm)
					continue
				crecall = len(cexplist) / term_info[cterm[1:]]['total_experiments']
			# if cterm == 'whale blow':
			# 	print('term: %s' % cterm)
			# 	print(cexplist)
			# 	print('term info: %s' % term_info[cterm])
			recall[cterm] += crecall / num_sequences
			debug(1, 'next term')
	debug(1, 'recall contains %d terms' % len(recall))
	return recall


def get_term_info(terms):
	'''
	Get the statistics about each term in annotations

	Parameters
	----------
	terms: list of str
		the terms to get the info about

	Returns
	-------
	term_info: dict of {term (str): details {"total_annotations": float, "total_sequences": float}} (see dbbact rest-api /ontology/get_term_stats)
		The statistics about each term
	'''
	debug(2, 'getting term_info for %d terms' % len(terms))
	res = requests.get(get_db_address() + '/ontology/get_term_stats', json={'terms': terms})
	if res.status_code != 200:
		debug(6, 'error encountered in get_term_stats: %s' % res.reason)
		return []
	ans = res.json()
	term_info = ans.get('term_info')
	return term_info


def get_precision(annotations, seqannotations, method='total-annotation', ignore_exp=[]):
	'''Calculate the precision (how many of the sequences contain the term) for each term in annotations.

	Parameters
	----------
	annotations: dict of {annotationid (str): annotation(dict)}
	seqannotations: list of (seqid, [annotation ids])
	method: str, optional
		the method to calculate the precision. options are:
		'per-sequence': what fraction of the sequences contain this term at least once
		'total-annotation': what fraction of all sequences annotations contain this term (annotation can be counted more than once since iterating over all seqannotations)

	Returns
	-------
	dict of {term (str): precision(float)}
	'''
	# get the sequences where each term appears (at least once in their annotations)
	if method == 'per-sequence':
		term_seqs = defaultdict(set)
		for cseqid, cseq_annotations in seqannotations:
			for cannotationid in cseq_annotations:
				cannotation = annotations[str(cannotationid)]
				if cannotation['expid'] in ignore_exp:
					continue
				for cterm in get_terms(cannotation):
					term_seqs[cterm].add(cseqid)
		# and calculate the precision (what fraction of the sequences have this term)
		precision = {}
		total_seqs = len(seqannotations)
		for cterm, cterm_seqs in term_seqs.items():
			precision[cterm] = len(cterm_seqs) / total_seqs

	elif method == 'total-annotation':
		term_counts = defaultdict(float)
		for cseqid, cseq_annotations in seqannotations:
			cseq_term_counts = defaultdict(float)
			cseq_total_annotations = 0
			for cannotationid in cseq_annotations:
				cannotation = annotations[str(cannotationid)]
				if cannotation['expid'] in ignore_exp:
					continue
				cseq_total_annotations += 1
				for cterm in get_terms(cannotation):
					# we weigh each annotation by the number of annotations for this sequence (since we want mean over all sequences)
					cseq_term_counts[cterm] += 1
					# if we use the annotation type score - must fix normalization!!!!! need to do
					# term_counts[cterm] += get_annotation_type_score(cannotation) / cseq_total_annotations
			if cseq_total_annotations == 0:
				continue
			for cterm in cseq_term_counts.keys():
				term_counts[cterm] += cseq_term_counts[cterm] / cseq_total_annotations
		precision = {}
		total_seqs = len(seqannotations)
		for cterm, cterm_counts in term_counts.items():
			precision[cterm] = cterm_counts / total_seqs

	else:
		raise ValueError('method %s unknown' % method)

	return precision


def get_annotation_type_score(annotation):
	'''Get the score factor associated with an annotation type.
	Score is based on the annotation type (i.e. "common/highfreq/diffexp/contamination/other")

	Parameters
	----------
	annotation: dict
		as from dbbact rest-api annotations/get_annotation.
		should contain at least:
			"annotationtype"
	Returns
	-------
	float: the annotation score factor
	'''
	score = 1
	anntationtype = annotation['anntationtype']
	if anntationtype == 'highfreq':
		score = 2
	elif anntationtype == 'common':
		score = 1
	elif anntationtype == 'other':
		score = 1
	elif anntationtype == 'diffexp':
		score = 1
	elif anntationtype == 'contamination':
		score = 1
	else:
		debug(4, 'unknown annotation type %s' % anntationtype)
	return score


def get_terms(annotation):
	'''Get a list of terms present in the annotation. terms that are "lower in" are preceded by a "-"

	Parameters
	----------
	annotation: dict
		as from dbbact rest-api annotations/get_annotation.
		should contain at least:
			"annotationid" (str)
			"annotationtype"
			"details" (list of [detail_type, term])

	Returns
	-------
	list of str - the terms in the annotation
	'''
	terms = []
	details = annotation['details']
	for cdetail in details:
		cterm = cdetail[1]
		ctype = cdetail[0]
		if ctype == 'low':
			cterm = '-' + cterm
		terms.append(cterm)

	# handle the contamination annotation as well
	if annotation['annotationtype'] == 'contamination':
		terms.append('contamination')

	debug(1, 'found %d terms for annotation %s' % (len(terms), annotation['annotationid']))
	return terms


def get_enriched_term_pairs(annotations, seqannotations=None, min_exp=3, get_pairs=False, get_singles=True, sequences=None):
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
	sequences: list of str or None, optional
		the list of sequences for which the scores are calculated (used to normalize the recall)
	'''
	term_pair_score, term_pair_exps, term_precision = get_term_pairs_score(annotations, seqannotations=seqannotations, min_exp=min_exp, get_pairs=get_pairs, get_singles=get_singles)
	term_pairs = list(term_pair_score.keys())
	debug(2, 'found %d terms/term pairs' % len(term_pairs))
	res = requests.get(get_db_address() + '/ontology/get_term_pair_count', json={'term_pairs': term_pairs})
	term_count = res.json()['term_count']

	if sequences is None:
		num_sequences = len(seqannotations)
	else:
		num_sequences = len(sequences)

	new_score = {}
	for cterm_pair, cterm_score in term_pair_score.items():
		if cterm_pair in term_count:
			cscore = cterm_score / term_count[cterm_pair]
		else:
			debug(2, 'term %s not in term_count' % cterm_pair)
			cscore = 0
		new_score[cterm_pair] = cscore / num_sequences

	fscores = {}
	for cterm, cv in new_score.items():
		if cterm not in term_precision:
			debug(6, 'term %s is in new_score but not in term_total_counts' % cterm)
			continue
		fscores[cterm] = cv * term_precision[cterm] / (cv + term_precision[cterm])
		debug(2, 'term %s, recall %f, precision %f, fscore %f' % (cterm, cv, term_precision[cterm], fscores[cterm]))
	return fscores

	# # keep only the best term pair from each set of experiments
	# used_exps = set()
	# filtered_score = {}
	# sorted_score = sorted(new_score.items(), key=operator.itemgetter(1), reverse=True)
	# for cpair, cscore in sorted_score:
	# 	cexps = tuple(sorted(list(term_pair_exps[cpair])))
	# 	if cexps not in used_exps:
	# 		used_exps.add(cexps)
	# 		filtered_score[cpair] = new_score[cpair]
	# 		debug(2, 'keeping %s: %f (exps: %s)' % (cpair, cscore, cexps))
	# 	# else:
	# 	# 	print('already used: %s: %f (exps: %s)' % (cpair, cscore, cexps))
	# # return new_score
	# return filtered_score


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
		key is the term / term-pair string ("homo sapiens+feces")
		value is the score (sum over all experiments of the fraction of annotations (in the experiment) containing this term pair where it appears)
	dict of {term pair(str): exp_list(set)}
		the experiments in which each term pair appears in
	dict of {term/term_pair (str): precision (fraction of sequence annotations that contain this term) (float)}
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
	annotation_dict = {}
	for cannotation in annotations:
		annotation_dict[cannotation['annotationid']] = cannotation
	term_total_counts = defaultdict(float)
	total_annotations = 0
	if seqannotations is not None:
		for cseq, cseqannotations in seqannotations:
			for cannotation in cseqannotations:
				cterm_pairs = get_annotation_term_pairs(annotation_dict[cannotation], get_pairs=get_pairs, get_singles=get_singles)
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
		debug(1, 'getting term singles for annotation %s' % cann['annotationid'])
		for cdetail in details:
			cterm = cdetail[1]
			ctype = cdetail[0]
			if ctype == 'low':
				cterm = '-' + cterm
			term_pairs.append(cterm)
		debug(1, 'found %d single terms' % len(term_pairs))

	# add term pairs
	if get_pairs:
		debug(1, 'getting term pairs for annotation %s' % cann['annotationid'])
		current_term_len = len(term_pairs)
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
		else:
			debug(1, 'too many terms (%d). skipping term pairs' % len(details))
			# print('new details: %d' % len(details))
		debug(1, 'found %d term pairs' % (len(term_pairs) - current_term_len))
	return term_pairs
