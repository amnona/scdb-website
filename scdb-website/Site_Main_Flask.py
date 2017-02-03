from flask import Blueprint, request, render_template, make_response, redirect, url_for
import urllib.parse
from collections import defaultdict
from io import TextIOWrapper
import os
import requests
import operator
from utils import debug

Site_Main_Flask_Obj = Blueprint('Site_Main_Flask_Obj', __name__, template_folder='templates')


def get_db_address():
    '''
    Get the database address based on the environment variable SCDB_WEBSITE_TYPE
    (use export SCDB_WEBSITE_TYPE="local" / "main"(default) / "develop")

    Parameters
    ----------

    Returns
    -------
    server_address : str
        the supercooldb server web address based on the env. variable
    '''
    if 'SCDB_WEBSITE_TYPE' in os.environ:
        servertype = os.environ['SCDB_WEBSITE_TYPE'].lower()
        if servertype == 'local':
            debug(1, 'servertype is local')
            server_address = 'http://127.0.0.1:5000'
        elif servertype == 'main':
            debug(1, 'servertype is main')
            server_address = 'http://amnonim.webfactional.com/scdb_main'
        elif servertype == 'develop':
            debug(1, 'servertype is develop')
            server_address = 'http://amnonim.webfactional.com/scdb_develop'
        else:
            raise ValueError('unknown server type %s in SCDB_WEBSITE_TYPE' % servertype)
    else:
        server_address = 'http://amnonim.webfactional.com/scdb_main'
        debug(1, 'using default server main (use env. variable SCDB_WEBSITE_TYPE to set)')

    return server_address


scbd_server_address = get_db_address()


@Site_Main_Flask_Obj.route('/', methods=['POST', 'GET'])
def landing_page():
    '''
    Redirect to the main search page
    '''
    # TODO: fix to non hard-coded
    return redirect('main')


@Site_Main_Flask_Obj.route('/main', methods=['POST', 'GET'])
def main_html():
    """
    Title: the main dbBact page and search tool
    URL: site/main_html
    Method: GET
    """
    httpRes = requests.get(scbd_server_address + '/stats/stats')
    # NumOntologyTerms = 0
    NumAnnotation = 0
    NumSequences = 0
    NumSequenceAnnotation = 0
    NumExperiments = 0
    if httpRes.status_code == 200:
        jsonRes = httpRes.json()
        # NumOntologyTerms = jsonRes.get("stats").get('NumOntologyTerms')
        NumAnnotation = jsonRes.get("stats").get('NumAnnotations')
        NumSequences = jsonRes.get("stats").get('NumSequences')
        NumSequenceAnnotation = jsonRes.get("stats").get('NumSeqAnnotations')
        NumExperiments = jsonRes.get("stats").get('NumExperiments')

    webPage = render_template('searchpage.html',
                              numAnnot=(str(NumAnnotation).replace('.0', '')),
                              numSeq=(str(NumSequences).replace('.0', '')),
                              numExp=(str(NumExperiments).replace('.0', '')),
                              numSeqAnnot=(str(NumSequenceAnnotation).replace('.0', '')))
    return webPage


@Site_Main_Flask_Obj.route('/search_results', methods=['POST', 'GET'])
def search_results():
    """
    Title: Search results page
    URL: site/search_results
    Method: POST
    """

    if request.method == 'GET':
        sequence = request.args['sequence']
    else:
        sequence = request.form['sequence']

    # if we have a fasta file attached, process it
    if sequence == '':
        if 'fasta file' in request.files:
            debug(1, 'Fasta file uploaded, processing it')
            file = request.files['fasta file']
            textfile = TextIOWrapper(file)
            seqs = get_fasta_seqs(textfile)
            if seqs is None:
                return('Error: Uploaded file not recognized as fasta', 400)
            err, webpage = draw_sequences_annotations_compact(seqs)
            return webpage

    # if it is short, try if it is an ontology term
    if len(sequence) < 80:
        err, webPage = get_ontology_info(sequence)
        if not err:
            return webPage
        err, webPage = get_taxonomy_info(sequence)
        if not err:
            return webPage
        return('term %s not found in ontology or taxonomy' % sequence, 400)

    webPage = sequence_annotations(sequence)
    return webPage


@Site_Main_Flask_Obj.route('/sequence_annotations/<string:sequence>')
def sequence_annotations(sequence):
    # long, so probably a sequence
    rdata = {}
    rdata['sequence'] = sequence
    httpRes = requests.get(scbd_server_address + '/sequences/get_annotations', json=rdata)
    webPage = render_template('info_header.html', title='dbBact Sequence annotations')

    webPage += render_template('seqinfo.html', sequence=sequence.upper(), taxonomy='na')

    if httpRes.status_code != requests.codes.ok:
        debug(6, "Error code:" + str(httpRes.status_code))
        webPage += "Failed to get annotations for sequence:\n%s" % sequence
    else:
        annotations = httpRes.json().get('annotations')
        for cannotation in annotations:
            cannotation['website_sequences'] = [0]
        annotations = sorted(annotations, key=lambda x: x.get('num_sequences', 0), reverse=False)
        webPage += draw_annotation_details(annotations)
    webPage += render_template('info_end.html')
    return webPage


def draw_sequences_annotations(seqs):
    '''Draw the webpage for annotations for a list of sequences

    Parameters
    ----------
    seqs : list of str
        list of DNA sequences sequences to get annotations for

    Returns
    -------
    err : str
        the error encountered or '' if ok
    webpage : str
        the webpage for the annotations of these sequences
    '''
    res = requests.get(get_db_address() + '/sequences/get_list_annotations', json={'sequences': seqs})
    if res.status_code != 200:
        msg = 'error getting annotations for sequences : %s' % res.content
        debug(6, msg)
        return msg, msg
    seqannotations = res.json()['seqannotations']
    if len(seqannotations) == 0:
        msg = 'no sequences found'
        return msg, msg
    annotations = []
    for cseqannotation in seqannotations:
        if len(cseqannotation) == 0:
            continue
        for cannotation in cseqannotation:
            annotations.append(cannotation)

    webPage = render_template('info_header.html')
    webPage += '<h2>Annotations for sequence list:</h2>'
    webPage += draw_annotation_details(annotations)
    webPage += render_template('info_end.html')
    return '', webPage


def draw_sequences_annotations_compact(seqs):
    '''Draw the webpage for annotations for a set of sequences

    Parameters
    ----------
    seqs : list of str sequences (ACGT)

    Returns
    -------
    err : str
        the error encountered or '' if ok
    webpage : str
        the webpage for the annotations of these sequences
    '''
    res = requests.get(get_db_address() + '/sequences/get_fast_annotations', json={'sequences': seqs})
    if res.status_code != 200:
        msg = 'error getting annotations for sequences : %s' % res.content
        debug(6, msg)
        return msg, msg

    res = res.json()
    dict_annotations = res['annotations']
    seqannotations = res['seqannotations']
    if len(seqannotations) == 0:
        msg = 'no sequences found'
        return msg, msg
    term_info = res['term_info']

    # convert to dict of key=annotationid, value=list of sequences with this annotation
    annotation_seqs = defaultdict(list)
    annotation_counts = defaultdict(int)
    for cseqannotation in seqannotations:
        cseqid = cseqannotation[0]
        annotationids = cseqannotation[1]
        for cid in annotationids:
            annotation_seqs[cid].append(cseqid)
            annotation_counts[cid] += 1

    # get the sorted annotations list
    annotations = []
    sorted_annotations = sorted(annotation_counts.items(), key=operator.itemgetter(1), reverse=True)
    for csan in sorted_annotations:
        # note we need str as json dict keys are stored as strings :(
        cannotation = dict_annotations[str(csan[0])]
        cannotation['website_sequences'] = annotation_seqs[csan[0]]
        annotations.append(cannotation)

    annotations = sorted(annotations, key=lambda x: x.get('num_sequences', 0), reverse=False)
    annotations = sorted(annotations, key=lambda x: len(x.get('website_sequences', [])), reverse=True)

    webPage = render_template('info_header.html')
    # webPage = render_template('ontologyterminfo.html', term='lala')
    webPage += '<h2>Annotations for sequence list:</h2>'
    webPage += draw_annotation_details(annotations, term_info=term_info)
    webPage += render_template('info_end.html')
    return '', webPage


def getannotationstrings(cann):
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


@Site_Main_Flask_Obj.route('/annotation_info/<int:annotationid>')
def annotation_info(annotationid):
    """
    get the information about an annotation
    input:
    annotationid : int
        the annotationid to get the info for
    """
    # get the experiment info for the annotation
    rdata = {}
    rdata['annotationid'] = annotationid
    # get the experiment annotations
    res = requests.get(get_db_address() + '/annotations/get_annotation', params=rdata)
    if res.status_code != 200:
        return('AnnotationID %d not found' % annotationid, 400)
    annotation = res.json()

    # get the experiment details
    rdata = {}
    expid = annotation['expid']
    rdata['expId'] = expid
    webPage = render_template('info_header.html')
    webPage += render_template('annotationinfo.html', annotationid=annotationid)
    res = requests.get(scbd_server_address + '/experiments/get_details', json=rdata)
    if res.status_code == 200:
        webPage += draw_experiment_info(expid, res.json()['details'])
    else:
        webPage += 'Error getting experiment details'
    webPage += '<h2>Annotations Details</h2>'
    webPage += draw_annotations_table([annotation])

    print(annotation)
    webPage += render_template('annotationsubdetails.html')
    webPage += '<tr><td>%s</td><td>%s</td></tr>' % ('description', annotation['description'])
    webPage += '<tr><td>%s</td><td>%s</td></tr>' % ('type', annotation['annotationtype'])
    webPage += '<tr><td>%s</td><td>%s</td></tr>' % ('agent', annotation['agent'])
    webPage += '<tr><td>%s</td><td>%s</td></tr>' % ('method', annotation['method'])
    webPage += '<tr><td>%s</td><td>%s</td></tr>' % ('num_sequences', annotation['num_sequences'])
    webPage += '<tr><td>%s</td><td>%s</td></tr>' % ('date', annotation['date'])
    webPage += '<tr><td>%s</td><td>%s</td></tr>' % ('username', annotation['username'])
    webPage += '<tr><td>%s</td><td>%s</td></tr>' % ('private', annotation['private'])
    annotationdetails = annotation['details']

    for cad in annotationdetails:
        webPage += "<tr>"
        webPage += '<td>' + str(cad[0]) + '</td>'
        webPage += '<td><a href=' + urllib.parse.quote('../ontology_info/' + str(cad[1])) + '>' + str(cad[1]) + '</a></td></tr>'

    webPage += '</table>'

    webPage += '<h2>Sequences</h2>'
    webPage += draw_download_fasta_button(annotationid)

    # add the ontology parent terms for the annotation
    webPage += '<h2>Ontology terms</h2>'
    res = requests.get(get_db_address() + '/annotations/get_annotation_ontology_parents', json={'annotationid': annotationid})
    if res.status_code != 200:
        debug(6, 'no ontology parents found for annotationid %d' % annotationid)
        parents = []
    else:
        parents = res.json().get('parents')
        debug(1, 'found %d parent groups for annotationid %d' % (len(parents), annotationid))
    for ctype, cparents in parents.items():
        cparents = list(set(cparents))
        webPage += ctype + ':'
        for cparentname in cparents:
            webPage += '<a href=' + urllib.parse.quote('../ontology_info/' + str(cparentname)) + '>' + cparentname + '</a> '
        webPage += '<br>'
    return webPage


def draw_download_fasta_button(annotationid):
    '''
    Draw a button with a link to download the fasta sequences of the annotation

    Parameters
    ----------
    annotationid : int
        the annotationid for which to download the sequences

    Returns
    -------
    webPage : str
        html for the download button with the link to the fasta file download page
    '''
    webPage = '<input type="button" onclick="location.href=\'%s\';" value="Download fasta" />' % url_for('.annotation_seq_download', annotationid=annotationid)
    return webPage


@Site_Main_Flask_Obj.route('/ontology_info/<string:term>')
def ontology_info(term):
    """
    get the information all studies containing an ontology term (exact or as parent)
    input:
    term : str
        the ontology term to look for
    """
    err, webpage = get_ontology_info(term)
    return webpage


def get_ontology_info(term):
    """
    get the information all studies containing an ontology term (exact or as parent)
    input:
    term : str
        the ontology term to look for
    """
    # get the experiment annotations
    res = requests.get(get_db_address() + '/ontology/get_annotations', params={'term': term})
    if res.status_code != 200:
        msg = 'error getting annotations for ontology term %s: %s' % (term, res.content)
        debug(6, msg)
        return msg, msg
    annotations = res.json()['annotations']
    if len(annotations) == 0:
        debug(1, 'ontology term %s not found' % term)
        return 'term not found', 'term not found'

    for cannotation in annotations:
        cannotation['website_sequences'] = [0]

    webPage = render_template('info_header.html')
    webPage += render_template('ontologyterminfo.html', term=term)
    webPage += '<h2>Annotations:</h2>'
    webPage += draw_annotation_details(annotations)

    return '', webPage


@Site_Main_Flask_Obj.route('/annotations_list')
def annotations_list():
    debug(1, 'annotations_list')
    res = requests.get(get_db_address() + '/annotations/get_all_annotations')
    if res.status_code != 200:
        msg = 'error getting annotations list: %s' % res.content
        debug(6, msg)
        return msg, msg
    webPage = render_template('info_header.html', title='dbBact annotations list')
    webPage += '<h1>dbBact Annotations List</h1>'
    annotations = res.json()['annotations']
    for cannotation in annotations:
        cannotation['website_sequences'] = [-1]
    annotations = sorted(annotations, key=lambda x: x.get('date', 0), reverse=True)
    webPage += draw_annotation_details(annotations)
    return webPage


@Site_Main_Flask_Obj.route('/experiments_list')
def experiments_list():
    err, webpage = get_experiments_list()
    if err:
        return err, 400
    return webpage


def get_experiments_list():
    '''Get the list of experiments in the database and the details about each one
    Parameters
    ----------

    Returns
    -------
    webpage : str
        the webpage for the experiment list
    '''
    # get the experiments list
    debug(1, 'get_experiments_list')
    res = requests.get(get_db_address() + '/experiments/get_experiments_list')
    if res.status_code != 200:
        msg = 'error getting experiments list: %s' % res.content
        debug(6, msg)
        return msg, msg
    explist = res.json().get('explist', [])
    if len(explist) == 0:
        msg = 'no experiments found.'
        debug(3, msg)
        return msg, msg
    webPage = render_template('info_header.html')
    webPage += '<h1>dbBact Experiments List</h1>'
    webPage += render_template('experimentslist.html')
    for cexp in explist:
        cid = cexp[0]
        cexpname = ''
        for cdetail in cexp[1]:
            cname = cdetail[0]
            cval = cdetail[1]
            if cname != 'name':
                continue
            if len(cexpname) > 0:
                cexpname += '<br>'
            cexpname += cval
        webPage += '<tr><td><a href=%s>%d</a></td>' % (url_for('.experiment_info', expid=cid), cid)
        webPage += '<td>%s</td>' % cexpname
        # webPage += '<td><a href=exp_info/' + str(cid) + ">" + str(cid) + "</a></td>"
        # webPage += '<td>' + cval + '</td>'
        webPage += "</tr>"
    webPage += "</table>"
    webPage += render_template('info_end.html')
    return '', webPage


@Site_Main_Flask_Obj.route('/taxonomy_info/<string:taxonomy>')
def taxonomy_info(taxonomy):
    '''
    get the information all studies containing any bacteria with taxonomy as substring

    Parameters
    ----------
    taxonomy : str
        the ontology term to look for

    Returns
    -------
    err : str
        empty ('') if found, none empty if error encountered
    webPage : str
        the html of the resulting table
    '''
    err, webpage = get_taxonomy_info(taxonomy)
    return webpage


def get_taxonomy_info(taxonomy):
    '''
    get the information all studies containing any bacteria with taxonomy as substring

    Parameters
    ----------
    taxonomy : str
        the ontology term to look for

    Returns
    -------
    err : str
        empty ('') if found, none empty if error encountered
    webPage : str
        the html of the resulting table
    '''
    # get the taxonomy annotations
    res = requests.get(get_db_address() + '/sequences/get_taxonomy_annotations', json={'taxonomy': taxonomy})
    if res.status_code != 200:
        msg = 'error getting taxonomy annotations for %s: %s' % (taxonomy, res.content)
        debug(6, msg)
        return msg, msg
    annotations_counts = res.json()['annotations']
    if len(annotations_counts) == 0:
        msg = 'no annotations found for taxonomy %s' % taxonomy
        debug(1, msg)
        return msg, msg

    # convert to list of annotations with counts as a key/value
    annotations = []
    for cann in annotations_counts:
        cannotation = cann[0]
        cannotation['website_sequences'] = [-1] * cann[1]
        annotations.append(cannotation)

    annotations = sorted(annotations, key=lambda x: x.get('num_sequences', 0), reverse=False)
    annotations = sorted(annotations, key=lambda x: len(x.get('website_sequences', [])), reverse=True)

    webPage = render_template('info_header.html')
    webPage += render_template('taxonomyterminfo.html', taxonomy=taxonomy)
    webPage += draw_annotation_details(annotations)
    webPage += render_template('info_end.html')
    return '', webPage


@Site_Main_Flask_Obj.route('/exp_info/<int:expid>')
def experiment_info(expid):
    """
    get the information about a given study dataid
    input:
    dataid : int
        The dataid on the study (DataID field)

    output:
    info : list of (str,str,str)
        list of tuples for each entry in the study:
        type,value,descstring about dataid
        empty if dataid not found
    """

    # get the experiment details
    webPage = render_template('info_header.html')
    res = requests.get(scbd_server_address + '/experiments/get_details', json={'expId': expid})
    if res.status_code == 200:
        webPage += draw_experiment_info(expid, res.json()['details'])
    else:
        webPage += 'Error getting experiment details'

    # get the experiment annotations
    res = requests.get(scbd_server_address + '/experiments/get_annotations', json={'expId': expid})
    annotations = res.json()['annotations']
    for cannotation in annotations:
        cannotation['website_sequences'] = [-1]
    annotations = sorted(annotations, key=lambda x: x.get('num_sequences', 0), reverse=True)
    webPage += '<h2>Annotations for experiment:</h2>'
    webPage += draw_annotation_details(annotations)
    webPage += render_template('info_end.html')
    return webPage


def draw_experiment_info(expid, exp_details):
    '''
    Draw the table with all experiment details

    Parameters
    ----------
    expid : int
        the experiment id
    exp_details : list of (str,str)
        list of experiment detail name and type ('details' from REST API /experiments/get_details/ )

    Returns
    -------
    webPage : str
        the html of the experiment info table
    '''
    webPage = render_template('expinfo.html', expid=expid)
    for cres in exp_details:
        webPage += "<tr>"
        webPage += '<td>' + cres[0] + '</td>'
        webPage += '<td>' + cres[1] + '</td><tr>'
    webPage += '</table>'
    return webPage


@Site_Main_Flask_Obj.route('/annotation_seqs/<int:annotationid>')
def annotation_seqs(annotationid):
    '''
    get the information about all sequences in a given annotation
    input:
    annotationid : int
        The annotation for which to show the sequence info

    returns
    -------
    webPage : str
        the html page for the annotation sequences
    '''

    # get the annotation details
    res = requests.get(scbd_server_address + '/annotations/get_annotation', json={'annotationid': annotationid})
    if res.status_code != 200:
        msg = 'error encountered when getting annotation info for annotationid %d: %s' % (annotationid, res.content)
        debug(6, msg)
        return msg, 600
    annotation = res.json()
    shortdesc = getannotationstrings(annotation)
    webPage = render_template('info_header.html')
    webPage += render_template('annotationsequences.html', annotationid=annotationid)
    webPage += shortdesc

    webPage += '<h2>Download</h2>'
    webPage += draw_download_fasta_button(annotationid)

    # get the sequence information for the annotation
    res = requests.get(scbd_server_address + '/annotations/get_full_sequences', json={'annotationid': annotationid})
    if res.status_code != 200:
        msg = 'error encountered when getting annotation sequence info for annotationid %d: %s' % (annotationid, res.content)
        debug(6, msg)
        return msg, 600
    sequences = res.json()['sequences']
    webPage += draw_sequences_info(sequences)

    return webPage


def draw_sequences_info(sequences):
    webPage = render_template('sequenceslist.html')
    # sort the sequences based on taxonomy
    sequences = sorted(sequences, key=lambda x: x.get('taxonomy', ''))
    for cseqinfo in sequences:
        cseqinfo['seq'] = cseqinfo['seq'].upper()
        webPage += "<tr>"
        webPage += '<td>' + cseqinfo['taxonomy'] + '</td>'
        webPage += '<td><a href=%s>%s</a></td>' % (url_for('.sequence_annotations', sequence=cseqinfo['seq']), cseqinfo['seq'])
        webPage += '<td>' + 'na' + '</td><tr>'
    webPage += '</table>'
    return webPage


@Site_Main_Flask_Obj.route('/forgot_password_submit', methods=['POST', 'GET'])
def forgot_password_submit():
    """
    this page will send the forgoten password to the user via mail
    input:
    dataid : string
        user email

    output:
    """

    usermail = ''
    if request.method == 'GET':
        usermail = request.args['useremail']
    else:
        usermail = request.form['useremail']

    json_user = {'user': usermail}
    httpRes = requests.post(scbd_server_address + '/users/forgot_password', json=json_user)
    if httpRes.status_code == 200:
        webpage = render_template('done_success.html')
    else:
        webpage = render_template('done_fail.html', mes='Failed to reset password', error=httpRes.text)
    return webpage


@Site_Main_Flask_Obj.route('/user_info/<int:userid>')
def user_info(userid):
    """
    get the information about a user
    input:
    dataid : int
        the user id

    output:
    """
    rdata = {}
    rdata['userid'] = userid
    if userid < 0:
        return "Error: Invalid user"

    debug(1, 'get user info for user %d' % userid)
    # get the experiment details
    httpRes = requests.post(scbd_server_address + '/users/get_user_public_information', json=rdata)
    if httpRes.status_code == 200:
        userInfo = httpRes.json()
        username = userInfo.get('name', '')
        name = userInfo.get('username', '')
        desc = userInfo.get('description', '')
        email = userInfo.get('email', '-')

        webPage = render_template('info_header.html')
        webPage += render_template('userinfo.html', userid=userid, name=name, username=username, desc=desc, email=email)

        # get user annotation
        forUserId = {'foruserid': userid}
        httpRes = requests.get(scbd_server_address + '/users/get_user_annotations', json=forUserId)
        if httpRes.status_code == 200:
            webPage += draw_annotation_details(httpRes.json().get('userannotations'))
        webPage += "</body></html>"
    else:
        webPage = "Failed to get user information<br>"
        webPage += '%s' % httpRes.content
    return webPage


def draw_annotation_details(annotations, term_info=None):
    '''
    create table entries for a list of annotations

    input:
    annotations : list of dict of annotation details (from REST API)
    term_info : dict of dict or None (optiona)
        None (default) to skip relative word cloud.
        Otherwise need to have information about all ontology terms to be drawn
        dict of {term: dict} where
            term : ontology term (str)
            dict: pairs of:
                'total_annotations' : int
                'total_sequences' : int

    output:
    wpart : str
        html code for the annotations table
    '''
    # The output webpage part
    wpart = ''
    # draw the wordcloud
    wpart += draw_wordcloud(annotations, term_info)

    # draw the annotations table
    wpart += draw_annotations_table(annotations)

    # draw the ontlogy terms list
    common_terms = get_common_terms(annotations)
    for cterm in common_terms:
        wpart += '<a href=%s>%s</a>: %d<br>' % (url_for('.ontology_info', term=cterm[0]), cterm[0], cterm[1])
        # wpart += '<a href=' + urllib.parse.quote(relpath + 'ontology_info/' + cterm[0]) + '>%s</a>: %d <br>' % (cterm[0], cterm[1])

    # draw the ontology term relative frequencies
    if term_info is not None:
        for cterm, cinfo in term_info.items():
            wpart += '%s : %d, %d<br>' % (cterm, cinfo['total_annotations'], cinfo['total_sequences'])

    return wpart


def draw_wordcloud(annotations, term_info=None):
    '''
    draw the wordcloud (image embedded in the html)

    Parameters
    ----------
    annotations : annotation
        The list of annotations to process for annotation ontology terms
    term_info: dict or None
        a dict with the total annotations per ontology term or None to skip relative abundance word cloud

    Returns
    -------
    wpart : str
        an html webpage part with the wordcloud embedded
    '''
    wpart = ''

    # draw the wordcloud
    # termstr = ''
    # total frequencies of each term (not only leaves) in dbbact
    num_term = defaultdict(int)
    num_low_term = defaultdict(int)
    num_high_term = defaultdict(int)
    for cannotation in annotations:
        for cdetail in cannotation['details']:
            if cdetail[0] == 'all' or cdetail[0] == 'high':
                orig_term = cdetail[1]
                if 'website_sequences' in cannotation:
                    # if it's high freq. it's worth more
                    if cannotation['annotationtype'] == 'highfreq':
                        factor = 2
                    else:
                        factor = 1
                    num_to_add = factor * len(cannotation['website_sequences'])
                    num_high_term[orig_term] += num_to_add
                    num_term[orig_term] += num_to_add
            elif cdetail[0] == 'low':
                orig_term = cdetail[1]
                if 'website_sequences' in cannotation:
                    num_low_term[orig_term] += len(cannotation['website_sequences'])
                    num_term[orig_term] += len(cannotation['website_sequences'])
    term_frac = None
    if term_info is not None:
        debug(1, 'drawing relative frequencies wordcloud')
        # do the relative freq. word cloud
        term_frac = {}
        for cterm in num_term:
            if cterm not in term_info:
                debug(2, 'term %s not in term_info!' % cterm)
                continue
            term_frac[cterm] = num_term[cterm] / term_info[cterm]['total_annotations']
        # wordcloud_image = draw_cloud(term_frac, num_high_term=num_high_term, num_low_term=num_low_term)
        # wordcloud_image = draw_cloud(num_term, num_high_term=num_high_term, num_low_term=num_low_term, term_frac=term_frac)
        # wpart += render_template('testimg.html', wordcloudimage=urllib.parse.quote(wordcloud_image))

    # do the absolute number word cloud
    debug(1, 'drawing absolute count wordcloud')
    wordcloud_image = draw_cloud(num_term, num_high_term=num_high_term, num_low_term=num_low_term, term_frac=term_frac)
    wpart += render_template('testimg.html', wordcloudimage=urllib.parse.quote(wordcloud_image))

    return wpart


def draw_annotations_table(annotations):
    wpart = ''

    # the table header and css
    wpart += render_template('annotations_table.html')
    for dataRow in annotations:
        wpart += "<tr>"
        # add the experimentid info+link
        expid = dataRow.get('expid', 'not found')
        wpart += "<td><a href=%s>%s</a></td>" % (url_for('.experiment_info', expid=expid), expid)

        # add user name+link
        userid = dataRow.get('userid', 'not found')
        username = dataRow.get('username', 'not found')
        wpart += "<td><a href=%s>%s</a></td>" % (url_for('.user_info', userid=userid), username)

        # add the annotation description
        cdesc = getannotationstrings(dataRow)
        annotationid = dataRow.get('annotationid', -1)
        wpart += "<td><a href=%s>%s</a></td>" % (url_for('.annotation_info', annotationid=annotationid), cdesc)

        # add the annotation date
        wpart += '<td>%s</td>' % dataRow['date']

        # add the sequences
        annotationid = dataRow.get('annotationid', -1)
        num_sequences = dataRow.get('num_sequences', '?')
        if 'website_sequences' in dataRow:
            observed_sequences = len(dataRow['website_sequences'])
            sequences_string = '%s / %s' % (observed_sequences, num_sequences)
        else:
            observed_sequences = '?'
            sequences_string = '%s' % num_sequences
        wpart += "<td><a href=%s>%s</a></td>" % (url_for('.annotation_seqs', annotationid=annotationid), sequences_string)
        wpart += '</tr>'
    wpart += "</table>"
    return wpart


def draw_annotation_terms(annotations):
    wpart = ''
    common_terms = get_common_terms(annotations)
    for cterm in common_terms:
        wpart += '<a href=%s>%s</a>: %d<br' % (url_for('.ontology_info', term=cterm[0]), cterm[0], cterm[1])
    return wpart


@Site_Main_Flask_Obj.route('/annotation_seq_download/<int:annotationid>')
def annotation_seq_download(annotationid):
    '''return a download of the sequences of the annotation as fasta
    '''
    # get the experiment annotations
    res = requests.get(get_db_address() + '/annotations/get_full_sequences', json={'annotationid': annotationid})
    annotation = res.json()
    seqs = annotation.get('sequences')
    if seqs is None:
        debug(6, 'No sequences found')
        return('No sequences found', 400)
    output = ''
    for idx, cseq in enumerate(seqs):
        output += '>%d %s\n%s\n' % (idx, cseq.get('taxonomy', ''), cseq['seq'])
    response = make_response(output)
    response.headers["Content-Disposition"] = "attachment; filename=annotation-%d-sequences.fa" % annotationid
    return response


def get_common_terms(annotations):
    '''
    Get the terms most common to all the annotations

    Parameters
    ----------
    annotations : list of annotations

    Resturns
    --------
    common_terms: sorted list of (term, count)
    '''
    terms = defaultdict(int)
    for cannotation in annotations:
        for cdetail in cannotation['details']:
            if cdetail[0] == 'all' or cdetail[0] == 'high':
                if 'website_sequences' in cannotation:
                    numseqs = len(cannotation['website_sequences'])
                else:
                    numseqs = 1
                terms[cdetail[1]] += numseqs
    common_terms = []
    for k, v in terms.items():
        common_terms.append([k, v])
    common_terms = sorted(common_terms, reverse=True, key=lambda x: x[1])
    return common_terms


def get_color(word, font_size, position, orientation, font_path, random_state, num_high_term=None, num_low_term=None, term_frac=None):
    # debug(1,**kwargs)
    clevel = hex(200)[2:]
    if term_frac is not None:
        if word in term_frac:
            clevel = format(int(155 + 100 * term_frac[word]), '02x')
    if num_high_term is None:
        return '#0000' + clevel
        return '#0000ff'
    if num_low_term is None:
        return '#0000' + clevel
        return '#0000ff'
    if word in num_high_term:
        if word in num_low_term:
            return '#' + clevel + clevel + '00'
            return '#999900'
        return '#00' + clevel + '00'
        return '#00ff00'
    else:
        if word in num_low_term:
            return '#' + clevel + '0000'
            return '#ff0000'
        return '#0000' + clevel
        return '#0000ff'


def draw_cloud(words, num_high_term=None, num_low_term=None, term_frac=None):
    '''
    Draw a wordcloud for a list of terms

    Parameters
    ----------
    words : str or dict
        If str, the terms (each replicated according to it's frequency).
        If dict, key is term and value is the frequency.

    Returns
    -------
    BytesIO file with the image
    '''
    from wordcloud import WordCloud
    import matplotlib.pyplot as plt
    from io import BytesIO

    debug(1, 'draw_cloud for %d words' % len(words))
    if len(words) == 0:
        debug(2, 'no words for wordcloud')
        return ''
    # wc = WordCloud(background_color="white", width=200, height=100)

    # normalize the fractions to a scale max=1
    if term_frac is not None:
        maxval = max(term_frac.values())
        debug(1, 'maxval is %f' % maxval)
        for ckey, cval in term_frac.items():
            term_frac[ckey] = term_frac[ckey] / maxval

    wc = WordCloud(background_color="white", relative_scaling=0.5, stopwords=set(), color_func=lambda *x, **y: get_color(*x, **y, num_high_term=num_high_term, num_low_term=num_low_term, term_frac=term_frac))
    if isinstance(words, str):
        debug(1, 'generating from words list')
        wordcloud = wc.generate(words)
    elif isinstance(words, dict):
        debug(1, 'generating from frequency dict')
        wordcloud = wc.generate_from_frequencies(words)
    else:
        debug(4, 'unknown type for generate_wordcloud!')
    fig = plt.figure()
    plt.imshow(wordcloud)
    plt.axis("off")
    fig.tight_layout()
    figfile = BytesIO()
    fig.savefig(figfile, format='png', bbox_inches='tight')
    figfile.seek(0)  # rewind to beginning of file
    import base64
    figdata_png = base64.b64encode(figfile.getvalue())
    figfile.close()
    return figdata_png


def get_fasta_seqs(file):
    '''Get sequences from a fasta file

    Parameters
    ----------
    file : text file
        the text fasta file to process

    Returns
    -------
    seqs : list of str sequences (ACGT)
        the sequences in the fasta file
    '''
    debug(1, 'reading fasta file')
    seqs = []
    cseq = ''
    isfasta = False
    for cline in file:
        cline = cline.strip()
        if cline[0] == '>':
            isfasta = True
            if cseq:
                seqs.append(cseq)
            cseq = ''
        else:
            cseq += cline
    # process the last sequence
    if cseq:
        seqs.append(cseq)

    # test if we encountered '>'
    if not isfasta:
        debug(2, 'not a fasta file')
        return None

    debug(1, 'read %d sequences' % len(seqs))
    return seqs


@Site_Main_Flask_Obj.route('/reset_password', methods=['POST', 'GET'])
def reset_password():
    """
    Title: Reset password via mail
    URL: /reset password
    Method: POST
    """
    webpage = render_template('reset_password.html')
    return webpage


@Site_Main_Flask_Obj.route('/about', methods=['POST', 'GET'])
def about():
    """
    Title: About us
    URL: /about
    Method: POST
    """
    webpage = render_template('about.html')
    return webpage
