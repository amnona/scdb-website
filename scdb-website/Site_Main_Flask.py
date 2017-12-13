from flask import Blueprint, request, render_template, make_response, redirect, url_for, Markup
import urllib.parse
from collections import defaultdict
from io import TextIOWrapper
import os
import json
import requests
import operator
from utils import debug, get_fasta_seqs

import enrichment

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
    if 'OPENU_FLAG' in os.environ:
            debug(1, 'servertype is openu')
            server_address = 'http://0.0.0.0:5001'
    elif 'SCDB_WEBSITE_TYPE' in os.environ:
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


@Site_Main_Flask_Obj.route('/enrichment', methods=['POST', 'GET'])
def test_enrichment():
    '''
    Redirect to the main search page
    '''
    # TODO: fix to non hard-coded
    webPage = render_template('enrichment.html')
    return webPage

@Site_Main_Flask_Obj.route('/add_data_results', methods=['POST', 'GET'])
def add_data_results():
    """
    Title: Add data processing
    URL: site/add_data_results
    Method: POST
    """
    
    webpage = "<html></br>"
    
    if 'fastaFileTb' in request.files:
        debug(1, 'Fasta file uploaded, processing it')
        file1 = request.files['fastaFileTb']
        textfile1 = TextIOWrapper(file1)
        seqs1 = get_fasta_seqs(textfile1)
        if seqs1 is None:
            webpage += "<h2>Error: Invalid fasta file</h2><br>"
            webpage += "<br><a href=\'main\'>Back to main page</a></html>"
            return webpage    
    
    #Prepare all exp data in array
    methodName = request.form.get('methodNameTb')
    if methodName is None or len(methodName.strip()) == 0:
        methodName = 'na'
    #print(">>>>>>>>>>>>>><<<<<<<<<<<<<<<<<method name" + methodName)
    
    hiddenExpName = request.form.get('hiddenExpName')
    hiddenExpValue = request.form.get('hiddenExpValue')
    hiddenOntName = request.form.get('hiddenOntName')
    hiddenOntType = request.form.get('hiddenOntType')
    hiddenOntDetType = request.form.get('hiddenOntDetType')
    
    if hiddenOntType is None or len(hiddenOntType.strip()) == 0:
        webpage += "<h2>Error: Invalid input 111</h2><br>"
        webpage += "<br><a href=\'main\'>Back to main page</a></html>"
        return webpage
    
    #in case one of the parameters is missing
    if hiddenExpName is None or len(hiddenExpName.strip()) == 0 or hiddenExpValue is None or len(hiddenExpValue.strip()) == 0 or hiddenOntName is None or len(hiddenOntName.strip()) == 0 or hiddenOntDetType is None or len(hiddenOntDetType.strip()) == 0:
        webpage += "<h2>Error: Invalid input 222</h2><br>"
        webpage += "<br><a href=\'main\'>Back to main page</a></html>"
        return webpage
    
    expDataNameArr = hiddenExpName.split(';') #split string into a list
    expDataValueArr = hiddenExpValue.split(';') #split string into a list
    ontDataNameArr = hiddenOntName.split(';') #split string into a list
    ontDataTypeArr = hiddenOntType.split(';') #split string into a list
    ontDataDetTypeArr = hiddenOntDetType.split(';') #split string into a list
    
    #####################################################      
    # Get expirement id or -1 if doesn't exist
    #####################################################      
    rdata = {}
    rdata['nameStrArr'] = expDataNameArr
    rdata['valueStrArr'] = expDataValueArr    
    httpRes = requests.get(scbd_server_address + '/experiments/get_id_by_list',json=rdata)
    if httpRes.status_code == 200:
        jsonRes = httpRes.json()
        expId = jsonRes.get("expId")
        errorCode = jsonRes.get("errorCode")
        errorText = jsonRes.get("errorText")
        if expId >= 0 : 
            webpage += "<h2>Existing expirement ID : " + str(expId) + "</h2><br>" 
        else:
            # identification appears in more than one expirement 
            if errorCode == -2 :
                webpage += "<h2>Error : More than one experiments was found</h2><br>" 
                webpage += "<br><a href=\'main\'>Back to main page</a></html>"
                return webpage
            # expirement was not found, try to find
            elif errorCode == -1:
                rdataExp = {}
                test = []
                
                for i in range(len(expDataNameArr)):
                    test.append((expDataNameArr[i],expDataValueArr[i]))
                
                rdataExp = {'expId': -1, 'private': False, 'details' : test}

                httpRes = requests.post(scbd_server_address + '/experiments/add_details',json=rdataExp)
                if httpRes.status_code == 200:
                    jsonRes = httpRes.json()
                    expId = jsonRes.get("expId")
                    webpage += "<h2>Created new expirement ID : " + str(expId) + "</h2><br>" 
                else: 
                    webpage += "<h2>Failed to get expirement</h2><br>" 
                    webpage += "<br><a href=\'main\'>Back to main page</a></html>"
                    return webpage
    else:
        webpage += "<h2>Failed to get expirement id</h2><br>" 
        webpage += "<br><a href=\'main\'>Back to main page</a></html>"
        return webpage
    #####################################################      
    
    #####################################################      
    # Add sequences if they are missing
    #####################################################      
    rdata = {}
    rdata['sequences'] = seqs1
    rdata['primer'] = 'V4'
    
    httpRes = requests.post(scbd_server_address + '/sequences/add',json=rdata)
    if httpRes.status_code == 200:
        jsonRes = httpRes.json()
        seqList = jsonRes.get("seqIds")    
        if len(seqList) != len(seqs1):
            webpage += "<h2>Error : Failed to retrieve all sequneces IDs</h2><br>" 
            webpage += "<br><a href=\'main\'>Back to main page</a></html>"
            return webpage
        webpage += "number of sequences: " + str(len(seqs1)) + "<br>"
        webpage += "number of ids: " + str(len(seqs1)) + "<br>"
    else:
        webpage += "<h2>Failed to retrieve sequneces IDs</h2><br>" 
        webpage += "<br><a href=\'main\'>Back to main page</a></html>"
        return webpage
    #####################################################      
    
    #####################################################      
    # Get ontologies list
    #####################################################      
    rdata = {}
    rdata['ontologies'] = ontDataNameArr
    
    httpRes = requests.post(scbd_server_address + '/ontology/get',json=rdata)
    if httpRes.status_code == 200:
        jsonRes = httpRes.json()
        ontList = jsonRes.get("ontIds")    
        if len(ontDataNameArr) != len(ontList):
            webpage += "<h2>Failed to get ontologies IDs</h2><br>" 
            return webpage
        webpage += "number of ontologies: " + str(len(ontList)) + "<br>"
    else:
        webpage += "<h2>Failed to retrieve ontologies IDs</h2><br>" 
        webpage += "<br><a href=\'main\'>Back to main page</a></html>"
        return webpage
    #####################################################      
    
    for i in range(len(ontList)):
        webpage += str(ontList[i]) + "<br>" 
    
    annotationListArr = []
                
    for i in range(len(ontDataNameArr)):
        annotationListArr.append((ontDataDetTypeArr[i],ontDataNameArr[i]))
    
    rannotation = {}
    rannotation['expId'] = expId
    rannotation['sequences'] = seqs1
    rannotation['region'] = 'V4'
    rannotation['annotationType'] = ontDataTypeArr[0]
    rannotation['method'] = methodName
    rannotation['agentType'] = 'DBBact website submission'
    rannotation['description'] = ''
    rannotation['annotationList'] = annotationListArr    
    #Everything is ready to add the data
    httpRes = requests.post(scbd_server_address + '/annotations/add',json=rannotation)
    if httpRes.status_code == 200:
        jsonRes = httpRes.json()
        annotId = jsonRes.get("annotationId")    
        webpage += "Added annotations with id: " + str(annotId) + "<br>"
    else:
        webpage += "<h2>Failed to add annotations</h2><br>" 
        webpage += "<br><a href=\'main\'>Back to main page</a></html>"
        return webpage
    
    webpage += "<br><a href=\'main\'>Back to main page</a></html>"
    return webpage
    
@Site_Main_Flask_Obj.route('/enrichment_results', methods=['POST', 'GET'])
def enrichment_results():
    """
    Title: Search results page
    URL: site/search_results
    Method: POST
    """
    if 'seqs1' in request.files:
        debug(1, 'Fasta file uploaded, processing it')
        file1 = request.files['seqs1']
        textfile1 = TextIOWrapper(file1)
        seqs1 = get_fasta_seqs(textfile1)
        if seqs1 is None:
            return('Error: Uploaded file1 not recognized as fasta', 400)
    if 'seqs2' in request.files:
        debug(1, 'Fasta file uploaded, processing it')
        file2 = request.files['seqs2']
        textfile2 = TextIOWrapper(file2)
        seqs2 = get_fasta_seqs(textfile2)
        if seqs2 is None:
            return('Error: Uploaded file1 not recognized as fasta', 400)
    webpage = render_template('header.html')
    # webpage = render_template('info_header.html')
    for term_type in ['term', 'annotation']:
        webpage += "<h2>%s enrichment</h2>" % term_type
        webpage += render_template('enrichment_results.html')
        err, terms, pval, odif = enrichment.enrichment(seqs1, seqs2, term_type=term_type)
        if err:
            return err
        for idx, cterm in enumerate(terms):
            if odif[idx] < 0:
                ccolor = 'red'
            else:
                ccolor = 'blue'
            webpage += '<tr><td><span style="color:%s">%s</span></td>' % (ccolor, cterm)
            webpage += '<td>%f</td>' % odif[idx]
            webpage += '<td>%f</td>' % pval[idx]
            webpage += "</tr>"
            webpage += '</span>'
        webpage += "</table>"
    return webpage


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
    if len(sequence) < 50:
        # if number assume it is a greengenes id
        if sequence.isdigit():
            webPage = sequence_annotations(sequence)
            return webPage
        # try is it an ontology term
        err, webPage = get_ontology_info(sequence)
        if not err:
            return webPage
        # or maybe a taxonomy term
        err, webPage = get_taxonomy_info(sequence)
        if not err:
            return webPage
        return error_message('Not found', 'Keyword <b>%s</b> was not found in dbBact ontology '
                             'or taxonomy.' % sequence)

    if len(sequence) < 100:
        return('Sequences must be at least 100bp long.', 400)
    webPage = sequence_annotations(sequence)
    return webPage


@Site_Main_Flask_Obj.route('/sequence_annotations/<string:sequence>')
def sequence_annotations(sequence):
    # long, so probably a sequence
    rdata = {}
    rdata['sequence'] = sequence
    httpRes = requests.get(scbd_server_address + '/sequences/get_annotations', json=rdata)
    webPage = render_template('header.html', title='dbBact sequence annotation')
    webPage += render_template('seqinfo.html', sequence=sequence.upper(), taxonomy='na')

    if httpRes.status_code != requests.codes.ok:
        debug(6, "Error code:" + str(httpRes.status_code))
        webPage += "Failed to get annotations for sequence:\n%s" % sequence
    else:
        annotations = httpRes.json().get('annotations')
        if len(annotations) == 0:
            webPage += '<br><br><h1>No annotations for sequence found in dbBact</h1>'
            webPage += '<h2>Are you using >100bp sequences?</h2>'
            webPage += 'Note dbBact is populated mostly with V4 (EMP primer) based sequences<br>'
            webPage += 'so no results could mean you are not using a EMP-V4 derived sequence<br>'
            webPage += 'These sequences usually begin with a "TACG"'
        else:
            for cannotation in annotations:
                cannotation['website_sequences'] = [0]
            annotations = sorted(annotations, key=lambda x: x.get('num_sequences', 0), reverse=False)
            term_info = get_term_info_for_annotations(annotations)
            webPage += draw_annotation_details(annotations, term_info, show_relative_freqs=True)
    webPage += render_template('footer.html')
    return webPage


def get_annotations_terms(annotations):
    '''
    Get a list of terms present in the annotations

    Parameters
    ----------
    annotations : list of annotations

    Returns
    -------
    terms : list of str
    list of terms from all annotations (each term appears once)
    '''
    terms = set()
    for cannotation in annotations:
        details = cannotation['details']
        for cdetail in details:
            terms.add(cdetail[1])
    terms = list(terms)
    return terms


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
    res = requests.get(get_db_address() + '/sequences/get_list_annotations',
                       json={'sequences': seqs})
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

    webPage = render_template('header.html')
    webPage += '<h2>Annotations for sequence list:</h2>'
    webPage += draw_annotation_details(annotations)
    webPage += render_template('footer.html')
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
    res = requests.get(get_db_address() + '/sequences/get_fast_annotations',
                       json={'sequences': seqs})
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

    webPage = render_template('header.html')
    webPage += '<h2>Annotations for sequence list relative freq:</h2>'
    # webPage += draw_annotation_details(annotations, term_info=term_info)
    webPage += draw_annotation_details(annotations, term_info=term_info, show_relative_freqs=True)
    webPage += render_template('footer.html')
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
        message = Markup('Annotation ID <b>%d</b> was not found.' % annotationid)
        return(render_template('header.html', title='Not found') +
               render_template('error.html', title='Not found',
                               message=message) +
               render_template('footer.html'), 400)
    annotation = res.json()

    # get the experiment details
    rdata = {}
    expid = annotation['expid']
    rdata['expId'] = expid
    webPage = render_template('header.html', title='Annotation %s' % annotationid)
    webPage += render_template('annotinfo.html', annotationid=annotationid)
    res = requests.get(scbd_server_address + '/experiments/get_details', json=rdata)
    if res.status_code == 200:
        webPage += draw_experiment_info(expid, res.json()['details'])
    else:
        message = Markup('Error getting experiment details.')
        return(render_template('header.html', title='Not found') +
               render_template('error.html', title='Not found',
                               message=message) +
               render_template('footer.html'), 400)
    webPage += '<h2>Annotations Details</h2>'
    webPage += draw_annotation_table([annotation])

    print(annotation)
    webPage += render_template('annotdetail.html')
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
    webPage += '<div style="margin: 20px"><blockquote style="font-size: 1em;">'
    for ctype, cparents in parents.items():
        cparents = list(set(cparents))
        webPage += '<p>%s: ' % ctype
        for cparentname in cparents:
            webPage += '<a href=' + urllib.parse.quote('../ontology_info/' + str(cparentname)) + '>' + cparentname + '</a> '
        webPage += '</p>'
    webPage += '</blockquote></div>'
    webPage += render_template('footer.html')
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
    webPage = '<div style="margin: 20px"><button class="btn btn-default" onclick="location.href=\'%s\';"><i class="glyphicon glyphicon-download-alt"></i> Download FASTA</button></div>' % url_for('.annotation_seq_download', annotationid=annotationid)
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

    webPage = render_template('header.html', title='dbBact taxonomy')
    webPage += '<h1>Summmary for ontology term: %s</h1>\n' % term
    webPage += '<h2>Annotations:</h2>'
    webPage += draw_annotation_details(annotations)
    webPage += render_template('footer.html')
    return '', webPage


@Site_Main_Flask_Obj.route('/annotations_list')
def annotations_list():
    debug(1, 'annotations_list')
    res = requests.get(get_db_address() + '/annotations/get_all_annotations')
    if res.status_code != 200:
        msg = 'error getting annotations list: %s' % res.content
        debug(6, msg)
        return msg, msg
    webPage = render_template('header.html', title='dbBact annotation list')
    webPage += '<h2>dbBact Annotation List</h2>'
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
    webPage = render_template('header.html', title='dbBact experiment List')
    webPage += render_template('explist.html')
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
    webPage += render_template('footer.html')
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

    webPage = render_template('header.html', title='dbBact ontology')
    webPage += render_template('taxinfo.html', taxonomy=taxonomy)
    webPage += draw_annotation_details(annotations)
    webPage += render_template('footer.html')
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
    webPage = render_template('header.html')
    res = requests.get(scbd_server_address + '/experiments/get_details', json={'expId': expid})
    if res.status_code == 200:
        webPage += draw_experiment_info(expid, res.json()['details'])
    else:
        message = Markup('Experiment ID <b>%s</b> was not found.' % expid)
        return(render_template('header.html', title='Not found') +
               render_template('error.html', title='Not found',
                               message=message) +
               render_template('footer.html'), 400)

    # get the experiment annotations
    res = requests.get(scbd_server_address + '/experiments/get_annotations', json={'expId': expid})
    annotations = res.json()['annotations']
    for cannotation in annotations:
        cannotation['website_sequences'] = [-1]
    annotations = sorted(annotations, key=lambda x: x.get('num_sequences', 0), reverse=True)
    webPage += '<h2>Annotations for experiment:</h2>'
    webPage += draw_annotation_details(annotations)
    webPage += render_template('footer.html')
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
        msg = 'Error encountered when getting info for annotation ID %d: %s' % (annotationid, res.content)
        debug(6, msg)
        return(render_template('header.html', title='Not found') +
               render_template('error.html', title='Not found',
                               message=msg) +
               render_template('footer.html'), 600)
    annotation = res.json()
    shortdesc = getannotationstrings(annotation)
    webPage = render_template('header.html')
    webPage += render_template('annotseqs.html', annotationid=annotationid)
    webPage += '<div style="margin: 20px;"><blockquote style="font-size: 1em;"><p>%s</p></blockquote></div>\n' % shortdesc

    webPage += '<h2>Download</h2>'
    webPage += draw_download_fasta_button(annotationid)

    # get the sequence information for the annotation
    res = requests.get(scbd_server_address + '/annotations/get_full_sequences', json={'annotationid': annotationid})
    if res.status_code != 200:
        msg = 'Error encountered when getting sequences for annotation ID %d: %s' % (annotationid, res.content)
        debug(6, msg)
        return(render_template('header.html', title='Not found') +
               render_template('error.html', title='Not found',
                               message=msg) +
               render_template('footer.html'), 600)
    sequences = res.json()['sequences']
    webPage += draw_sequences_info(sequences)

    return webPage


def draw_sequences_info(sequences):
    webPage = render_template('seqlist.html')
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

    json_user={'user':usermail}
    httpRes=requests.post(scbd_server_address +'/users/forgot_password',json=json_user)
    if httpRes.status_code==200:
        webpage = render_template('recover_form.html')
    else:
        webpage = render_template('done_fail.html', mes='Failed to reset password', error=httpRes.text)
    return webpage


@Site_Main_Flask_Obj.route('/recover_user_password',methods=['POST','GET'])
def recover_user_password():
    """
    this function will update new user password in the db
    input:
    dataid : string
        user email

    output:
    """

    usermail = ''
    if request.method=='GET':
        usermail=request.args['user']
        recoverycode=request.args['recoverycode']
        newpassword=request.args['newpassword']
    else:
        usermail = request.form['user']
        recoverycode = request.form['recoverycode']
        newpassword = request.form['newpassword']

    json_user={}
    json_user['user']=usermail
    json_user['recoverycode']=recoverycode
    json_user['newpassword']=newpassword

    httpRes=requests.post(scbd_server_address +'/users/recover_password',json=json_user)
    if httpRes.status_code==200:
        webpage = render_template('done_success.html')
    else:
        webpage = render_template('done_fail.html',mes='Failed to reset password',error=httpRes.text)
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

        webPage = render_template('header.html', title=username)
        webPage += render_template('userinfo.html', userid=userid, name=name, username=username, desc=desc, email=email)

        # get user annotation
        forUserId = {'foruserid': userid}
        httpRes = requests.get(scbd_server_address + '/users/get_user_annotations', json=forUserId)
        if httpRes.status_code == 200:
            webPage += draw_annotation_details(httpRes.json().get('userannotations'))
        webPage += render_template('footer.html')
        return webPage
    else:
        message = Markup('Failed to get user information:<br><br><blockquote><code>%s</code></blockquote>'
                         % httpRes.content)
        return(render_template('header.html', title='Not found') +
               render_template('error.html', title='Not found',
                               message=message) +
               render_template('footer.html'))


def draw_annotation_details(annotations, term_info=None, show_relative_freqs=False):
    '''
    Create table entries for a list of annotations

    Parameters
    ----------
    annotations : list of dict of annotation details (from REST API)
    term_info : dict of dict or None (optional)
        None (default) to skip relative word cloud.
        Otherwise need to have information about all ontology terms to be drawn
        dict of {term: dict} where
            term : ontology term (str)
            dict: pairs of:
                'total_annotations' : int
                'total_sequences' : int
    show_relative_freqs: bool (optional)
        False to draw absolute term abundance word cloud
        (i.e. term size based on how many times we see the term in the annotations)
        True to draw relative term abundance word cloud
        (i.e. term size based on how many times we see the term in the annotations divided by total times we see the term in the database)

    Returns
    -------
    wpart : str
        html code for the annotations table
    '''
    # The output webpage part
    wpart = ''

    # draw the wordcloud
    wpart += draw_wordcloud(annotations, term_info, show_relative_freqs=show_relative_freqs)

    wpart += render_template('tabs.html')

    # draw the annotation table
    wpart += draw_annotation_table(annotations)

    # wpart += '<div style="-webkit-column-count: 3; -moz-column-count: 3; column-count: 3;">\n'

    # draw the ontology term list
    wpart += draw_ontology_list(annotations, term_info)

    wpart += '    </div>\n'
    wpart += '  </div>\n'
    return wpart


def draw_wordcloud(annotations, term_info=None, show_relative_freqs=False):
    '''
    draw the wordcloud (image embedded in the html)

    Parameters
    ----------
    annotations : annotation
        The list of annotations to process for annotation ontology terms
    term_info: dict or None
        a dict with the total annotations per ontology term or None to skip relative abundance word cloud
    show_relative_freqs: bool (optional)
        False to draw absolute term abundance word cloud
        (i.e. term size based on how many times we see the term in the annotations)
        True to draw relative term abundance word cloud
        (i.e. term size based on how many times we see the term in the annotations divided by total times we see the term in the database)

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

    # calculate the relative enrichment of each term (if term_info is supplied)
    term_frac = None
    if term_info is not None:
        debug(1, 'drawing relative frequencies wordcloud')
        # do the relative freq. word cloud
        term_frac = {}
        for cterm in num_term:
            if cterm not in term_info:
                debug(2, 'term %s not in term_info!' % cterm)
                continue
            # if we don't have enough statistics about the term, ignore it
            # so we need at least 4 annotations with this term.
            # Otherwise we get a lot of discretization effect (i.e. 100% of the times we observe this term are fish,
            # but we have 1 fish annotation)
            if term_info[cterm]['total_annotations'] < 4:
                debug(2, 'term %s has <4 (%d) total annotations' % (cterm, term_info[cterm]['total_annotations']))
                continue
            if num_term[cterm] == 0:
                debug(4,'numterm for %s is 0' % cterm)
                continue
            # we use -2 to give lower weight to low. num
            term_frac[cterm] = num_term[cterm] / term_info[cterm]['total_annotations']
        # wordcloud_image = draw_cloud(term_frac, num_high_term=num_high_term, num_low_term=num_low_term)
        # wordcloud_image = draw_cloud(num_term, num_high_term=num_high_term, num_low_term=num_low_term, term_frac=term_frac)
        # wpart += render_template('wordcloud.html', wordcloudimage=urllib.parse.quote(wordcloud_image))

    if show_relative_freqs:
        debug(1, 'drawing relative freq. wordcloud')
        # draw relative frequencies
        if len(term_frac) == 0:
            debug(4, 'not enough info for relative freq - switching to absolute')
            wordcloud_image = draw_cloud(num_term, num_high_term=num_high_term, num_low_term=num_low_term, term_frac=term_frac)
        else:
            wordcloud_image = draw_cloud(term_frac, num_high_term=num_high_term, num_low_term=num_low_term, term_frac=term_frac)
    else:
        debug(1, 'drawing absolute count wordcloud')
        # draw absolute frequencies
        wordcloud_image = draw_cloud(num_term, num_high_term=num_high_term, num_low_term=num_low_term, term_frac=term_frac)

    wordcloudimage=urllib.parse.quote(wordcloud_image)
    if wordcloudimage:
        wpart += render_template('wordcloud.html', wordcloudimage=urllib.parse.quote(wordcloud_image))
    else:
        wpart += '<p></p>'
    return wpart


def draw_annotation_table(annotations):
    wpart = '<div id="annot-table" class="tab-pane in active" style="margin-top: 20px; margin-bottom: 20px;">\n'

    # the table header and css
    wpart += render_template('annottable.html')
    for dataRow in annotations:
        wpart += '  <tr>'
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
        wpart += '</tr>\n'
    wpart += '</table>\n'
    wpart += '</div>\n'
    return wpart


def draw_ontology_list(annotations, term_info=None):
    '''
    create table entries for a list of ontology terms

    input:
    annotations : list of dict of annotation details (from REST API)
    term_info : dict of dict or None (optional)
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
    wpart = '<div id="onto-list" class="tab-pane" style="margin-top: 20px; margin-bottom: 20px;">\n'

    # draw the ontology terms list
    common_terms = get_common_terms(annotations)
    wpart += '<table style="width: 100%;">\n'
    wpart += '<col><col width="100px">\n'
    wpart += '<tr><th>Term</th><th>No.</th></tr>\n'
    for cterm in common_terms:
        wpart += '<tr><td><a href=%s>%s</a></td><td>%d</td></tr>\n' % (url_for('.ontology_info', term=cterm[0]), cterm[0], cterm[1])
        # wpart += '<a href=' + urllib.parse.quote(relpath + 'ontology_info/' + cterm[0]) + '>%s</a>: %d <br>' % (cterm[0], cterm[1])
    wpart += '</table>\n'

    # draw the ontology term relative frequencies
    if term_info is not None:
        wpart += '<table style="width: 100%;">\n'
        for cterm, cinfo in term_info.items():
            if cinfo['total_annotations'] is None:
                debug(4, 'missing info total_annotations for %s' % cterm)
                continue
            if cinfo['total_sequences'] is None:
                debug(4, 'missing info total_sequences for %s' % cterm)
                continue
            wpart += '<tr><td>%s : %d, %d</td></tr>\n' % (cterm, cinfo['total_annotations'], cinfo['total_sequences'])
        wpart += '</table>\n'
    wpart += '</div>\n'
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


def _get_color(word, font_size, position, orientation, font_path, random_state, num_high_term=None, num_low_term=None, term_frac=None):
    '''Get the color for a given ontology term in the word cloud
    The color depends on whether the term was observed in 'higher in' or 'lower in' (or both or none)
    if only higher in - green
    if only lower in - red
    if both - yellow
    if none (i.e. all annotations with this term in our annotation set are 'all') - black

    Intensity depends on the tern_frac for this tern

    This function is used by WordCloud

    Parameters
    ----------
    word : str
        The term (word) to determine the color for
    num_high_term : dict of {str: int} (optional)
        The number of times the term (key) was observed as 'higher in'
    num_low_term : dict of {str: int} (optional)
        The number of times the term (key) was observed as 'lower in'
    term_frac : dict of {str: float} (optional)
        The fraction of times we observed this term (key) in our annotations
        out of total annotations containing the term (key) in the database

    Returns
    -------
    str of hex formatted color (i.e. '#00FFAC')
    '''
    # debug(1,**kwargs)
    clevel = hex(200)[2:]
    if term_frac is not None:
        if word in term_frac:
            clevel = format(int(155 + 100 * term_frac[word]), '02x')
    if num_high_term is None:
        return '#0000' + clevel
        # return '#0000ff'
    if num_low_term is None:
        return '#0000' + clevel
        # return '#0000ff'
    if word in num_high_term:
        if word in num_low_term:
            return '#' + clevel + clevel + '00'
            # return '#999900'
        return '#00' + clevel + '00'
        # return '#00ff00'
    else:
        if word in num_low_term:
            return '#' + clevel + '0000'
            # return '#ff0000'
        return '#0000' + clevel
        # return '#0000ff'


def draw_cloud(words, num_high_term=None, num_low_term=None, term_frac=None):
    '''
    Draw a wordcloud for a list of terms

    Parameters
    ----------
    words : str or dict
        If str, the terms (each replicated according to it's frequency).
        If dict, key is term and value is the number of times we observed this term (or relative frequency)
    num_high_term : dict of {str: int} (optional)
        The number of times the term (key) was observed as 'higher in'
    num_low_term : dict of {str: int} (optional)
        The number of times the term (key) was observed as 'lower in'
    term_frac : dict of {str: float} (optional)
        The fraction of times we observed this term (key) in our annotations
        out of total annotations containing the term (key) in the database

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

    wc = WordCloud(background_color="white", relative_scaling=0.5, stopwords=set(), color_func=lambda *x, **y: _get_color(*x, **y, num_high_term=num_high_term, num_low_term=num_low_term, term_frac=term_frac))
    if isinstance(words, str):
        debug(1, 'generating from words list')
        wordcloud = wc.generate(words)
    elif isinstance(words, dict):
        debug(1, 'generating from frequency dict')
        print(words)
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


def get_term_info_for_annotations(annotations):
    '''
    Get the statistics about each term in annotations

    Parameters
    ----------
    annotations: list of annotations

    Returns
    -------
    term_info: dict of XXX
    The statistics about each term
    '''
    terms = get_annotations_terms(annotations)
    res = requests.get(get_db_address() + '/ontology/get_term_stats', json={'terms': terms})
    if res.status_code != 200:
        debug(6, 'error encountered in get_term_stats: %s' % res.reason)
        return []
    ans = res.json()
    term_info = ans.get('term_info')
    return term_info


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


"""
Auto complete tests
"""
@Site_Main_Flask_Obj.route('/add_data', methods=['POST', 'GET'])
def add_data():
    """
    Title: About us
    URL: /about
    Method: POST
    """
    
    res = requests.get(get_db_address() + '/ontology/get_all_descriptions')
    if res.status_code != 200:
           debug(6, 'failed to get list of ontologies')
           parents = []
    else:
           import json
           list_of_ont = json.dumps(res.json())
    
    res = requests.get(get_db_address() + '/ontology/get_all_synonyms')
    if res.status_code != 200:
           debug(6, 'failed to get list of synonyms')
           parents = []
    else:
           import json
           list_of_synonym = json.dumps(res.json())
    
    webpage = render_template('add_data.html',syn_list=list_of_synonym,ont_list=list_of_ont,display='{{display}}',group='{{group}}',query='{{query}}')
    return webpage



def error_message(title, message):
    '''
    '''
    return(render_template('header.html', title=title) +
           render_template('error.html', title=title,
                           message=Markup(message)) +
           render_template('footer.html'))


"""
Auto complete tests
"""
@Site_Main_Flask_Obj.route('/auto_complete_test', methods=['POST', 'GET'])
def auto_complete_test():
    """
    Title: About us
    URL: /about
    Method: POST
    """
    
    res = requests.get(get_db_address() + '/ontology/get_all_descriptions')
    if res.status_code != 200:
           debug(6, 'failed to get list of ontologies')
           parents = []
    else:
           import json
           list_of_ont = json.dumps(res.json())
    
    res = requests.get(get_db_address() + '/ontology/get_all_synonyms')
    if res.status_code != 200:
           debug(6, 'failed to get list of synonyms')
           parents = []
    else:
           import json
           list_of_synonym = json.dumps(res.json())
    
    webpage = render_template('demo-autocomplete.html',syn_list=list_of_synonym,ont_list=list_of_ont,display='{{display}}',group='{{group}}',query='{{query}}')
    return webpage


def error_message(title, message):
    '''
    '''
    return(render_template('header.html', title=title) +
           render_template('error.html', title=title,
                           message=Markup(message)) +
           render_template('footer.html'))
