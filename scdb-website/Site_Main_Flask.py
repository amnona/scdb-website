from flask import Blueprint,request,render_template, make_response
import urllib.parse
import os
import json
import requests
from utils import debug,getdoc

Site_Main_Flask_Obj = Blueprint('Site_Main_Flask_Obj', __name__,template_folder='templates')


def get_db_address():
	'''
	Get the database address based on the environment variable SCDB_WEBSITE_TYPE
	(use export SCDB_WEBSITE_TYPE="local" / "main"(default) / "develop")

	input:

	output:
	server_address : str
		the supercooldb server web address based on the env. variable
	'''
	if 'SCDB_WEBSITE_TYPE' in os.environ:
		servertype=os.environ['SCDB_WEBSITE_TYPE'].lower()
		if servertype=='local':
			print('servertype is local')
			server_address='http://127.0.0.1:5000'
		elif servertype=='main':
			print('servertype is main')
			server_address='http://amnonim.webfactional.com/scdb_main'
			#####!!!!!!!!!!!!!!!!!!!!!!!!!!####################################################
			server_address='http://amnonim.webfactional.com/scdb_develop'
		elif servertype=='develop':
			print('servertype is develop')
			server_address='http://amnonim.webfactional.com/scdb_develop'
		else:
			raise ValueError('unknown server type %s in SCDB_WEBSITE_TYPE' % servertype)
	else:
		server_address='http://amnonim.webfactional.com/scdb_main'
		print('using default server main (use env. variable SCDB_WEBSITE_TYPE to set')
		server_address='http://amnonim.webfactional.com/scdb_develop'
		#####!!!!!!!!!!!!!!!!!!!!!!!!!!####################################################

	return server_address


scbd_server_address = get_db_address()


@Site_Main_Flask_Obj.route('/main',methods=['POST','GET'])
def main_html():
	"""
	Title: Test Html
	URL: site/main_html
	Method: GET
	"""
	httpRes=requests.get(scbd_server_address + '/stats/stats')
	# NumOntologyTerms = 0
	NumAnnotation = 0
	NumSequences = 0
	NumSequenceAnnotation = 0
	if httpRes.status_code==200:
		jsonRes = httpRes.json()
		# NumOntologyTerms = jsonRes.get("stats").get('NumOntologyTerms')
		NumAnnotation = jsonRes.get("stats").get('NumAnnotations')
		NumSequences = jsonRes.get("stats").get('NumSequences')
		NumSequenceAnnotation = jsonRes.get("stats").get('NumSeqAnnotations')

	webPage = render_template('searchpage.html',numAnnot=(str(NumAnnotation).replace('.0','')),numSeq=(str(NumSequences).replace('.0','')),numSeqAnnot=(str(NumSequenceAnnotation).replace('.0','')))
	return webPage


@Site_Main_Flask_Obj.route('/reset_password',methods=['POST','GET'])
def reset_password():
	"""
	Title: Reset password via mail
	URL: /reset password
	Method: POST
	"""
	webpage=render_template('reset_password.html')
	return webpage


@Site_Main_Flask_Obj.route('/about',methods=['POST','GET'])
def about():
	"""
	Title: About us
	URL: /about
	Method: POST
	"""
	webpage=render_template('about.html')
	return webpage


@Site_Main_Flask_Obj.route('/search_results',methods=['POST','GET'])
def search_results():
	"""
	Title: Search results page
	URL: site/search_results
	Method: POST
	"""

	if request.method=='GET':
		sequence=request.args['sequence']
	else:
		sequence = request.form['sequence']

	# if it is short, try if it is an ontology term
	if len(sequence)<80:
		return(getontologyinfo(sequence, relpath=''))

	# long, so probably a sequence
	rdata = {}
	rdata['sequence'] = sequence
	httpRes=requests.get(scbd_server_address + '/sequences/get_annotations',json=rdata)
	webPage=render_template('seqinfo.html',sequence=sequence.upper())

	if httpRes.status_code != requests.codes.ok:
		debug(6,"Error code:" + str(httpRes.status_code))
		webPage += "Failed to get annotations for the given sequence"
	else:
		webPage += draw_annotation_details(httpRes.json().get('annotations'),'')
		# jsonResponse = httpRes.json()
		# # webPage += "<table>"
		# # webPage += "<col width='10%'>"
		# # webPage += "<col width='30%'>"
		# # webPage += "<col width='60%'>"
		# # webPage += "<tr>"
		# # webPage += "<th>Expirment id</th>"
		# # webPage += "<th>Description</th>"
		# # webPage += "<th>Details</th>"
		# # webPage += "</tr>"
		# strDetails = ""
		# for dataRow in jsonResponse.get('annotations'):
		# 	webPage += "<tr>"
		# 	webPage += "<td><a href=exp_info/"+str(dataRow.get('expid','not found'))+">" + str(dataRow.get('expid','not found')) + "</a></td>"
		# 	cdesc = getannotationstrings(dataRow)
		# 	# webPage += "<td>" + str(dataRow.get('description','not found')) + "</td>"
		# 	webPage += '<td>' + cdesc + "</td>"
		# 	#webPage += "<td>" + str(dataRow) + "</td>"
		# 	strDetails = ''
		# 	for detailesRow in dataRow.get('details'):
		# 		strDetails += str(detailesRow)
		# 	webPage += "<td>" + str(strDetails) + "</td>"
		# 	webPage += "</tr>"
		# webPage += "</table>"
	webPage += "</body>"
	webPage += "</html>"
	return webPage


def getannotationstrings(cann):
	"""
	get a nice string summary of a curation

	input:
	cann : dict from /sequences/get_annotations (one from the list)
	output:
	cdesc : str
		a short summary of each annotation
	"""
	cdesc=''
	if cann['description']:
		cdesc+=cann['description']+' ('
	if cann['annotationtype']=='diffexp':
		chigh=[]
		clow=[]
		call=[]
		for cdet in cann['details']:
			if cdet[0]=='all':
				call.append(cdet[1])
				continue
			if cdet[0]=='low':
				clow.append(cdet[1])
				continue
			if cdet[0]=='high':
				chigh.append(cdet[1])
				continue
		cdesc+=' high in '
		for cval in chigh:
			cdesc+=cval+' '
		cdesc+=' compared to '
		for cval in clow:
			cdesc+=cval+' '
		cdesc+=' in '
		for cval in call:
			cdesc+=cval+' '
	elif cann['annotationtype']=='isa':
		cdesc+=' is a '
		for cdet in cann['details']:
			cdesc+='cdet,'
	elif cann['annotationtype']=='contamination':
		cdesc+='contamination'
	else:
		cdesc+=cann['annotationtype']+' '
		for cdet in cann['details']:
			cdesc=cdesc+' '+cdet[1]+','

	if len(cdesc) >= 1 and cdesc[-1] == ',':
		cdesc = cdesc[:-1]
	return cdesc


@Site_Main_Flask_Obj.route('/annotation_info/<int:annotationid>')
def getannotationinfo(annotationid):
	"""
	get the information about an annotation
	input:
	annotationid : int
		the annotationid to get the info for
	"""
	# get the experiment info for the annotation
	rdata={}
	rdata['annotationid']=annotationid
	# get the experiment annotations
	res=requests.get(get_db_address() +'/annotations/get_annotation',params=rdata)
	if res.status_code != 200:
		return('AnnotationID %d not found' % annotationid, 400)
	annotation=res.json()

	# get the experiment details
	rdata={}
	expid=annotation['expid']
	rdata['expId']=expid
	res=requests.get(scbd_server_address +'/experiments/get_details',json=rdata)
	webPage = render_template('annotationinfo.html',expid=expid,annotationid=annotationid)
	if res.status_code==200:
		for cres in res.json()['details']:
			webPage += "<tr>"
			webPage += '<td>'+cres[0]+'</td>'
			webPage += '<td>'+cres[1]+'</td><tr>'
	else:
		webPage+='Error getting experiment details'
	webPage += '</table>'
	webPage += '<h2>Annotations Details</h2>'
	webPage += draw_annotation_details([annotation],'../')

	webPage += render_template('annotationsubdetails.html')
	annotationdetails = []
	for k,v in annotation.items():
		if type(v)==list:
			annotationdetails = v
		else:
			webPage += "<tr>"
			webPage += '<td>'+str(k)+'</td>'
			webPage += '<td>'+str(v)+'</td><tr>'
	for cad in annotationdetails:
			webPage += "<tr>"
			webPage += '<td>'+str(cad[0])+'</td>'
			webPage += '<td><a href='+urllib.parse.quote('../ontology_info/'+str(cad[1]))+'>'+str(cad[1])+'</td><tr>'

	webPage += '</table>'
	webPage += '<h2>Sequences</h2>'
	webPage += '<input type="button" onclick="location.href=\'../annotation_seq_download/%d\';" value="Download fasta" />' % annotationid

	# add the ontology parent terms for the annotation
	webPage += '<h2>Ontology terms</h2>'
	res=requests.get(get_db_address() +'/annotations/get_annotation_ontology_parents',json={'annotationid':annotationid})
	if res.status_code != 200:
		debug(6,'no ontology parents found for annotationid %d' % annotationid)
		parents=[]
	else:
		parents = res.json().get('parents')
		debug(1,'found %d parent groups for annotationid %d' % (len(parents),annotationid))
	for ctype,cparents in parents.items():
		webPage += ctype + ':'
		for cparentname in cparents:
			webPage += cparentname + ', '
		webPage += '<br>'
	return webPage



@Site_Main_Flask_Obj.route('/ontology_info/<string:term>')
def getontologyinfo(term, relpath='../'):
	"""
	get the information all studies containing an ontology term (exact or as parent)
	input:
	term : str
		the ontology term to look for
	"""
	rdata={}
	rdata['term']=term
	# get the experiment annotations
	res=requests.get(get_db_address() +'/ontology/get_annotations',params=rdata)
	webPage = render_template('ontologyterminfo.html',term=term)
	webPage += '<h2>Annotations for ontology term:</h2>'
	webPage += draw_annotation_details(res.json()['annotations'], relpath)

	return webPage



@Site_Main_Flask_Obj.route('/exp_info/<int:expid>')
def getexperimentinfo(expid):
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
	rdata={}
	rdata['expId']=expid

	# get the experiment details
	res=requests.get(scbd_server_address +'/experiments/get_details',json=rdata)
	webPage = render_template('expinfo.html',expid=expid)
	if res.status_code==200:
		for cres in res.json()['details']:
			webPage += "<tr>"
			webPage += '<td>'+cres[0]+'</td>'
			webPage += '<td>'+cres[1]+'</td><tr>'
	else:
		webPage+='Error getting experiment details'
	webPage += '</table>'
	# get the experiment annotations
	res=requests.get(scbd_server_address+'/experiments/get_annotations',json=rdata)
	webPage += '<h2>Annotations for experiment:</h2>'
	webPage += draw_annotation_details(res.json()['annotations'],'../')

	return webPage


@Site_Main_Flask_Obj.route('/forgot_password_submit',methods=['POST','GET'])
def forgot_password_submit():
	"""
	this page will send the forgoten password to the user via mail
	input:
	dataid : string
		user email

	output:
	"""

	usermail = ''
	if request.method=='GET':
		usermail=request.args['useremail']
	else:
		usermail = request.form['useremail']

	json_user={'user':usermail}
	httpRes=requests.post(scbd_server_address +'/users/forgot_password',json=json_user)
	if httpRes.status_code==200:
		webpage = render_template('done_success.html')
	else:
		webpage = render_template('done_fail.html',mes='Failed to reset password',error=httpRes.text)
	return webpage


@Site_Main_Flask_Obj.route('/user_info/<int:userid>')
def getuserid(userid):
	"""
	get the information about a user
	input:
	dataid : int
		the user id

	output:
	"""
	rdata={}
	rdata['userid']=userid
	if userid < 0:
		return "Error: Invalid user"

	# get the experiment details
	httpRes=requests.post(scbd_server_address +'/users/get_user_public_information',json=rdata)
	if httpRes.status_code==200:
		userInfo = httpRes.json()
		username = userInfo.get('name','')
		name = userInfo.get('username','')
		desc = userInfo.get('description','')
		email = userInfo.get('email','-')
		webPage = render_template('userinfo.html',userid=userid,name=name,username=username,desc=desc,email=email)

		# get user annotation
		forUserId={'foruserid':userid}
		httpRes=requests.get(scbd_server_address + '/users/get_user_annotations',json=forUserId)
		if httpRes.status_code==200:
			webPage += draw_annotation_details(httpRes.json().get('userannotations'),'../')
		webPage += "</body></html>"
	else:
		webPage = "Failed to get user information"
	return webPage


def draw_annotation_details(annotations,relpath):
	'''
	create table entries for a list of annotations

	input:
	annotations : list of dict of annotation details (from REST API)

	output:
	wpart : str
		html code for the annotations table
	'''
	wpart = render_template('annotations_table.html')

	for dataRow in annotations:
		wpart += "<tr>"
		wpart += "<td><a href=" + relpath + "exp_info/"+str(dataRow.get('expid','not found'))+">" + str(dataRow.get('expid','not found')) + "</a></td>"
		wpart += "<td><a href=" + relpath + "user_info/"+str(dataRow.get('userid',-1))+">" + str(dataRow.get('username','not found')) + "</a></td>"
		cdesc = getannotationstrings(dataRow)
		# webPage += "<td>" + str(dataRow.get('description','not found')) + "</td>"
		wpart +='<td><a href=' + relpath + 'annotation_info/'+str(dataRow.get('annotationid',-1))+'>'+cdesc+'</td>'
		wpart +='<td>'+dataRow['date']+'</td>'
		rdata={}
		rdata['annotationid']=dataRow['annotationid']
		res=requests.get(scbd_server_address+'/annotations/get_sequences',json=rdata)
		if res.status_code==200:
			wpart +='<td><a href=' + relpath + 'annotation_seq_download/'+str(dataRow.get('annotationid',-1))+'>%d</td>' % len(res.json()['seqids'])
		else:
			wpart +='<td>'+'NA'+'</td>'
		wpart += "</tr>"
	wpart += "</table>"
	return wpart


@Site_Main_Flask_Obj.route('/annotation_seq_download/<int:annotationid>')
def download_sequences(annotationid):
	'''return a download of the sequences of the annotation as fasta
	'''
	rdata={}
	rdata['annotationid']=annotationid
	# get the experiment annotations
	res=requests.get(get_db_address() +'/annotations/get_full_sequences',json=rdata)
	annotation=res.json()
	seqs = annotation.get('sequences')
	if seqs is None:
		debug(6,'No sequences found')
		return('No sequences found',400)
	output=''
	for idx,cseq in enumerate(seqs):
		output += '>%d\n%s\n' % (idx,cseq)
	response = make_response(output)
	response.headers["Content-Disposition"] = "attachment; filename=annotation-%d-sequences.fa" % annotationid
	return response

