from flask import Blueprint,request,render_template
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
		elif servertype=='develop':
			print('servertype is develop')
			server_address='http://amnonim.webfactional.com/scdb_develop'
		else:
			raise ValueError('unknown server type %s in SCDB_WEBSITE_TYPE' % servertype)
	else:
		print('using default server main (use env. variable SCDB_WEBSITE_TYPE to set')
		server_address='http://amnonim.webfactional.com/scdb_main'
	return server_address


scbd_server_address = get_db_address()


@Site_Main_Flask_Obj.route('/site/test_html',methods=['POST','GET'])
def test_html():
	"""
	Title: Test Html
	URL: site/test_html
	Method: GET
	"""

	return render_template('seqinfo.html',sequence='lala3.txt')


@Site_Main_Flask_Obj.route('/main',methods=['POST','GET'])
def main_html():
	"""
	Title: Test Html
	URL: site/main_html
	Method: GET
	"""

	webPage = ""
	webPage += "<html>"
	webPage += "<title>Seqeunce Search</title>"
	webPage += "<body>"
	webPage += "<a href='/reset_password'>Reset password</a><br>"
	webPage += "<a href='/about'>About us</a>"
	webPage += "<center>"
	webPage += "<div style='border-radius: 5px; background-color: #f2f2f2; padding: 20px;'>"
	webPage += "<form action='search_results' method='post'><h1>Sequence Search</h1><br>"
	webPage += "<input value='' style='width: 100%; font-size:20px; height: 30px; margin-bottom: 20px;' type='text' name='sequence'><br>"
	webPage += "<input style='height: 40px; width: 140px; font-size:20px;' type='submit'>"
	webPage += "</form><br><a href='search_results?sequence=tacggagggtgcgagcgttaatcggaataactgggcgtaaagggcacgcaggcggtgacttaagtgaggtgtgaaagccccgggcttaacctgggaattgcatttcatactgggtcgctagagtactttagggaggggtagaattccacg'>Or click here for example</a></div>"
	webPage += "</center></body></html>"

	cfunc=test_html
	if request.method=='POST':
		return(getdoc(cfunc))

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

	# style = "<style>table {margin:40px; border-collapse: collapse; width: 100%;} th, td {text-align: left; padding: 8px;}tr:nth-child(even){background-color: #f2f2f2}th {background-color: #4CAF50;color: white; margin-top:100px;}</style>"

	# webPage = ""
	# webPage += "<html>"
	# webPage += "<title>Seqeunce Search</title>"
	# webPage += "<head>" + style + "</head>"
	# webPage += "<body>Search results for sequence:<br>" + sequence + "<br>"

	rdata = {}
	rdata['sequence'] = sequence
	httpRes=requests.get(scbd_server_address + '/sequences/get_annotations',json=rdata)
	webPage=render_template('seqinfo.html',sequence=sequence.upper())

	if httpRes.status_code != requests.codes.ok:
		debug(6,"Error code:" + str(httpRes.status_code))
		webPage += "Failed to get annotations for the given sequence"
	else:
		webPage += draw_annotation_details(httpRes.json().get('annotations'))
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
	return cdesc


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
	webPage += draw_annotation_details(res.json()['annotations'])

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
		return "Error: Invalid user";
    
	# get the experiment details
	httpRes=requests.post(scbd_server_address +'/users/get_user_public_information',json=rdata)
	if httpRes.status_code==200:
		userInfo = httpRes.json()
		username = userInfo.get('name','')
		name = userInfo.get('username','')
		desc = userInfo.get('description','')
		email = userInfo.get('email','-')
		webPage = render_template('userinfo.html',userid=userid,name=name,username=username,desc=desc,email=email)
	else:
		webPage = "Failed to get user information"
	return webPage


def draw_annotation_details(annotations):
	'''
	create table entries for a list of annotations

	input:
	annotations : list of dict of annotation details (from REST)

	output:
	wpart : str
		html code for the annotations table
	'''
	wpart = render_template('annotations_table.html')

	for dataRow in annotations:
		wpart += "<tr>"
		wpart += "<td><a href=exp_info/"+str(dataRow.get('expid','not found'))+">" + str(dataRow.get('expid','not found')) + "</a></td>"
		wpart += "<td><a href=user_info/"+str(dataRow.get('userid',-1))+">" + str(dataRow.get('username','not found')) + "</a></td>"
		cdesc = getannotationstrings(dataRow)
		# webPage += "<td>" + str(dataRow.get('description','not found')) + "</td>"
		wpart += '<td>' + cdesc + "</td>"
		#webPage += "<td>" + str(dataRow) + "</td>"
		# wpart += "</tr>"
		wpart +='<td>'+dataRow['date']+'</td>'
		rdata={}
		rdata['annotationid']=dataRow['annotationid']
		res=requests.get(scbd_server_address+'/annotations/get_sequences',json=rdata)
		if res.status_code==200:
			wpart +='<td>%d</td>' % len(res.json()['seqids'])
			# wpart +='<td>na</td>'
		else:
			# print(res.content)
			wpart +='<td>'+'NA'+'</td>'
		wpart += "</tr>"
	wpart += "</table>"
	return wpart
