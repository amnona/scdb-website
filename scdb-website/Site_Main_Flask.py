from flask import Blueprint,request,g
from flask.ext.login import login_required
from flask.ext.login import current_user
import json
import dbuser
import requests
from utils import debug,getdoc

Site_Main_Flask_Obj = Blueprint('Site_Main_Flask_Obj', __name__,template_folder='templates')


@Site_Main_Flask_Obj.route('/site/test_html',methods=['POST','GET'])
def test_html():
	"""
	Title: Test Html
	URL: site/test_html
	Method: GET
    """
    
	webPage = "<html><body><center>this is supercool db!</center></body></html>"

	cfunc=test_html
	if request.method=='POST':
		return(getdoc(cfunc))

	return webPage

@Site_Main_Flask_Obj.route('/site/main',methods=['POST','GET'])
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
	webPage += "<center>"
	webPage += "<div style='border-radius: 5px; background-color: #f2f2f2; padding: 20px;'>"
	webPage += "<form action='search_results' method='post'><h1>Sequence Search</h1><br>"
	webPage += "<input value='tacggagggtgcgagcgttaatcggaataactgggcgtaaagggcacgcaggcggtgacttaagtgaggtgtgaaagccccgggcttaacctgggaattgcatttcatactgggtcgctagagtactttagggaggggtagaattccacg' style='width: 100%; font-size:20px; height: 30px; margin-bottom: 20px;' type='text' name='sequence'><br>"
	webPage += "<input style='height: 40px; width: 140px; font-size:20px;' type='submit'>"
	webPage += "</form></div>"
	webPage += "</center></body></html>"

	cfunc=test_html
	if request.method=='POST':
	   return(getdoc(cfunc))

	return webPage

@Site_Main_Flask_Obj.route('/site/search_results',methods=['POST'])
def search_results():
	"""
	Title: Search results page
	URL: site/search_results
	Method: POST
    """

	sequence = request.form['sequence']
    
	style = "<style>table {margin:40px; border-collapse: collapse; width: 100%;} th, td {text-align: left; padding: 8px;}tr:nth-child(even){background-color: #f2f2f2}th {background-color: #4CAF50;color: white; margin-top:100px;}</style>"

	webPage = ""
	webPage += "<html>"
	webPage += "<title>Seqeunce Search</title>"
	webPage += "<head>" + style + "</head>"
	webPage += "<body>Search results for sequence:" +  sequence + "<br>"
    
	rdata = {}
	rdata['sequence'] = sequence 
	httpRes=requests.get('http://127.0.0.1:5000/sequences/get_annotations',json=rdata)
	
	if httpRes.status_code != requests.codes.ok:
	   debug(6,"Error code:" + str(httpRes.status_code))
	   webPage += "Failed to get annotations for the given sequence"
	else:
	   jsonResponse = httpRes.json()
	   webPage += "<table>"   
	   webPage += "<col width='10%'>"
	   webPage += "<col width='30%'>"
	   webPage += "<col width='60%'>"
	   webPage += "<tr>" 
	   webPage += "<th>Expirment id</th>"
	   webPage += "<th>Description</th>"
	   webPage += "<th>Details</th>"
	   webPage += "</tr>" 
	   strDetails = "" 
	   for dataRow in jsonResponse.get('annotations'):
	       webPage += "<tr>" 
	       webPage += "<td><a href="">" + str(dataRow.get('expid','not found')) + "</a></td>"
	       webPage += "<td>" + str(dataRow.get('description','not found')) + "</td>"
	       #webPage += "<td>" + str(dataRow) + "</td>"
	       strDetails = ''
	       for detailesRow in dataRow.get('details'):    
                strDetails += str(detailesRow)
	       webPage += "<td>" + str(strDetails) + "</td>"
	       webPage += "</tr>" 
	   webPage += "</table>"     
	webPage += "</body>"
	webPage += "</html>"
    
	return webPage
