from flask import Blueprint,request
import json
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
	webPage += "<form><h1>Sequence Search</h1><br>"
	webPage += "<input style='width: 100%; font-size:20px; height: 30px; margin-bottom: 20px;' type='text' name='sequence'><br>"
	webPage += "<input style='height: 40px; width: 140px; font-size:20px;' type='submit'>"
	webPage += "</form></div>"
	webPage += "</center></body></html>"

	cfunc=test_html
	if request.method=='POST':
	   return(getdoc(cfunc))

	return webPage
