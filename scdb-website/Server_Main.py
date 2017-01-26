from flask import Flask
from Site_Main_Flask import Site_Main_Flask_Obj

from utils import debug, SetDebugLevel

dbDefaultUser = "na"  # anonymos user in case the field is empty
dbDefaultPwd = ""

recentLoginUsers = []

app = Flask(__name__)
app.register_blueprint(Site_Main_Flask_Obj)

# the following function will be called for every request autentication is required

if __name__ == '__main__':
	SetDebugLevel(0)
	debug(2,'starting server')
	app.run(debug=True)
