from flask import Flask
from .Site_Main_Flask import Site_Main_Flask_Obj
import os
from .utils import debug, SetDebugLevel

dbDefaultUser = "na"  # anonymos user in case the field is empty
dbDefaultPwd = ""

recentLoginUsers = []

app = Flask(__name__)
app.register_blueprint(Site_Main_Flask_Obj)


def gunicorn(debug_level=6):
    '''The entry point for running the api server through gunicorn (http://gunicorn.org/)
    to run dbbact rest server using gunicorn, use:

    gunicorn 'dbbact.Server_Main:gunicorn(debug_level=6)' -b 0.0.0.0:5001 --workers 4 --name=dbbact-rest-api


    Parameters
    ----------
    debug_level: int, optional
        The minimal level of debug messages to log (10 is max, ~5 is equivalent to warning)

    Returns
    -------
    Flask app
    '''
    SetDebugLevel(debug_level)
    app.debug = True
    debug(6, 'starting dbbact website server using gunicorn, debug_level=%d' % debug_level)
    return app


if __name__ == '__main__':
    SetDebugLevel(6)
    debug(2, 'starting server')
    if 'OPENU_FLAG' in os.environ:
        app.run(host='0.0.0.0', port=5000, use_reloader=False, threaded=True)
    else:
        app.run(use_reloader=False, threaded=True)
