import sys
import smtplib
import os
import inspect
import datetime

debuglevel = 6


def debug(level, msg):
    """
    print a debug message

    input:
    level : int
        error level (0=debug, 4=info, 7=warning,...10=critical)
    msg : str
        the debug message
    """
    global debuglevel

    if level >= debuglevel:
        try:
            cf = inspect.stack()[1]
            cfile = cf.filename.split('/')[-1]
            cline = cf.lineno
            cfunction = cf.function
        except:
            cfile = 'NA'
            cline = 'NA'
            cfunction = 'NA'
        print('[%s] [%d] [%s:%s:%s] %s' % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), level, cfile, cfunction, cline, msg), file=sys.stderr)


def SetDebugLevel(level):
    global debuglevel

    debuglevel = level


def getdoc(func):
    """
    return the json version of the doc for the function func
    input:
    func : function
            the function who's doc we want to return
    output:
    doc : str (html)
            the html of the doc of the function
    """
    print(func.__doc__)
    s = func.__doc__
    s = "<pre>\n%s\n</pre>" % s
    return(s)


def tolist(data):
    """
    if data is a string, convert to [data]
    if already a list, return the list
    input:
    data : str or list of str
    output:
    data : list of str
    """
    if isinstance(data, basestring):
        return [data]
    return data


def send_email(user, pwd, recipient, subject, body):
    import smtplib

    gmail_user = user
    gmail_pwd = pwd
    FROM = user
    TO = recipient if type(recipient) is list else [recipient]
    SUBJECT = subject
    TEXT = body

    # Prepare actual message
    message = """From: %s\nTo: %s\nSubject: %s\n\n%s
    """ % (FROM, ", ".join(TO), SUBJECT, TEXT)
    try:
        server = smtplib.SMTP("smtp.gmail.com", 587)
        server.ehlo()
        server.starttls()
        server.login(gmail_user, gmail_pwd)
        server.sendmail(FROM, TO, message)
        server.close()
        return ('successfully sent the mail')
    except:
        return ('failed to send mail')


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
            debug(0, 'servertype is openu')
            server_address = 'http://0.0.0.0:5001'
    elif 'SCDB_WEBSITE_TYPE' in os.environ:
        servertype = os.environ['SCDB_WEBSITE_TYPE'].lower()
        if servertype == 'local':
            debug(0, 'servertype is local')
            server_address = 'http://127.0.0.1:5001'
        elif servertype == 'main':
            debug(0, 'servertype is main')
            # server_address = 'http://amnonim.webfactional.com/scdb_main'
            server_address = 'http://api.dbbact.org'
        elif servertype == 'develop':
            debug(0, 'servertype is develop')
            server_address = 'http://amnonim.webfactional.com/scdb_develop'
        else:
            raise ValueError('unknown server type %s in SCDB_WEBSITE_TYPE' % servertype)
    else:
        # server_address = 'http://amnonim.webfactional.com/scdb_main'
        server_address = 'http://api.dbbact.org'
        debug(0, 'using default server main (use env. variable SCDB_WEBSITE_TYPE to set)')

    return server_address
