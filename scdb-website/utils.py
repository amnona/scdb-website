import sys
import smtplib


def debug(level, msg):
    """
    print a debug message
    input:
    level : int
            error level (0=debug...10=critical)
    msg : str
            the debug message
    """
    global debuglevel

    # if level>=debuglevel:
    if True:
        print(msg, file=sys.stderr)


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
