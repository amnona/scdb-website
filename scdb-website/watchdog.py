import os
import subprocess
import time
import datetime
import requests
from pathlib import Path


def removeFile(file_name):
    try:
        os.remove(file_name)
    except OSError:
        pass

def createSeqFile(file_name,seq_str):
    with open(file_name, "w") as text_file:
        text_file.write(">seq\n%s" % seq_str)

def saveStringToFile(file_name,content_str):
    with open(file_name, "w") as text_file:
        text_file.write(content_str)        
        
def readResultFromFile(file_name):
    ret = ""
    try:
        with open (file_name, "r") as myfile:
            ret=myfile.readlines()
    except:
        pass
    
    return ret

def isFileExist(fileName):
    my_file = Path("stop")
    if my_file.is_file():
        # file exists
        return True
    return False

def isProcessExist(processId):
    """ Check For the existence of a unix pid. """
    try:
        c=requests.post('http://127.0.0.1:5000/main')
    except:
        return False
    else:
        return True

if __name__ == '__main__':
    
    removeFile("stop")
    sleep_time = 5
    sleep_after_load = 10
    date_time_str = datetime.datetime.now().strftime("%Y-%m-%d--%H:%M:%S")
    count_execution = 0
    
    proc  = subprocess.Popen(["python","/Users/admin/scdb-website/scdb-website/Server_Main.py"])
    process_to_watch = proc.pid
    count_execution = count_execution + 1
    summary_str = "count_execution = %s\nprocess_id = %s\n" % (count_execution,process_to_watch)
    saveStringToFile("watchdog_summary_log" + date_time_str,summary_str)
    time.sleep(sleep_after_load)
    
    while isFileExist("stop") == False:
        if isProcessExist(process_to_watch) == False :
            proc  = subprocess.Popen(["python","/Users/admin/scdb-website/scdb-website/Server_Main.py"])
            process_to_watch = proc.pid
            count_execution = count_execution + 1
            summary_str = "count_execution = %s\nprocess_id = %s\n" % (count_execution,process_to_watch)
            saveStringToFile("watchdog_summary_log" + date_time_str,summary_str)
            time.sleep(sleep_after_load)
            
        time.sleep(sleep_time)
    

