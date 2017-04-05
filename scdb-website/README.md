# dbBact website server
## the web frontend for the dbBact bacterial knowledgebase
[dbBact](https://

# Installation
- create the conda environment

```
conda create --name dbbact-web python=3 numpy flask requests matplotlib
```

- activate it

```
source activate dbbact-web
```

- install the wordcloud module from github (the pip install version is old and not good enough)

```
pip install git+git://github.com/amueller/word_cloud
```

# running
```
export FLASK_APP=Server_Main.py
```

and then

```
flask run
```

The website should be available in:
http://127.0.0.1:5000/
