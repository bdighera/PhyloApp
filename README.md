# PhyloApp
Full Stack Web Application for PhyloPy
### Requirements
python 3.7+
Linux 64 bit
### Download PhyloApp Repository
```
cd desired/clone/path
git clone https://github.com/bdighera/PhyloApp.git
```
### Create python virtual environment
```
python3 -m venv PhyloEnv
```
### Activate Virtual Environment
```
source PhyloEnv/bin/activate
```

### Install packages
```
cd PhyloPy/
pip install -r requirements.txt
```
### Set Permissions for Linux Executibles
```
cd PhyloPy/
chmod 777 execs/clustalo-1.2.0
chmod 777 execs/FastTree
chmod 777 execs/spidey.linux.64
```

### Launch Webserver
Currently runs @ http://localhost:5000/
```
cd PhyloPy/
python app.py
```

