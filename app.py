from flask import Flask, render_template, request, jsonify
from flask_sqlalchemy import SQLAlchemy
from flask_restful import Api
import matplotlib.pyplot as plt, mpld3
import pandas as pd
import json, os


from src import parser, processor, collector


app = Flask(__name__, static_url_path='/static')
api = Api(app)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///records.db'
db = SQLAlchemy(app)

class SeqModel(db.Model):

	__tablename__ = 'PTBP'
	id = db.Column(db.Text, primary_key=True)
	accession = db.Column(db.Text, unique=True)
	seq = db.Column(db.Text)
	description = db.Column(db.Text)
	proteinID = db.Column(db.Integer)
	accessionCDS = db.Column(db.Text)
	seqCDS = db.Column(db.Text)
	descriptionCDS = db.Column(db.Text)
	accessionGenomic = db.Column(db.Text)
	seqGenomic = db.Column(db.Text)
	descriptionGenomic = db.Column(db.Text)
	geneID = db.Column(db.Integer)
	genomicContext = db.Column(db.Text)
	parentDomains = db.Column(db.Text)
	introns = db.Column(db.Text)
	exonLength = db.Column(db.Text)
	taxonomy = db.Column(db.Text)
	commonName = db.Column(db.Text)

@app.route('/InitialFigure', methods=['GET'])
def InitialFigure():
	if request.method == 'GET':

		if request.args.get("DisplaySeqs") == 'radio':

			if request.args['typeofrun'] == 'introns':
				args = request.args['name']

				P = parser.argparseJSON(args)

				P.parseInput()
				P.pullDBrecords()
				data = P.serialize()

				Phylo = processor.PhyloTreeConstruction(

					proteinAccession=data['proteinAccession'],
					proteinSeq=data['proteinSeq'],
					proteinDescription=data['proteinDescription'],
					GenomicContext=data['genomicContext'],
					ParentDomains=data['parentDomains'],
					Introns=data['introns'],
					ExonLenghts=data['exonLength'],
					commonNames=data['commonNames'],
					GeneID=data['geneID'],
					image_scaling=1,
					stretch=0
				)
				introns = Phylo.buildIntrons()

				# Figure for Displaying Introns
				df = pd.DataFrame.from_dict(introns)

				Seqs = df['Sequences']
				labels = [Seqs[i]['name'] for i in range(len(Seqs))]
				yticks = []
				y = 10

				fig = plt.figure()
				for i in range(len(Seqs)):
					xList = []
					for j in range(len(Seqs[i]['introns'])):
						x = Seqs[i]['introns'][j]['realLocation']
						color = Seqs[i]['introns'][j]['background']
						# print(x,y)
						plt.scatter(x, y, color=color, zorder=2, linestyle='-')
						xList.append(x)

					x1 = max(xList)
					x2 = 0
					plt.plot([x1, x2], [y, y], linestyle='solid', zorder=1)
					yticks.append(y)
					y += 10

				plt.yticks(yticks, labels=labels)
				plt.tight_layout()
				mpld3.fig_to_html(fig=fig)
				mpld3_html = mpld3.fig_to_html(fig=fig)

				return render_template('index.html', plot=mpld3_html)

			elif  request.args['typeofrun'] == 'genomicContext':
				args = request.args['name']

				P = parser.argparseJSON(args)

				P.parseInput()
				P.pullDBrecords()
				data = P.serialize()

				Phylo = processor.PhyloTreeConstruction(

					proteinAccession=data['proteinAccession'],
					proteinSeq=data['proteinSeq'],
					proteinDescription=data['proteinDescription'],
					GenomicContext=data['genomicContext'],
					ParentDomains=data['parentDomains'],
					Introns=data['introns'],
					ExonLenghts=data['exonLength'],
					commonNames=data['commonNames'],
					GeneID=data['geneID'],
					image_scaling=1,
					stretch=0
				)
				genomicContext = Phylo.buildGenomicContext()
				return jsonify(genomicContext)

			elif request.args['typeofrun'] == 'domains':
				args = request.args['name']

				P = parser.argparseJSON(args)

				P.parseInput()
				P.pullDBrecords()
				data = P.serialize()

				Phylo = processor.PhyloTreeConstruction(

					proteinAccession=data['proteinAccession'],
					proteinSeq=data['proteinSeq'],
					proteinDescription=data['proteinDescription'],
					GenomicContext=data['genomicContext'],
					ParentDomains=data['parentDomains'],
					Introns=data['introns'],
					ExonLenghts=data['exonLength'],
					commonNames=data['commonNames'],
					GeneID=data['geneID'],
					image_scaling=1,
					stretch=0
				)
				domains = Phylo.buildDomains()
				return jsonify(domains)

		elif request.args.get("CollectSeqs") == 'radio':
			args = request.args['name']

			P = parser.argparseJSON(args)
			input = P.parseInput()

			collector.collectSeqs(input)

		else:
			return{'THIS CODE':'IS NOT WORKING'}



@app.route('/')
def index():
	filename = os.path.join('../','static','images', 'orthologo.png')
	print(filename)
	return render_template("index.html", user_image=filename)


@app.route('/dB', methods=['GET'])
def dB():
	if request.method == 'GET':
		data = parser.get_all_users()
		return render_template('records.html', data=data)



if __name__ == '__main__':
	app.run(debug= True)
from flask import Flask, render_template, request, jsonify
from flask_sqlalchemy import SQLAlchemy
from flask_restful import Api
import matplotlib.pyplot as plt, mpld3
import pandas as pd
import json, os


from src import parser, processor, collector


app = Flask(__name__, static_url_path='/static')
api = Api(app)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///records.db'
db = SQLAlchemy(app)

class SeqModel(db.Model):

	__tablename__ = 'PTBP'
	id = db.Column(db.Text, primary_key=True)
	accession = db.Column(db.Text, unique=True)
	seq = db.Column(db.Text)
	description = db.Column(db.Text)
	proteinID = db.Column(db.Integer)
	accessionCDS = db.Column(db.Text)
	seqCDS = db.Column(db.Text)
	descriptionCDS = db.Column(db.Text)
	accessionGenomic = db.Column(db.Text)
	seqGenomic = db.Column(db.Text)
	descriptionGenomic = db.Column(db.Text)
	geneID = db.Column(db.Integer)
	genomicContext = db.Column(db.Text)
	parentDomains = db.Column(db.Text)
	introns = db.Column(db.Text)
	exonLength = db.Column(db.Text)
	taxonomy = db.Column(db.Text)
	commonName = db.Column(db.Text)

@app.route('/InitialFigure', methods=['GET'])
def InitialFigure():
	if request.method == 'GET':

		if request.args.get("DisplaySeqs") == 'radio':

			if request.args['typeofrun'] == 'introns':
				args = request.args['name']

				P = parser.argparseJSON(args)

				P.parseInput()
				P.pullDBrecords()
				data = P.serialize()

				Phylo = processor.PhyloTreeConstruction(

					proteinAccession=data['proteinAccession'],
					proteinSeq=data['proteinSeq'],
					proteinDescription=data['proteinDescription'],
					GenomicContext=data['genomicContext'],
					ParentDomains=data['parentDomains'],
					Introns=data['introns'],
					ExonLenghts=data['exonLength'],
					commonNames=data['commonNames'],
					GeneID=data['geneID'],
					image_scaling=1,
					stretch=0
				)
				introns = Phylo.buildIntrons()

				# Figure for Displaying Introns
				df = pd.DataFrame.from_dict(introns)

				Seqs = df['Sequences']
				labels = [Seqs[i]['name'] for i in range(len(Seqs))]
				yticks = []
				y = 10

				fig = plt.figure()
				for i in range(len(Seqs)):
					xList = []
					for j in range(len(Seqs[i]['introns'])):
						x = Seqs[i]['introns'][j]['realLocation']
						color = Seqs[i]['introns'][j]['background']
						# print(x,y)
						plt.scatter(x, y, color=color, zorder=2, linestyle='-')
						xList.append(x)

					x1 = max(xList)
					x2 = 0
					plt.plot([x1, x2], [y, y], linestyle='solid', zorder=1)
					yticks.append(y)
					y += 10

				plt.yticks(yticks, labels=labels)
				plt.tight_layout()
				mpld3.fig_to_html(fig=fig)
				mpld3_html = mpld3.fig_to_html(fig=fig)

				return render_template('index.html', plot=mpld3_html)

			elif  request.args['typeofrun'] == 'genomicContext':
				args = request.args['name']

				P = parser.argparseJSON(args)

				P.parseInput()
				P.pullDBrecords()
				data = P.serialize()

				Phylo = processor.PhyloTreeConstruction(

					proteinAccession=data['proteinAccession'],
					proteinSeq=data['proteinSeq'],
					proteinDescription=data['proteinDescription'],
					GenomicContext=data['genomicContext'],
					ParentDomains=data['parentDomains'],
					Introns=data['introns'],
					ExonLenghts=data['exonLength'],
					commonNames=data['commonNames'],
					GeneID=data['geneID'],
					image_scaling=1,
					stretch=0
				)
				genomicContext = Phylo.buildGenomicContext()
				return jsonify(genomicContext)

			elif request.args['typeofrun'] == 'domains':
				args = request.args['name']

				P = parser.argparseJSON(args)

				P.parseInput()
				P.pullDBrecords()
				data = P.serialize()

				Phylo = processor.PhyloTreeConstruction(

					proteinAccession=data['proteinAccession'],
					proteinSeq=data['proteinSeq'],
					proteinDescription=data['proteinDescription'],
					GenomicContext=data['genomicContext'],
					ParentDomains=data['parentDomains'],
					Introns=data['introns'],
					ExonLenghts=data['exonLength'],
					commonNames=data['commonNames'],
					GeneID=data['geneID'],
					image_scaling=1,
					stretch=0
				)
				domains = Phylo.buildDomains()
				return jsonify(domains)

		elif request.args.get("CollectSeqs") == 'radio':
			args = request.args['name']

			P = parser.argparseJSON(args)
			input = P.parseInput()

			collector.collectSeqs(input)

		else:
			return{'THIS CODE':'IS NOT WORKING'}



@app.route('/')
def index():
	filename = os.path.join('../','static','images', 'orthologo.png')
	print(filename)
	return render_template("index.html", user_image=filename)


@app.route('/dB', methods=['GET'])
def dB():
	if request.method == 'GET':
		data = parser.get_all_users()
		return render_template('records.html', data=data)



if __name__ == '__main__':
	app.run(debug= True)
