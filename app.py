from flask import Flask, render_template, request, jsonify, json
from flask_sqlalchemy import SQLAlchemy
from flask_restful import Api
from jinja2 import environment
import os, sqlite3


from src import parser, processor, collector, sqlite


app = Flask(__name__, static_url_path='/static')
api = Api(app)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///Sequences.db'
db = SQLAlchemy(app)

class SeqModel(db.Model):

	__tablename__ = 'Records'
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

@app.route('/InitialFigure', methods=['GET', 'POST'])
def InitialFigure():
	# ToDO: This is going to be deprecated when the visualization can only be run from dB page
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
				json = Phylo.buildIntrons()

				mpld3_html = processor.buildIntronFig(json)
				return render_template('intron.html', plot=mpld3_html)
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
				return render_template('index.html', data=domains)

		elif request.args.get("CollectSeqs") == 'radio':
			args = request.args['name']
			P = parser.argparseJSON(args)
			collector.collectSeqs(P.parseInput())

		else:
			return '<h1><center>404 ERROR - BROKEN PATH</center></h1>'
	#TODO: Make this the only method once the GET deprecated code above is removed
	elif request.method == 'POST':
		#File upload for collection of sequences
		seqs = ''
		if request.files['myfile']:
			file = request.files['myfile']
			seqs = str(file.read())
			args = seqs
			P = parser.argparseJSON(args)
			collector.collectSeqs(P.parseInput())


@app.route('/', methods=['GET', 'POST'])
def index():
	if request.method == 'GET':
		#GET method will happen with the refresh button and will also make new dB + table upon initializing program
		try:
			data = parser.get_all_users()
			return render_template('records.html', data=data)

		except sqlite3.OperationalError:
			os.remove('Sequences.db')
			create = sqlite.Create()
			create.NewTable()
			data = parser.get_all_users()
			return render_template('records.html', data=data)
	elif request.method == 'POST':
		#POST Job will submit what the user inputs as sequences to submit from dB page

		runtype = request.form.get('typeofrun')
		seqs = ''
		if request.form.getlist('entries'):
			seqs = request.form.getlist('entries')
		elif request.files['myfile']:
			file = request.files['myfile']
			seqs = str(file.read())
		elif request.form.getlist('entries') and request.files['myfile']:
			file = request.files['myfile']
			db = request.form.getlist('entries')
			seqs = db + file

		if runtype == 'introns':
			args = seqs

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
			json = Phylo.buildIntrons()
			mpld3_html = processor.buildIntronFig(json)
			data = parser.get_all_users()
			# ToDO: returns the index.html file in the iframe - the image displays in the iframe introns folder
			return render_template('records.html', plot=mpld3_html, data=data)
		elif runtype == 'genomicContext':
			args = seqs

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
		elif runtype == 'domains':
			args = seqs

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
			data = parser.get_all_users()
			return render_template('domain.html', domainData=domains, data=data)
	else:
		return render_template("index.html", Title='HomePage - PhyloApp')

if __name__ == '__main__':
	app.run(debug= True)

