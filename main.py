from flask import Flask, render_template, request, jsonify, json
from flask_sqlalchemy import SQLAlchemy
from flask_restful import Api
from jinja2 import environment
import os, sqlite3, logging


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

	if request.method == 'GET':

		#File upload for collection of sequences
		args = request.args['name']
		print('running the following sequences: %s' % str(args))
		P = parser.argparseJSON(args)
		seqList = P.parseInput()
		for seq in seqList:
			status = collector.collectSeqs([seq])
			print(status)

		data = parser.get_all_users()
		return render_template('records.html', data=data )
	else:
		return '<h1>ERROR</h1>'

@app.route('/', methods=['GET', 'POST'])
def index():
	if request.method == 'GET':
		#GET method will happen with the refresh button and will also make new dB + table upon initializing program
		try:
			data = parser.get_all_users()
			return render_template('records.html', data=data)

		except sqlite3.OperationalError:
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
			#intronData = processor.buildIntronFig(json)
			data = parser.get_all_users()
			return render_template('introns.html', intronData=json, data=data)
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
			data = parser.get_all_users()
			return render_template('genomicContext.html', gcData=genomicContext, data=data)
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

@app.errorhandler(500)
def server_error(e):
	logging.exception('An error occurred during a request.')
	return """
	An internal error occurred: <pre>{}</pre>
	See logs for full stacktrace.
	""".format(e), 500


if __name__ == '__main__':
	app.run(host='127.0.0.1', port=8080, debug= True)
