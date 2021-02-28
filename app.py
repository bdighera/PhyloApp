import sqlite3
from flask import Flask, render_template, request, jsonify, redirect, url_for
from flask_sqlalchemy import SQLAlchemy
from flask_restful import Api, Resource
from src import parser, processor
import psycopg2


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

			if request.args.get("intronsCheckBox") == 'checkbox':
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
				return jsonify(introns)

			elif request.args.get("genomicContextCheckBox") == 'checkbox':
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

			elif request.args.get("domainsCheckBox") == 'checkbox':
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
			print('This code still needs to be added')
			pass

		else:
			return{'THIS CODE:IS NOT WORKING'}

@app.route('/', methods=['POST', 'GET'])
def index():
	return render_template('index.html')

@app.route('/dB', methods=['GET'])
def dB():
	if request.method == 'GET':
		data = parser.get_all_users()
		return render_template('records.html', data=data)




if __name__ == '__main__':
	app.run(debug= True)
