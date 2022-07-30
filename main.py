from flask import Flask, render_template, request, jsonify, json, redirect, url_for, session
from flask_sqlalchemy import SQLAlchemy
from flask_restful import Api
from flask_cors import CORS
import secrets
from jinja2 import environment
import os, sqlite3, logging
from time import sleep
from flask_session import Session


from src import parser, processor, collector, sqlite


app = Flask(__name__)
app.secret_key = secrets.token_urlsafe(16)
app.config['CORS_HEADERS'] = 'Content-Type'

CORS(app)
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

@app.route('/deleteRow', methods=['GET', 'POST'])
#TODO: Delete a row from the records table.
def DeleteRow():
	print("hi")
	if request.method == 'DELETE' :
		for id in request.args['delete']:
			parser.deleteRow(id)
	

@app.route('/InitialFigure', methods=['GET', 'POST'])
#TODO: Update the name of this to reflect collection of sequences.
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
#TODO: Update the name of this to reflect the image displays
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
		print(request.form.get('action'))
		if(request.form.get('action') == 'Delete'):
			#msa = request.form.get('compared_motifs')
			#if msa:
				#msa = msa.replace('"', '')
				#print(msa)
			rowIds = request.form.getlist('rowId')
			entrieIds = request.form.getlist('entries')
			
			for entrie in entrieIds:
				rowId = request.form.get(entrie)
				parser.deleteRow(rowId)
			
			data = parser.get_all_users()
			return render_template('records.html', data=data )
		msa = request.form.get('compared_motifs')
		if msa:
			msa = msa.replace('"', '')
			seq1 = msa.split(',')[0]
			seq2 = msa.split(',')[1]
			alignment = processor.MSA(seq1, seq2)
			print(alignment)
			return alignment


		#POST Job will submit what the user inputs as sequences to submit from dB page
		runtype = request.form.get('typeofrun')
		seqs=''
		if runtype:
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
			#Deprecated code - makes intron figure using mpld3
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
		elif runtype == 'tree':
			#os.remove(path='static/images/tree_img.svg')
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
			#TODO: Need to figure out why the image is not regenerating, maybe need to figure out how to drop it in dynamically?
			Phylo.buildTree()
			data = parser.get_all_users()
			return render_template('tree.html', data=data)
		elif runtype == 'MSA':
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

			data = parser.get_all_users()

			msa = Phylo.collectMultipleSequencingAlignment()
			msa = {str(i.id).split('_')[0]+'_'+str(i.id).split('_')[1]+' '+str(i.id).split('_')[-1]:str(i.seq) for i in msa}

			return render_template('MSA.html', data=data, msa=msa)

		figs = request.form.getlist('msa_entries')
		if figs:
			figs_seqs = [i.split(' ')[0] for i in figs]
			P = parser.argparseJSON(figs_seqs)

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

			gc = Phylo.buildGenomicContext()
			domains = Phylo.buildDomains()
			introns = Phylo.buildIntrons()
			print(gc)
			print(domains)
			print(introns)

			return render_template('genomicContext.html', data=data, genomicContext=gc, domains=domains, introns=introns)


		else:
			return '<h1>ERROR</h1>'


@app.route('/msa', methods=['POST', 'GET'])
def msa():
	pass

@app.errorhandler(500)
def server_error(e):
	logging.exception('An error occurred during a request.')
	return """
	An internal error occurred: <pre>{}</pre>
	See logs for full stacktrace.
	""".format(e), 500

if __name__ == '__main__':
	app.run(host='127.0.0.1', port=8080, debug=True, threaded=True)


