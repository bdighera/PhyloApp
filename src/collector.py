from Bio import Entrez
import subprocess, sqlite3
from randomcolor import RandomColor
from time import sleep
from Bio import SeqIO
from src import handler
import uuid

def collectSeqs(formInput):

    timer = 6
    # sets the timer value between ping requests to the NCBI server

    # Iterates through each parsed record in the input file
    for accession in formInput:

        try:

            # Uses the record accession number to check if it is already in the SQL database, if not: continue.
            SQL = SQLiteChecker(accession)

            # If accession number is not located in the DB then it will run through the pipeline
            if SQL.checkRecords() != True:

                # Initialize Sequence Collector class with input of the working protein accession record
                Seq = SequenceCollector(accession)

                # collects the protein ID for the current working protein
                # protein ID is needed in order to establish relationship between NCBI databases ie. Gene, Protein, etc.
                # runs timer after the method is executed in order to not ping NCBI server too fast - or they will SHUT YOU DOWN
                id = Seq.collectProteinIDs()
                sleep(timer)

                # Runs the method to collect the CDS and taxonomy
                # TO understand this method better (along with other methods of class Sequence Collector) visit SequenceCOllector.py in CWD
                CDS, Taxonomy = Seq.collectCDS()
                sleep(timer)

                # Collects the entire sequence of the gene and the gene ID
                Genomic, GeneID = Seq.collectGenomicSeq()
                sleep(timer)

                # Collects the intron phases and exon lengths - this method requires executible files from ./exec folder
                # Internal function -does not ping NCBI server
                IntronPhases, ExonLengths = Seq.collectIntrons(CDS, Genomic)

                # Collects the domains for the protein
                Domains = Seq.collectProteinDomains()
                sleep(timer)

                # Initialize Genomic Context class
                # Collects neighboring genes to parent protein, coding direction, and domains of those genes
                GC = GenomicContext(Genomic)

                # Collects record for  +/- 50k basepairs up/downstream of parent gene
                gcRecord = GC.fetchRecord()
                sleep(timer)

                # parses raw genomic context data collected from NCBI
                parsedGCrecord = GC.parseRecord(gcRecord)

                commonName = Seq.collectCommonName(id)

                # Initialize the SQLite class with all of the data for the working protein
                SQL = SQliteRecordInput(Seq.proteinSeq,
                                        Seq.proteinAccession,
                                        Seq.proteinDescription,
                                        Seq.proteinID,
                                        CDS,
                                        Genomic,
                                        GeneID,
                                        parsedGCrecord,
                                        Domains,
                                        IntronPhases,
                                        ExonLengths,
                                        Taxonomy,
                                        commonName)

                # Uploads records to the DB
                # REMEMBER TO CHANGE THE NAME OF THE TABLE TO REFLECT WHERE YOU WANT THE RECORDS TOGO!
                SQL.uploadRecords()
                print('Sequence Uploaded!')

            # If protein accession number is in the database then it will skip over and go to the next record
            else:
                print('Record ' + accession + ' is already in database. Proceeding to next record')
                pass

        # Errors in some of the sequences do occur on the NCBI server end. For this errors the DNAJC_Errors file is populated with those accession #s
        except IndexError:
            print('Index Error Occurred on: ' + str(accession))

class SequenceCollector():

    def __init__(self, proteinRecord):
        self.proteinRecord = proteinRecord
        self.proteinID = ''
        self.proteinDescription = ''
        self.proteinName = ''
        self.proteinSeq = ''
        self.proteinAccession = ''

    def collectProteinIDs(self):
        Entrez.email = 'bdighera@csu.fullerton.edu'

        proteinAccession = self.proteinRecord
        print(proteinAccession)
        eSearch = Entrez.esearch(db='protein',term= proteinAccession, rettype='gb', api_key='8f07f1fa1ec8252ef42382cb722e1749b608' )
        result = Entrez.read(eSearch, validate=False)

        proteinID = result['IdList'][0]

        self.proteinID = proteinID

        return proteinID

    def collectCDS(self):
        Entrez.email = 'bdighera@csu.fullerton.edu'

        proteinID = self.proteinID

        generalEfetch= SeqIO.read(Entrez.efetch(db='protein', id=proteinID, rettype='gb', retmode='text', api_key='8f07f1fa1ec8252ef42382cb722e1749b608'), 'genbank')

        self.proteinDescription = generalEfetch.description
        self.proteinName = generalEfetch.name
        self.proteinSeq = generalEfetch.seq
        self.proteinAccession = generalEfetch.id


        Taxonomy = generalEfetch.annotations['taxonomy']

        result = [result.qualifiers['coded_by'] for result in generalEfetch.features if result.type == 'CDS'][0][0]

        accession, codingRegion = result.split(':')
        startseq, endseq = codingRegion.split('..')

        recordEfetch = SeqIO.read(Entrez.efetch(db='nuccore', id=accession, seq_start=startseq, seq_top=endseq, rettype='fasta', api_key='8f07f1fa1ec8252ef42382cb722e1749b608'), 'fasta')

        return recordEfetch, Taxonomy

    def collectGenomicSeq(self):
        Entrez.email = 'bdighera@csu.fullerton.edu'

        proteinID = self.proteinID

        elinkResult = Entrez.read(Entrez.elink(db='gene', dbfrom='protein', id=proteinID, api_key='8f07f1fa1ec8252ef42382cb722e1749b608'))
        geneID = elinkResult[0]['LinkSetDb'][0]['Link'][0]['Id']

        generalEfetch= Entrez.read(Entrez.efetch(db='gene', id=geneID, rettype='fasta', retmode='xml', api_key='8f07f1fa1ec8252ef42382cb722e1749b608'), validate=False)

        accession = generalEfetch[0]['Entrezgene_locus'][0]['Gene-commentary_accession'] + '.' + generalEfetch[0]['Entrezgene_locus'][0]['Gene-commentary_version']
        startseq = generalEfetch[0]['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from']
        endseq = generalEfetch[0]['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to']

        genomeEfetch = SeqIO.read(Entrez.efetch(db='nuccore', id=accession, seq_start=int(startseq), seq_stop=int(endseq), rettype='fasta', api_key='8f07f1fa1ec8252ef42382cb722e1749b608'),'fasta')

        return genomeEfetch, geneID

    def collectIntrons(self, CDS, Genomic):

        CDSseq = str('>' + CDS.description + '\n' + CDS.seq)
        GenomicSeq = str('>' + Genomic.description + '\n' + Genomic.seq)

        I = handler.IntronFileHandler()

        CDSfilePath = I.getCDSTmpPath()
        GenomicfilePath = I.getGenomicTmpPath()

        I.clearPreviousInput(CDSfilePath)
        I.clearPreviousInput(GenomicfilePath)

        I.writeNewInput(CDSfilePath, CDSseq)
        I.writeNewInput(GenomicfilePath, GenomicSeq)

        IntronCalculatorExec = I.getSpideyPath()

        #ToDO: This is removed for beta testing
        # intronProcessoutput = subprocess.Popen([IntronCalculatorExec,
        #                                    "-i", GenomicfilePath,
        #                                    "-m", CDSfilePath,
        #                                    "-p", "1"], stdout=subprocess.PIPE)
        # 
        # return I.spideyOutputParser(intronProcessoutput.stdout.readlines())
        
        return 'NONE'

    def collectProteinDomains(self):
        Entrez.email ='bdighera@csu.fullerton.edu'

        proteinAccession = self.proteinAccession

        e_fetch = SeqIO.parse(Entrez.efetch(db='protein', id=proteinAccession, retmax=1000, rettype='gb', retmode='fasta',api_key='8f07f1fa1ec8252ef42382cb722e1749b608'), 'gb')

        resultFeatures = [result.features for result in e_fetch][0]

        rawDomains = [result for result in resultFeatures if result.type == 'Region']

        domainLocation =[str(result.location).split('[')[1].split(']')[0] for result in rawDomains]
        domainName = [result.qualifiers['region_name'][0] for result in rawDomains]


        Domains = {domainName[i]:domainLocation[i] for i in range(len(domainLocation))}

        return Domains

    def collectCommonName(self, id):
        '''
        :param fileHandle: path to file containing protein accessions, will collect record from database
        :return: Does not return anything, populates database file with updated info
        '''

        try:

            elinkResult = Entrez.read(
                Entrez.elink(db='taxonomy', dbfrom='protein', id=id,
                             api_key='8f07f1fa1ec8252ef42382cb722e1749b608'))

            tax_id = elinkResult[0]['LinkSetDb'][0]['Link'][0]['Id']

            taxonomy = Entrez.read(
                Entrez.efetch(db='taxonomy', id=tax_id, api_key='8f07f1fa1ec8252ef42382cb722e1749b608'))

            GenBankCommonName = taxonomy[0]['OtherNames']['GenbankCommonName']
            return GenBankCommonName

        except KeyError:
            return 'Common Name Unknown'
        except IndexError:
            try:
                GenBankCommonName = taxonomy[0]['OtherNames']['GenbankCommonName']
                return GenBankCommonName
            except KeyError:
                return 'Common Name Unknown'

class GenomicContext():
    def __init__(self, geneRecord, kbps=50000):

        self.geneRecord = geneRecord
        self.timer = 0
        self.run_name = ''
        self.input_protein_accession_number = ['NP_005336.3']
        self.input_protein_sequence = []
        self.retrieved_protein_ids = []
        self.retrieved_full_cds = []
        self.parsed_cds_acession = []
        self.gene_accession_list = []
        self.gene_start_sequence_list = []
        self.gene_end_sequence_list = []
        self.gene = []
        self.intron_phase = []
        self.exon_lengths = []
        self.accession_dict_with_introns = {}
        self.number_protein_IDs = 0
        self.number_gene_IDs = 0
        self.number_CDS_IDs = 0
        self.msa_aligned_protein = []
        self.genomic_context = {}
        self.genomic_context_coords = ''
        self.genomic_context_motifs = {}
        self.trueGenestart = 0
        self.trueGeneend = 0
        self.trueGenedirection = ''
        self.genomic_context_gene_counter = 0
        self.protein_name = 'STRING'
        self.kbps = kbps

    # Collects the genomic context of the parent protein
    def fetchRecord(self):
        Entrez.email = 'bdighera@csu.fullerton.edu'

        geneRecord = self.geneRecord

        #Gets gene accession number
        accession = str(geneRecord.description.split(':')[0])
        #Gets the start sequence of the gene and moves 50k basepairs upstream
        startseq = int(geneRecord.description.split(' ')[0].split(':')[1].split('-')[0])

        #Gets the end sequence of the gene and moves 50k basepairs downstream
        endseq = int(geneRecord.description.split(' ')[0].split(':')[1].split('-')[0])

        # Begins the efetch call to acquire all the neighboring genes of the gene of interest (genomic context efetch call)

        GCfetch= Entrez.read(Entrez.efetch(db="nuccore",
                                                           id=accession,
                                                           seq_start=startseq,
                                                           seq_stop=endseq,
                                                           rettype='gb',
                                                           retmode='xml', validate=False,
                                                           api_key='8f07f1fa1ec8252ef42382cb722e1749b608'))

        return GCfetch

    # Goes through the Gene portion of the genomic context fetch, (Entrez Dictionary) and grabs the gene direction, name, and start/end point
    def parseRecord(self, GCrecord):

        #Working on rewriting code to be more pythonic, start is the next 3 lines:
        #rawCDSRecord = [record['GBFeature_quals'] for record in GCrecord[0]['GBSeq_feature-table'] if record['GBFeature_key'] == 'CDS'][0]
        #rawGeneRecord = [record for record in GCrecord[0]['GBSeq_feature-table'] if record['GBFeature_key'] == 'gene']
        #geneName = [record['GBQualifier_value'] for record in rawCDSRecord if record['GBQualifier_name'] == 'gene']

        gb_fetch = GCrecord[0]
        gene_name_list = []
        gene_start_list = []
        gene_end_list = []
        gene_id_list = []
        gene_direction_list = []
        protein_accession_dict = {}
        i = -1

        for gb_gene in gb_fetch['GBSeq_feature-table']:

            if gb_gene['GBFeature_key'] == 'CDS':

                for protein_info in gb_gene['GBFeature_quals']:

                    if protein_info['GBQualifier_name'] == 'gene':
                        name = protein_info['GBQualifier_value']

                    if protein_info['GBQualifier_name'] == 'protein_id':
                        protein_accession = protein_info['GBQualifier_value']

                        protein_accession_dict[name] = protein_accession

        for gb_gene in gb_fetch['GBSeq_feature-table']:

            if gb_gene['GBFeature_key'] == 'gene':

                for gb_info in gb_gene['GBFeature_quals']:

                    if gb_info['GBQualifier_name'] == 'gene':
                        name = str(gb_info['GBQualifier_value'])

                        try:
                            if protein_accession_dict[name]:

                                gene_name_list.append(name)
                                i += 1
                                start = int(gb_gene['GBFeature_intervals'][0]['GBInterval_from'])
                                gene_start_list.append(start)
                                end = int(gb_gene['GBFeature_intervals'][0]['GBInterval_to'])
                                gene_end_list.append(end)

                                if start < end:
                                    direction = "+"
                                    gene_direction_list.append(direction)
                                else:
                                    direction = "-"
                                    gene_direction_list.append(direction)

                                try:
                                    if 'GeneID:' in gb_info['GBQualifier_value']:
                                        gene_id = int(gb_info['GBQualifier_value'].replace('GeneID:', ''))
                                        gene_id_list.append(gene_id)
                                except:
                                    print("No GBQualifier_value")
                        except:
                            pass


                    if gb_info['GBQualifier_name'] == 'db_xref':
                        try:
                            if protein_accession_dict[name]:
                                if 'GeneID:' in gb_info['GBQualifier_value']:
                                    gene_id = str(gb_info['GBQualifier_value'])
                                    gene_id_list.append(gene_id.replace('GeneID:', ''))

                        except:
                            pass
        GC_List = []
        GC = GenomicContext(self.geneRecord)
        for i in range(len(gene_name_list)):

            completeDomains, domain = GC.fetchDomains(protein_accession_dict[gene_name_list[i]])

            GC_List.append({
                'gene_name':gene_name_list[i],
                'gene_start_seq':gene_start_list[i],
                'gene_end_seq':gene_end_list[i],
                'coding_direction':gene_direction_list[i],
                'protein_accession':protein_accession_dict[gene_name_list[i]],
                'domain': list(completeDomains[0].values())[0],
                'gene_id':gene_id_list[i]

            })

        return GC_List

    # Domain Fetch for Genomic Context
    def fetchDomains(self, gc_accession):

        Complete_Domains = []
        domainNameList = []

        e_fetch = SeqIO.parse(
            Entrez.efetch(db='protein', id="%s" % gc_accession, retmax=1000, rettype='gb', retmode='fasta',
                          api_key='8f07f1fa1ec8252ef42382cb722e1749b608'), 'gb')

        for seq_record in e_fetch:
            domain_list = []
            accession_number = seq_record.id

            for i in range(len(seq_record.features)):
                if seq_record.features[i].type == 'Region':
                    domain_location = str(seq_record.features[i].location).split('[')[1].split(']')[0]
                    domain_name = str(seq_record.features[i].qualifiers['region_name'][0])
                    domainNameList.append(domain_name)

                    domain_list.append([domain_name, domain_location])

            Complete_Domains.append(dict([(accession_number, domain_list)]))

        rand_color = RandomColor()
        Domains = [domain for domain in set(domainNameList)]

        return Complete_Domains, Domains

    def domain_colors(self, completed_domains):

        parent_protein_domains = completed_domains

        # pprint.pprint(parent_protein_domains)

        domain_list = []

        # Goes in and pulls out all of the individual domains for all of the proteins in the genomic context. List is cleaned removing duplicate domains from GC proteins
        for child_protein_domains in parent_protein_domains:

            for each_domain_list in child_protein_domains:

                if each_domain_list['domain'] != []:

                    for i in range(len(each_domain_list['domain'])):
                        domain_list.append(each_domain_list['domain'][i])

        domain_list = list(dict.fromkeys(domain_list))

        rand_color = RandomColor()

        domains_dict_colors = {domain: rand_color.generate()[0] for domain in set(domain_list)}

        # Assigns the correct color to the domains. This makes it so all gc domains are assigned the same color based on their name
        for child_protein_domains in parent_protein_domains:

            for each_domain_list in child_protein_domains:

                if each_domain_list['domain'] != []:

                    color_list = []
                    for i in range(len(each_domain_list['domain'])):
                        domain = each_domain_list['domain'][i]
                        correct_color = domains_dict_colors.get(domain)
                        domain_color_pair = {domain: correct_color}

                        color_list.append(domain_color_pair)

                        each_domain_list['color'] = color_list

                else:

                    each_domain_list['color'] = []

        return parent_protein_domains

class SQLiteChecker():
    def __init__(self, proteinAccession, dbfile='Sequences.db'):
        self.dbfile= dbfile
        self.connect = sqlite3.connect(dbfile)
        self.proteinAccessionList = proteinAccession

    # Checks to make sure that the record isn't already in the table before it runs it. Checks against the accession number of the protein fasta
    def checkRecords(self):
        #Will check all of the records in all of the tables

        dataList = []

        C = self.connect.cursor()

        for accession in self.proteinAccessionList:

            C.execute('SELECT * FROM Records WHERE ProteinAccession= (?)''', (accession,))

            data = C.fetchall()

            if data != []:
                dataList.append(data)

            #ToDO: Make a switch for all tables, if uncommented then it will search all tables. Cant have repeats in diff tables tho
            # #SQL statement getting all of the table names ie. DNAJC, HSP70
            # C.execute("SELECT name FROM sqlite_master WHERE type='table'")
            #
            # for tablename in C.fetchall():
            #     tablename = tablename[0]
            #
            #     #Checking all of the tables for the accession number
            #     C.execute('SELECT * FROM {t} WHERE ProteinAccession= (?)'''.format(t=tablename), (accession,))
            #
            #     #Checking only HSP70s table
            #     #C.execute('SELECT * FROM HSP70s WHERE ProteinAccession= (?)'''.format(t=tablename), (accession,))
            #
            #     data = C.fetchall()
            #
            #     if data != []:
            #         dataList.append(data)

        return dataList

class SQliteRecordInput():

    def __init__(self, Seq, Accession, Description, ProteinID, CDS, Genomic, GeneID, GC, Domains, IntronPhase, ExonLength, Taxonomy, CommonName):
        self.conn = sqlite3.connect('Sequences.db')
        self.ProteinSeq = str(Seq)
        self.ProteinAccession = Accession
        self.ProteinDescription =Description
        self.ProteinID = int(ProteinID)
        self.CDS = CDS
        self.CDSSeq = str(CDS.seq)
        self.CDSAccession = CDS.id
        self.CDSDescription = CDS.description
        self.Genomic = Genomic
        self.GenomicSeq = str(Genomic.seq)
        self.GenomicAccession = str(Genomic.id)
        self.GenomicDescription = Genomic.description
        self.GeneID = int(GeneID)
        self.GC = str(GC)
        self.Domains = str(Domains)
        self.Introns = str(IntronPhase)
        self.ExonLength = str(ExonLength)
        self.uuid = str(uuid.uuid4())
        self.taxonomy = str(Taxonomy)
        self.CommonName = CommonName

    def uploadRecords(self):

        C = self.conn.cursor()

        # C.execute('''CREATE TABLE HSP70 (UUID TEXT PRIMARY KEY ,ProteinAccession TEXT UNIQUE, ProteinSequence TEXT, ProteinDescription TEXT,
        #                           ProteinID INTEGER, CDSAccession TEXT, CDSSeq TEXT, CDSDescription TEXT, GenomicAccession TEXT,
        #                         GenomicSeq TEXT, GenomicDescription TEXT, GeneID INTEGER, GenomicContext TEXT, Introns TEXT, ExonLength TEXT)''')

        try:
            C.execute('''INSERT INTO Records VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',(self.uuid, self.ProteinAccession, self.ProteinSeq, self.ProteinDescription, self.ProteinID,
                                                                                  self.CDSAccession, self.CDSSeq, self.CDSDescription, self.GenomicAccession,
                                                                                  self.GenomicSeq, self.GenomicDescription,
                                                                                  self.GeneID,
                                                                                  self.GC,
                                                                                    self.Domains,
                                                                                  self.Introns,
                                                                                  self.ExonLength,
                                                                                    self.taxonomy,
                                                                                          self.CommonName))
        except sqlite3.IntegrityError:
            pass
        self.conn.commit()
        self.conn.close()



