import sqlite3, ast, json, argparse
from Bio import SeqIO

class argparseJSON():

    def __init__(self, seqs):
        self.seqs = seqs.split(',')

    def parseInput(self):

        if self.seqs != None:
            return self.seqs

    def pullDBrecords(self, dbfile= 'Records.db'):
        SQL = SQLiteChecker(records= self.seqs, dbfile= dbfile)

        self.records = SQL.checkRecords()
        return self.records

    def serialize(self):
        struct = {}
        struct.update({'proteinAccession':[record[0][1] for record in self.records]})
        struct.update({'proteinSeq': [record[0][2] for record in self.records]})
        struct.update({'proteinDescription': [record[0][3] for record in self.records]})
        struct.update({'geneID': [record[0][11] for record in self.records]})
        struct.update({'genomicContext': [ast.literal_eval(record[0][12]) for record in self.records]})
        struct.update({'parentDomains': [ast.literal_eval(record[0][13]) for record in self.records]})
        struct.update({'introns': [record[0][14] for record in self.records]})
        struct.update({'exonLength': [record[0][15] for record in self.records]})
        struct.update({'taxonomy': [record[0][16] for record in self.records]})
        struct.update({'commonNames': [record[0][17] for record in self.records]})
        return struct

    def cmdline(self):

        # Initializes the parsing function
        parser = argparse.ArgumentParser()

        # Adds the parsing arguments, sets default value, type, and help response

        parser.add_argument('-infile', default=bool(False), type=str, help='FASTA formatted protein list, .fa')
        parser.add_argument('-dir', default=os.getcwd(), type=str,
                            help='Working Directory, Default current directory')
        parser.add_argument('-rec', default=bool(False), type=str,
                            help='Recursively runs for protein files in other directories - BETA')
        parser.add_argument('-outfile', default=os.getcwd(), type=str, help='Path to output file directory')
        parser.add_argument('-genomic', default=bool(False), type=str, help='FASTA formatted genomic list, .fa')
        parser.add_argument('-CDS', default=bool(False), type=str, help='FASTA formatted CDS list, .fa')
        parser.add_argument('-tree', default=bool(False), type=str, help='Newick Tree object, .nwk')

        args = vars(parser.parse_args())

        return args['infile'], args['dir'], args['rec'], args['outfile'], args['genomic'], args['CDS'], args[
            'tree']

class argparseFile():

    def __init__(self, path):
        self.path = path
        self.seqs = [seqs]
        self.input = ''

    def parseInput(self):

        if self.path != None:
            file = SeqIO.parse(self.path, 'fasta')
            self.input = [item.id for item in file]
            return self.input


    def pullDBrecords(self, dbfile= 'Records.db'):
        SQL = SQLiteChecker(records= self.input, dbfile= dbfile)

        self.records = SQL.checkRecords()
        return self.records

    def serialize(self):
        struct = {}
        struct.update({'proteinAccession':[record[0][1] for record in self.records]})
        struct.update({'proteinSeq': [record[0][2] for record in self.records]})
        struct.update({'proteinDescription': [record[0][3] for record in self.records]})
        struct.update({'geneID': [record[0][11] for record in self.records]})
        struct.update({'genomicContext': [ast.literal_eval(record[0][12]) for record in self.records]})
        struct.update({'parentDomains': [ast.literal_eval(record[0][13]) for record in self.records]})
        struct.update({'introns': [record[0][14] for record in self.records]})
        struct.update({'exonLength': [record[0][15] for record in self.records]})
        struct.update({'taxonomy': [record[0][16] for record in self.records]})
        struct.update({'commonNames': [record[0][17] for record in self.records]})
        return struct

    def cmdline(self):

        # Initializes the parsing function
        parser = argparse.ArgumentParser()

        # Adds the parsing arguments, sets default value, type, and help response

        parser.add_argument('-infile', default=bool(False), type=str, help='FASTA formatted protein list, .fa')
        parser.add_argument('-dir', default=os.getcwd(), type=str,
                            help='Working Directory, Default current directory')
        parser.add_argument('-rec', default=bool(False), type=str,
                            help='Recursively runs for protein files in other directories - BETA')
        parser.add_argument('-outfile', default=os.getcwd(), type=str, help='Path to output file directory')
        parser.add_argument('-genomic', default=bool(False), type=str, help='FASTA formatted genomic list, .fa')
        parser.add_argument('-CDS', default=bool(False), type=str, help='FASTA formatted CDS list, .fa')
        parser.add_argument('-tree', default=bool(False), type=str, help='Newick Tree object, .nwk')

        args = vars(parser.parse_args())

        return args['infile'], args['dir'], args['rec'], args['outfile'], args['genomic'], args['CDS'], args[
            'tree']

class SQLiteChecker():
    def __init__(self, records, dbfile):
        self.dbfile= dbfile
        self.connect = sqlite3.connect(dbfile)
        self.records = records

    def checkRecords(self):


        dataList = []

        C = self.connect.cursor()

        for accession in self.records:

            C.execute('SELECT * FROM PTBP WHERE ProteinAccession= (?)''', (accession,))

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

def JSONtofile(data, filename):
    #clear previous output from file
    with open(filename, 'w') as clearfile:
        clearfile.write('')

    #write the json output to file
    for i in range(len(data)):
        with open(filename, 'a') as outfile:
            json.dump(data[i], outfile)


class SQLiteGetAll():
    def __init__(self, dbfile):
        self.dbfile = dbfile

    def getAllRecords(self):
        C = self.connect.cursor()

        C.execute('SELECT * FROM PTBP')

        self.records = C.fetchall()

        return self.records

    def serialize(self):
        struct = {}
        struct.update({'proteinAccession':[record[0][1] for record in self.records]})
        struct.update({'proteinSeq': [record[0][2] for record in self.records]})
        struct.update({'proteinDescription': [record[0][3] for record in self.records]})
        struct.update({'geneID': [record[0][11] for record in self.records]})
        struct.update({'genomicContext': [ast.literal_eval(record[0][12]) for record in self.records]})
        struct.update({'parentDomains': [ast.literal_eval(record[0][13]) for record in self.records]})
        struct.update({'introns': [record[0][14] for record in self.records]})
        struct.update({'exonLength': [record[0][15] for record in self.records]})
        struct.update({'taxonomy': [record[0][16] for record in self.records]})
        struct.update({'commonNames': [record[0][17] for record in self.records]})
        return struct