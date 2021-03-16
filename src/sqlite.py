import sqlite3

class Create():
    def __init__(self):
        self.connect = sqlite3.connect('Sequences.db')
    def NewTable(self):
        C = self.connect.cursor()

        '''
        Insert the desired table name where it says insert table name
        '''
        C.execute('''CREATE TABLE Records (UUID TEXT PRIMARY KEY ,ProteinAccession TEXT UNIQUE, ProteinSequence TEXT, ProteinDescription TEXT,
                                  ProteinID INTEGER, CDSAccession TEXT, CDSSeq TEXT, CDSDescription TEXT, GenomicAccession TEXT,
                                GenomicSeq TEXT, GenomicDescription TEXT, GeneID INTEGER, GenomicContext TEXT, ParentDomains TEXT,
                                Introns TEXT, ExonLength TEXT, Taxonomy TEXT, CommonName TEXT)''')

        self.connect.commit()
