from Bio.Align.Applications import ClustalOmegaCommandline
import itertools, os, subprocess, sys, math, dendropy, pprint, ast, uuid
from randomcolor import RandomColor
from src import handler
from ete3 import Tree, TreeStyle, TextFace

class PhyloTreeConstruction(object):

    def __init__(self, proteinAccession, proteinSeq, proteinDescription, GenomicContext, ParentDomains, Introns, ExonLenghts,
				 commonNames, GeneID, image_scaling, stretch):
        self.proteinAccessions = proteinAccession
        self.proteinSeqs = proteinSeq
        self.proteinDescs = proteinDescription
        self.GenomicContexts = GenomicContext
        self.parentDomains = ParentDomains
        self.Introns = Introns
        self.exonLengths = ExonLenghts
        self.commonNames = commonNames
        self.collectMultipleSequencingAlignment()
        self.rootedTreeConstruction()
        self.geneID = GeneID
        self.image_scaling = image_scaling
        self.stretch = stretch


    def rootedTreeConstruction(self):

        in_file = os.path.join('execs', 'tmp', 'aligned.tmp')
        out_file = os.path.join('execs', 'tmp', 'unrooted_tree.nwk')

        subprocess.call(["./execs/FastTree", "-out", out_file, in_file])

        rooted_tree = dendropy.Tree.get_from_path(out_file, schema='newick')
        rooted_tree.reroot_at_midpoint()
        rooted_tree.write_to_path(os.path.join('execs', 'tmp', "rooted_tree.nwk"), schema='newick')

        newick_rooted_file = open(os.path.join('execs', 'tmp', "rooted_tree.nwk"), 'r')
        read_edit_newick = newick_rooted_file.read()

        # The tree generated sometimes has the [&R] region - if not stripped it will throw error. Try except handles if the [&R] is not generated
        try:

            stripped_tree = read_edit_newick.strip('\[&R\] ')
            with open(os.path.join('execs', 'tmp', "rooted_tree.nwk"), 'w') as writeStrippedTree:
                writeStrippedTree.write('')
                writeStrippedTree.write(stripped_tree)

                with open(os.path.join('execs', 'tmp', "rooted_tree.nwk"), 'w') as writeStrippedTree:
                    writeStrippedTree.write('')
                    writeStrippedTree.write(stripped_tree)

        except AttributeError:
            pass

    def collectMultipleSequencingAlignment(self):

        MSA = handler.MSAfileHandler()

        in_file = MSA.getUnalignedTmpPath()
        out_file = MSA.getAlignedTmpPath()

        MSA.clearPreviousInput(in_file)
        MSA.clearPreviousInput(out_file)

        protein_description = self.proteinDescs
        protein_sequence = self.proteinSeqs
        common_names = self.commonNames

        self.msa = []

        with open(in_file, 'a') as msaprotein_writeFile:

            for i in range(len(protein_description)):


                protein = '\n'+ '>' + str(protein_description[i]) + str(common_names[i]).replace(' ', '_') + '\n' + str(protein_sequence[i])
                msaprotein_writeFile.write(protein)
                self.msa.append(str(protein_sequence[i]))

        clustalomega_cline = ClustalOmegaCommandline(cmd=os.path.join('execs', "clustalo-1.2.0"),
                                                     infile=in_file,
                                                     outfile=out_file, verbose=True, auto=True, force=True)
        clustalomega_cline()
        MSA.msa_FileCorrection()

    def constructTreeObj(self):

        tree = handler.treeOBjFileHandler()

        MSAfilePath = tree.getTreeInputPath()
        unrootedTreePath = tree.getTreeOutputPath()

        subprocess.call(["./execs/FastTree", "-out", MSAfilePath, unrootedTreePath])

        rootedTree = dendropy.Tree.get_from_path(unrootedTreePath, schema='newick')

        rootedTree.reroot_at_midpoint()

    def assignDomainColors(self, Domains):

        #Iterates through the domains taken from the database, removes the duplicates, and assigns a random color which will be used in the final figure
        rawDomainNames = [domain.keys() for domain in Domains]
        domainNames = []
        for domain in rawDomainNames:
            for eachDomain in domain:
                domainNames.append(eachDomain)

        domainNames = list(dict.fromkeys(domainNames))

        randomColor = RandomColor()
        domainNameColors = {domain:randomColor.generate()[0] for domain in domainNames}

        return domainNameColors

    def buildDomains(self):

        proteinAccessions = self.proteinAccessions
        parentDomains = self.parentDomains

        treeObj = handler.treeOBjFileHandler()
        treeObj.getRootedTreePath()

        with open(treeObj.getRootedTreePath()) as nwkTreeFile:
            nwkTree = nwkTreeFile.read()
            dt = Tree(nwkTree)

        dts = TreeStyle()
        dts.title.add_face(TextFace('PhyloPy - Protein Ortholog Finding Tool by Bryan Dighera: Protein Domains', fsize= 16,), column= 0)
        dts.allow_face_overlap = True
        dts.show_leaf_name = True
        dts.show_branch_support = True

        leafNames = dt.get_leaf_names()

        accessionDomains = {proteinAccessions[i]: parentDomains[i] for i in range(len(leafNames))}

        domainColors = self.assignDomainColors(accessionDomains.values())

        domainMotifs = []

        # The leaf names contain the description so the accession must be stripped in order to index with protein accessions from db
        leafAccessionExtracted = treeObj.getProteinAccession([leaf for leaf in leafNames])

        for leaf in leafAccessionExtracted:

            domains = accessionDomains[leaf]
            domainLen = len(accessionDomains[leaf])

            domainName = [list(domains.keys())[i] for i in range(domainLen)]
            domainStart = [int(list(domains.values())[i].split(':')[0].strip('<').strip('>'))/self.image_scaling for i in range(domainLen)]
            domainEnd = [int(list(domains.values())[i].split(':')[1].strip('<').strip('>'))/self.image_scaling for i in range(domainLen)]

            domainColor = [domainColors[domain] for domain in domains]

            leafMotifs = []

            for i in range(domainLen):

                leafMotifs.append({'startLocation': int(domainStart[i]),
                                                'endLocation': int(domainEnd[i]),
                                                'shape': '<>',
                                                'width': None,
                                                'height': 12,
                                                'background': 'Black',
                                                'foreground': domainColor[i],
                                                'name':domainName[i]})


            domainMotifs.append({'domains':leafMotifs, 'name':leaf})

        return {'Sequences':domainMotifs}

    def buildIntrons(self):

        proteinAccessions = self.proteinAccessions
        introns = self.Introns
        exonLengths = self.exonLengths

        treeObj = handler.treeOBjFileHandler()
        treeObj.getRootedTreePath()

        with open(treeObj.getRootedTreePath()) as nwkTreeFile:
            nwkTree = nwkTreeFile.read()
            t = Tree(nwkTree)
            nwkTreeFile.close()

        ts = TreeStyle()
        ts.title.add_face(
            TextFace('PhyloPy - Protein Ortholog Finding Tool by Bryan Dighera: Intron Location and Phases',
                     fsize=16, ), column=0)
        ts.allow_face_overlap = True
        ts.show_leaf_name = True
        ts.show_branch_support = True

        leafNames = t.get_leaf_names()

        accessionIntrons = {proteinAccessions[i]: [introns[i], exonLengths[i]] for i in range(len(leafNames))}

        dummyIntronMotif = [[0, 400, "-", None, 12, "Black", "Black", None]]
        intronMotifs = []


        #The leaf names contain the description so the accession must be stripped in order to index with protein accessions from db
        leafAccessionExtracted = treeObj.getProteinAccession([leaf for leaf in leafNames])
        lengths = []
        count = 1
        for leaf in leafAccessionExtracted:  # Corrects introns, and builds intron motifs


            intronPhases = accessionIntrons[leaf][0]
            exonLengths = accessionIntrons[leaf][1]


            if intronPhases and exonLengths != 'NONE':

                IPH = handler.ImageProcessingHandler()

                intronPhases = ast.literal_eval(intronPhases)
                exonLengths = ast.literal_eval(exonLengths)

                exonLength = [math.floor(int(exonLengths[i][0].split('-')[1]) / 3) for i in range(len(exonLengths))]
                lengths.append(exonLength)
                recordMotifs = []

                exonLocation, MSASeqlen = IPH.intron_fix(leaf, exonLength)

                for i in range(len(exonLocation)):

                    intronPhase = int(intronPhases[i]) % 3


                    if intronPhase == 0:
                        recordMotifs.append({'startLocation': exonLocation[i] - 1,
                                             'endLocation': exonLocation[i] + 1,
                                             'shape': '[]',
                                             'width': None,
                                             'height': 5,
                                             'background': 'Blue',
                                             'foreground': 'Blue',
                                             'realLocation': exonLocation[i]})

                    elif intronPhase == 1:
                        recordMotifs.append({'startLocation': exonLocation[i] - 1,
                                             'endLocation': exonLocation[i] + 1,
                                             'shape': '[]',
                                             'width': None,
                                             'height': 5,
                                             'background': 'Black',
                                             'foreground': 'Black',
                                             'realLocation': exonLocation[i]})


                    elif intronPhase == 2:
                        recordMotifs.append({'startLocation': exonLocation[i] - 1,
                                                           'endLocation': exonLocation[i] + 1,
                                                           'shape': '[]',
                                                           'width': None,
                                                           'height': 5,
                                                           'background': 'Grey',
                                                           'foreground': 'Grey',
                                                           'realLocation': exonLocation[i]})

                intronMotifs.append({'name':leaf, 'id':'000'+str(count), 'introns':recordMotifs})
                count +=1


        return {'Sequences':intronMotifs}

    def buildGenomicContext(self):

        #Set parent proteins and genomic context retrieved from SQLite into a variable
        proteinAccessions = self.proteinAccessions
        parentGC = self.GenomicContexts

        accessionGCdict = {proteinAccessions[i]:parentGC[i] for i in range(len(proteinAccessions))}

        #Largest length of the genomic context
        maxGClength = max([len(parentGC[i]) for i in range(len(parentGC))])


        #Strip all the domains from the entire datastructure so that all domains are stored in a single list
        GCdomains = itertools.chain(*[[parentGC[i][j]['domain'] for j in range(len(parentGC[i]))] for i in range(len(proteinAccessions))])
        GCdomains = list(itertools.chain(*GCdomains))
        GCdomains = list(itertools.chain(*GCdomains))[::2]


        #Assign each domain a color, as key value pair (dict), which will be assigned during motif construction
        rand_color = RandomColor()
        GCcolors = {GCdomains[i]:rand_color.generate()[0] for i in range(len(GCdomains))}

        treeObj = handler.treeOBjFileHandler()
        treeObj.getRootedTreePath()

        with open(treeObj.getRootedTreePath()) as nwkTreeFile:
            nwkTree = nwkTreeFile.read()
            t = Tree(nwkTree)
            nwkTreeFile.close()

        ts = TreeStyle()
        ts.title.add_face(
            TextFace('PhyloPy - Protein Ortholog Finding Tool by Bryan Dighera: Genomic Context',
                     fsize=16, ), column=0)
        ts.allow_face_overlap = True
        ts.show_leaf_name = True
        ts.show_branch_support = True

        leafNames = t.get_leaf_names()

        GCMotifs = {}

        # The leaf names contain the description so the accession must be stripped in order to index with protein accessions from db
        leafAccessionExtracted = treeObj.getProteinAccession([leaf for leaf in leafNames])

        geneIDs = {self.proteinAccessions[i]:self.geneID[i] for i in range(len(leafNames))}

        for j, leaf in enumerate(leafAccessionExtracted):

            try:
                # Function which takes the current leaf genomic context and aligns all working genes in same direction
                record = treeObj.fix_coding_direction(accessionGCdict[leaf], geneIDs[leaf])
                record = treeObj.fix_coding_alignment(record, geneIDs[leaf], maxGClength)

                coding_direction = [record[i]['coding_direction'] for i in range(len(record))]
                geneName = [record[i]['gene_name'] for i in range(len(record))]
                numberofDomains = [record[i]['domain'] for i in range(len(record))]

                flip = [record[i]['flip'] for i in range(len(record))]

                numberofGenes = len([math.floor(record[i]['img_start']) for i in range(len(record))])
                start_gene_location = [math.floor(record[i]['img_start']-self.stretch) for i in range(len(record))]
                end_gene_location = [math.floor(record[i]['img_end']-self.stretch) for i in range(len(record))]

                recordMotifs = {}

                for i in range(numberofGenes):

                    if i != None:

                        if coding_direction[i] == '-' and flip[i] == False:

                            if numberofDomains[i] != []:
                                recordMotifs[str(uuid.uuid4())] = {'startLocation': int(start_gene_location[i]),
                                                                   'endLocation': int(end_gene_location[i]),
                                                                   'shape': '[]',
                                                                   'width': None,
                                                                   'height': 12,
                                                                   'foreground': 'Black',
                                                                   'background': 'White',
                                                                   'direction': '<',
                                                                   'gene name': geneName[i],
                                                                   'domainSize': {j[0]: j[1] for j in
                                                                                  numberofDomains[i]},
                                                                   'domainColor': {j[0]: GCcolors[j[0]] for j in
                                                                                   numberofDomains[i]},
                                                                   'domainName': [j[0] for j in numberofDomains[i]]}



                            else:
                                recordMotifs[str(uuid.uuid4())] = {'startLocation': int(start_gene_location[i]),
                                                             'endLocation': int(end_gene_location[i]),
                                                             'shape': '[]',
                                                             'width': None,
                                                             'height': 12,
                                                             'foreground': 'Black',
                                                             'background': 'White',
                                                            'direction': '<',
                                                               'gene name':geneName[i],
                                                                   'domainSize': None,
                                                                   'domainColor': None,
                                                                   'domainName': None}

                        elif coding_direction[i] == '+' and flip[i] == True:

                            if numberofDomains[i] != []:
                                recordMotifs[str(uuid.uuid4())] = {'startLocation': int(start_gene_location[i]),
                                                                   'endLocation': int(end_gene_location[i]),
                                                                   'shape': '[]',
                                                                   'width': None,
                                                                   'height': 12,
                                                                   'foreground': 'Black',
                                                                   'background': 'White',
                                                                   'direction': '<',
                                                                   'gene name': geneName[i],
                                                                   'domainSize': {j[0]:j[1] for j in numberofDomains[i]},
                                                                   'domainColor':{j[0]:GCcolors[j[0]] for j in numberofDomains[i]},
                                                                   'domainName':[j[0] for j in numberofDomains[i]]}



                            else:
                                recordMotifs[str(uuid.uuid4())] = {'startLocation': int(start_gene_location[i]),
                                                             'endLocation': int(end_gene_location[i]),
                                                             'shape': '[]',
                                                             'width': None,
                                                             'height': 12,
                                                             'foreground': 'Black',
                                                             'background': 'White',
                                                            'direction': '<',
                                                               'gene name':geneName[i],
                                                                   'domainSize': None,
                                                                   'domainColor': None,
                                                                   'domainName': None}

                        elif coding_direction[i] == '-' and flip[i] == True:

                            if numberofDomains[i] != []:
                                recordMotifs[str(uuid.uuid4())] = {'startLocation': int(start_gene_location[i]),
                                                                   'endLocation': int(end_gene_location[i]),
                                                                   'shape': '[]',
                                                                   'width': None,
                                                                   'height': 12,
                                                                   'foreground': 'Black',
                                                                   'background': 'White',
                                                                   'direction': '<',
                                                                   'gene name': geneName[i],
                                                                   'domainSize': {j[0]:j[1] for j in numberofDomains[i]},
                                                                   'domainColor':{j[0]:GCcolors[j[0]] for j in numberofDomains[i]},
                                                                   'domainName':[j[0] for j in numberofDomains[i]]}



                            else:
                                recordMotifs[str(uuid.uuid4())] = {'startLocation': int(start_gene_location[i]),
                                                             'endLocation': int(end_gene_location[i]),
                                                             'shape': '[]',
                                                             'width': None,
                                                             'height': 12,
                                                             'foreground': 'Black',
                                                             'background': 'White',
                                                            'direction': '<',
                                                               'gene name':geneName[i],
                                                                   'domainSize': None,
                                                                   'domainColor': None,
                                                                   'domainName': None}

                        elif coding_direction[i] == '+' and flip[i] == False:

                            if numberofDomains[i] != []:
                                if numberofDomains[i] != []:
                                    recordMotifs[str(uuid.uuid4())] = {'startLocation': int(start_gene_location[i]),
                                                                       'endLocation': int(end_gene_location[i]),
                                                                       'shape': '[]',
                                                                       'width': None,
                                                                       'height': 12,
                                                                       'foreground': 'Black',
                                                                       'background': 'White',
                                                                       'direction': '<',
                                                                       'gene name': geneName[i],
                                                                       'domainSize': {j[0]: j[1] for j in
                                                                                      numberofDomains[i]},
                                                                       'domainColor': {j[0]: GCcolors[j[0]] for j in
                                                                                       numberofDomains[i]},
                                                                       'domainName': [j[0] for j in numberofDomains[i]]}



                                else:
                                    recordMotifs[str(uuid.uuid4())] = {'startLocation': int(start_gene_location[i]),
                                                                       'endLocation': int(end_gene_location[i]),
                                                                       'shape': '[]',
                                                                       'width': None,
                                                                       'height': 12,
                                                                       'foreground': 'Black',
                                                                       'background': 'White',
                                                                       'direction': '<',
                                                                       'gene name': geneName[i],
                                                                       'domainSize': None,
                                                                       'domainColor': None,
                                                                       'domainName': None}
                    else:
                        recordMotifs[leaf] = None

                GCMotifs[leaf]= recordMotifs

            except IndexError:
                print('Genomic Context Index Error at Sequence: %s' % leaf)
            except TypeError:
                pass

        return GCMotifs

