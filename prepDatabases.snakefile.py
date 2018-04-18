import gzip
import time
import csv
from Bio import Entrez
from Bio.Entrez.Parser import ValidationError
from Bio import SeqIO

configfile: "config.json"


rule db_all:
    input: "%(dbFolder)s/58S_derep.fasta" % config, "%(dbFolder)s/58S_tax.tsv" % config

rule db_getRfamFile:
    output: "%(dbFolder)s/RF00002_%(rfam_version)s.fasta.gz" % config
    shell:
        "wget -O %(dbFolder)s/RF00002_%(rfam_version)s.fasta.gz \"ftp://ftp.ebi.ac.uk/pub/databases/Rfam/%(rfam_version)s/fasta_files/RF00002.fa.gz\"" % config
        
rule db_makeRfamFiles:
    input: fasta="%(dbFolder)s/RF00002_%(rfam_version)s.fasta.gz" % config
    output: fasta="%(dbFolder)s/RF00002_%(rfam_version)s_dna.fasta" % config, tax="%(dbFolder)s/RF00002_%(rfam_version)s_tax.tsv" % config
    log: "%(dbFolder)s/logs/rfam_tax_log.txt" % config
    run:
        ranks = ["superkingdom","kingdom", "phylum", "class", "order", "family", "genus", "species"]
        if not config["email"] or not "@" in config["email"]:
            raise ValueError("'%s' is not a valid email. Set email in config file for NCBI querries." % config["email"])
        nuc2tax = NuclId2TaxIdMap(config["email"], cachePath="%(dbFolder)s/nuc2tax.csv" % config, retry=3)
        tax2lin = LineageMap(config["email"], cachePath="%(dbFolder)s/tax2lin.csv" % config, retry=3)
        noId = 0
        noLin = 0
        with open(output.fasta, "w") as out, open(output.tax, "w") as taxOut, open(log[0], "w") as logFile:
            for rec in SeqIO.parse(gzip.open(input.fasta, "rt"), "fasta"):
                
                
                sId, reg = rec.id.split("/")
                rec.id="%s|%s" % (sId, reg)
                try:
                    tId = nuc2tax[sId]
                except:
                    logFile.write("no ncbi tax ID\t%s\n" % sId)
                    noId += 1
                    continue
                try:
                    lin = tax2lin[tId]
                except:
                    logFile.write("no ncbi lineage\t%s\n" % Id)
                    noLin += 1
                    continue
                linStr = ";".join([l[2] for l in lin if l[0] in ranks])
                taxOut.write("%s\t%s\n" % (rec.id, linStr))
                out.write(rec.format("fasta"))
        sys.stderr.write("No taxon ID: %i, no lineage: %i\n" % (noId, noLin))
        nuc2tax.save()
        tax2lin.save()


rule db_getUniteFile:
    output: "%(dbFolder)s/sh_general_release_dynamic_%(unite_version)s.fasta" % config
    shell:
        "cd %(dbFolder)s;"\
        "wget %(uniteUrl)s;"\
        "unzip *.zip;"\
        "rm *.zip" % config

rule db_makeUniteFiles:
    input: "%(dbFolder)s/sh_general_release_dynamic_%(unite_version)s.fasta" % config
    output: fasta="%(dbFolder)s/unite_%(unite_version)s.fasta" % config, tax="%(dbFolder)s/unite_%(unite_version)s.tsv" % config, sh2gId="%(dbFolder)s/unite_%(unite_version)s_gIds.tsv" % config
    run:
        with open(output.fasta, "w") as fasta, open(output.tax, "w") as tax, open(output.sh2gId, "w") as sh2gId:
            for line in open(input[0]):
                if line[0] == ">":
                    name, gId, sh, typ, lin = line[1:].strip("\n").split("|")
                    fasta.write(">%s %s\n" % (sh, name))
                    tax.write("%s\t%s\n" % (sh, lin))
                    sh2gId.write("%s\t%s\n" % (sh, gId))
                else:
                    fasta.write(line)

rule db_extract58S:
    input: "%(dbFolder)s/unite_%(unite_version)s.fasta" % config
    output: "%(dbFolder)s/ITSx/unite_%(unite_version)s.5_8S.fasta" % config
    log: "%(dbFolder)s/logs/db_itsx.log"
    threads: 6
    shell:
        "%(itsx)s -t . -i {input} -o %(dbFolder)s/ITSx/unite_%(unite_version)s --save_regions 5.8S --cpu {threads} --graphical F &> {log}" % config

rule cat58S:
    input: "%(dbFolder)s/ITSx/unite_%(unite_version)s.5_8S.fasta" % config, "%(dbFolder)s/RF00002_%(rfam_version)s_dna.fasta" % config
    output: "%(dbFolder)s/58S.fasta" % config
    shell:
        "cat {input} > {output}"

#rule derepById:
#    input: fasta="%(dbFolder)s/58S.fasta" % config, sh2gId="%(dbFolder)s/unite_%(unite_version)s_gIds.tsv" % config
#    output: "%(dbFolder)s/58S_derepById.fasta" % config
#    log: "%(dbFolder)s/logs/derepById58S.log" % config
#    run:
#        sh2gId = {}
#        gId2sh = {}
#        for line in open(input.sh2gId):
#            sh, gId = line.strip("\n").split("\t")
#            sh2gId[sh] = gId
#            gId2sh[gId] = sh
#        with open(output[0], "w") as out:
#            for rec in SeqIO.parse(open(input.fasta), "fasta"):
#                if rec.id[:2] != "SH":
#                    if rec.id.split(".", 1)[0] not in gId2sh:
#                        #this genbank sequence is not represented by a UNITE sequence
#                        out.write(rec.format("fasta"))
#                else:
#                    out.write(rec.format("fasta"))

rule derep:
    input: "%(dbFolder)s/58S.fasta" % config
    output: fasta="%(dbFolder)s/58S_derep.fasta" % config, uc="%(dbFolder)s/58S_derep.uc.txt" % config
    log: "%(dbFolder)s/logs/derep58S.log" % config
    run:
        shell("%(vsearch)s --derep_fulllength {input} --output {output.fasta} --uc {output.uc} --sizeout --log {log} &> /dev/null" % config)

rule createTax:
    input: uc="%(dbFolder)s/58S_derep.uc.txt" % config, uTax="%(dbFolder)s/unite_%(unite_version)s.tsv" % config, rTax="%(dbFolder)s/RF00002_%(rfam_version)s_tax.tsv" % config
    output: "%(dbFolder)s/58S_tax.tsv" % config
    run:
        clu={}
        for line in open(input.uc):
            lType, cNum, sLen, ident, strand, _1, _2, aln, query, target = line.strip("\n").split("\t")
            if lType == "S":
                clu[query] = [query]
            elif lType == "H":
                clu[target].append(query)
            elif lType == "C":
                assert int(sLen) == len(clu[query])
            else:
                raise ValueError("Unknown record type: %s" % lType)
        uTax = {}
        for line in open(input.uTax):
            tId, tLin = line.strip().split("\t")
            uTax[tId] = tLin
        rTax = {}
        for line in open(input.rTax):
            tId, tLin = line.strip().split("\t")
            rTax[tId] = tLin
        with open(output[0], "w") as out:
            for rep, memList in clu.items():
                linStrs = []
                for m in memList:
                    if m[:2] == "SH":
                        linStrs.append(uTax[m.split("|")[0]])
                    else:
                        lin = rTax[m]
                        if lin.split(";")[1] == "Fungi":
                            #ignore everything RFAM has to say about fungi taxonomy
                            linStrs.append("Eukaryota;Fungi")
                        else:
                            linStrs.append(lin)
                lcaLin = lca(linStrs, 0.95)
                out.write("%s\t%s\n" % (rep, lcaLin))

rule db_creat58SIndex:
    input: "%(dbFolder)s/58S_derep.fasta" % config
    output: touch("%(dbFolder)s/58S_derep.fasta.lambdaIndexCreated" % config)
    threads: 6
    shell:
        "%(lambdaFolder)s/lambda_indexer -d {input} -p blastn -t {threads}" % config

rule db_creatUniteIndex:
    input: "%(dbFolder)s/unite_%(unite_version)s.fasta" % config
    output: touch("%(dbFolder)s/unite_%(unite_version)s.fasta.lambdaIndexCreated" % config)
    threads: 6
    shell:
        "%(lambdaFolder)s/lambda_indexer -d {input} -p blastn -t {threads}" % config


class NcbiMap(dict):
    """Base class for a dictonary relying on the BioPython NCBI interface for 
    assignments.
    
    This is an abstract base class for a dictonary that uses the BioPython NCBI 
    interface to represent relations between certain NCBI governed objects like
    genes, taxonomic objects or genomes.
    Key requests will be looked up from a dictonary. If the key is not available
    a request to NCBI will be send via the BioPython interface. The result of 
    this lookup will be stored in the local dictonary. The detaileds of the 
    request have to be implemeneted by inheriting calsses.
    The local dictonary can be saved to the hard drive as csv file and loaded
    when an object is constructed.
    """
    def __init__(self, email, indict={}, cachePath=None, retry=0, useCache=True):
        """Constuctor for mapping object.
        
        A dictonary of already known assignments can be passed via the indict
        parameter. The cachePath is the path were the save and load function
        will search for a csv file. retry gives the number of times a request 
        should be send to NCBI after again before giving up (0 means only send 
        once).
        """
        dict.__init__(self, indict)
        Entrez.email = email
        self.useCache = useCache
        if useCache:
            if not cachePath:
                #if no cache path was given default to the class name
                cachePath = self.__class__.__name__+".csv"
            self.cachePath = cachePath
            try:
                self.load()
            except IOError:
                pass
                #this happens if no cache was saved previously
        self.retry = retry
       
        
    def __getitem__(self, key):
        if key not in self:
            for tries in range(self.retry+1):
                try:
                    handle = self.requestFunction(key)
                    response = Entrez.read(handle, validate=False)
                except ValidationError as e:
                    sys.stderr.write("Problem while parsing response for key '%s'. "
                                     "This might be due to outdatetd DTD files in "
                                     "Biopython. Try updating from http://www.ncbi."
                                     "nlm.nih.gov/data_specs/dtd/" % key)
                    raise e
                except Exception as e:
                    if tries == self.retry:
                        raise KeyError("Problem with key: '%s'. "
                                       "It raised an error: %s" 
                                       % (key, str(e)))
                    else:
                        time.sleep(1)
                        sys.stderr.write("Problem with key: '%s'. "
                                         "It raised an error: %s\n "
                                         "Will try again...\n" % (key, str(e)) )
                else:
                    #if no exception occured no retry is needed
                    break
            self.readResponse(response, key)
        return dict.__getitem__(self, key)
            
    def requestFunction(key):
        raise NotImplemented("This base class does not implement any request. "
                             "Use a more specific subclass.")
            
    def load(self):
        if not self.useCache:
            raise CacheNotUsedError()
        with open(self.cachePath, "r") as inFile:
            for row in csv.reader(inFile, delimiter=",", quotechar="\""):
                self[row[0]] = row[1]
    
    def save(self):
        if not self.useCache:
            raise CacheNotUsedError()
        with open(self.cachePath, "w") as out:
            writer = csv.writer(out, delimiter=",", quotechar="\"", 
                                quoting=csv.QUOTE_MINIMAL)
            for key, value in self.items():
                writer.writerow([str(key), str(value)])

class LineageMap(NcbiMap):
    """Map NCBI taxonomy IDs to full lineage information from NCBI taxonomy.
    
    Returns a list of tuples of the form (<Rank>, <Taxonomy Node ID>, 
    <Taxonomy Node Name>).
    """
    
    def requestFunction(self, key):
        return Entrez.efetch(db="taxonomy", id=str(key))
    
    def readResponse(self, resp, key):
        if len(resp) == 0:
            raise KeyError("No result for lineage of species '%s'" % key)
        m = []
        if "LineageEx" not in resp[0]:
            if resp[0]["ScientificName"] != "root":
                raise ValueError("Wired NCBi reponse without lineage info: %s" 
                                 % str(resp))
        else:
            for r in resp[0]["LineageEx"]:
                m.append((r["Rank"], r["TaxId"], r["ScientificName"]))
        m.append((resp[0]["Rank"], resp[0]["TaxId"], resp[0]["ScientificName"]))
        self[key] = m
        
    def save(self):
        if not self.useCache:
            raise CacheNotUsedError()
        tab = []
        for tax, m in self.items():
            for rank, lTax, lName in m:
                tab.append([tax, rank, lTax, '"%s"' % lName])
        with open(self.cachePath, "w") as out:
            for row in tab:
                out.write(",".join([str(field) for field in row])+"\n")
    
    def load(self):
        if not self.useCache:
            raise CacheNotUsedError()
        for row in csv.reader(open(self.cachePath, "r")):
            tax, rank, lTax, lName = row
            if tax not in self:
                self[tax] = []
            self[tax].append((rank, lTax, lName))

class NuclId2TaxIdMap(NcbiMap):
    """Map NCBI nucleotide ID to the taxonomy ID of the species they come from.
    
    """
    
    def requestFunction(self, key):
        return Entrez.efetch(db="nucleotide", id=str(key), retmode="xml")
    
    def readResponse(self, resp, key):
        if len(resp) < 1:
            raise KeyError("'%s' is not in the dictionary. "
                           "NCBI response did not contain taxonomy information.")
        if len(resp) > 1:
            self[key] = None
            raise ValueError("Problem with key: %s. "
                             "It got multiple answers." % key)
        
        #if there is only one feature for some reason 
        # the feature list is not actually a list
        #here is a work around:
        try:
            feature_list = resp[0]["GBSeq_feature-table"]
            for feature in feature_list:
                if feature["GBFeature_key"] == "source":
                    for qual in feature["GBFeature_quals"]:
                        if qual["GBQualifier_name"] == "db_xref":
                            match = re.match("taxon:(\d+)", 
                                             qual["GBQualifier_value"])
                            if match:
                                self[key] = match.group(1)
                                return
        except Exception as e:
            sys.stderr.write(str(resp))
            raise KeyError("'%s' is not in the dictonary. "
                           "Reading NCBI response caused an exception."
                           "NCBI response was:\n%s" % (key, str(resp)))
        raise KeyError("'%s' is not in the dictonary. "
                       "NCBI response did not contain taxonomy information. "
                       "NCBI response was:\n%s" % (key, str(resp)[:100]))

