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

