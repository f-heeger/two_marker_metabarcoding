class TaxTreeNode(object):
    def __init__(self, name):
        self.name = name
        self.children = {}
        self.reads = {}
        self.otus = {}
        
    def loadSubTree(self, inDict, sample):
        for key, value in inDict.items():
            if key == "$":
                try:
                    self.reads[sample] += sum(value)
                except KeyError:
                    self.reads[sample] = sum(value)
                try:
                    self.otus[sample] += len(value)
                except KeyError:
                    self.otus[sample] = len(value)
            else:
                if key not in self.children:
                    newChild = TaxTreeNode(key)
                    self.children[key] = newChild
                self.children[key].loadSubTree(value, sample)
    
            
    def kronaXml(self, order, indent=""):
        lines=[]
        lines.append('%s<reads>' % indent)
        for sample in order:
            lines.append('%s    <val sample="%s">%i</val>' % (indent, sample, 
                                                              self.reads.get(sample, 0)))
        lines.append('%s</reads>' % indent)
        lines.append('%s<otus>' % indent)
        for sample in order:
            lines.append('%s    <val sample="%s">%i</val>' % (indent, sample, 
                                                              self.otus.get(sample, 0)))
        lines.append('%s</otus>' % indent)
        for name, child in self.children.items():
            lines.append('%s<node name="%s">' % (indent, name.split(";")[-1]))
            lines.extend(child.kronaXml(order, indent+"    "))
            lines.append('%s</node>' % indent)
        return lines
