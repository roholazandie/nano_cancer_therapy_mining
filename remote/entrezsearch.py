from Bio import Entrez
from Bio.Entrez import efetch


class EntrezSearch():

    def __init__(self):
        pass


    def fetch(self, pubmedid):
        handle = efetch(db='pubmed', id=pubmedid, retmode='xml', rettype='abstract')
        return handle.read()

    def search(self, query):
        Entrez.email = 'your.email@example.com'
        handle = Entrez.esearch(db='pubmed',
                                sort='relevance',
                                retmax='100000',
                                retmode='xml',
                                term=query)
        results = Entrez.read(handle)
        return results


    def rule1_query(self, name1, name2):
        '''
        In case you would like to perform an exact phrase search (logical AND),
         you can do so by specifying double quotes in the search text
        :return:
        '''
        # TODO need to do something to handle multipart words! how to search them. like the first one in cancer.txt file
        phrase_to_search = "(\""+str(name1)+"\"[Title/Abstract])" + " AND " + str(name2)+ "[Title/Abstract]"
        results = self.search(phrase_to_search)
        pubmed_ids = results["IdList"]
        return pubmed_ids



if __name__ == "__main__":
    entrez_search = EntrezSearch()
    results = entrez_search.rule1_query("Electrochemical Biosensor", "cancer")
    print(len(results))
    #entrez_search.fetch(24366930)