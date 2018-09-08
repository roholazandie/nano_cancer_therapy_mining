from Bio import Entrez

class EntrezSearch():

    def __init__(self):
        pass


    def search(self, query):
        Entrez.email = 'your.email@example.com'
        handle = Entrez.esearch(db='pubmed',
                                sort='relevance',
                                retmax='100000',
                                retmode='xml',
                                term=query + "[Title/Abstract]")
        results = Entrez.read(handle)
        return results


    def rule1_query(self, name1, name2):
        '''
        In case you would like to perform an exact phrase search (logical AND),
         you can do so by specifying double quotes in the search text
        :return:
        '''
        # TODO need to do something to handle multipart words! how to search them. like the first one in cancer.txt file
        phrase_to_search = str(name1) + " AND " + str(name2)
        results = self.search(phrase_to_search)
        pubmed_ids = results["IdList"]
        return pubmed_ids



if __name__ == "__main__":
    entrez_search = EntrezSearch()
    results = entrez_search.rule1_query("tumor", "nanoparticle")
    print(len(results))