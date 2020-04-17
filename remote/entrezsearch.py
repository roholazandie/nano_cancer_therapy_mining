from collections import OrderedDict

from Bio import Entrez
from Bio.Entrez import efetch
import xmltodict, json


class EntrezSearch():

    def __init__(self):
        pass


    def fetch(self, pubmedid):
        handle = efetch(db='pubmed', id=pubmedid, retmode='xml', rettype='abstract')
        return handle.read()

    def get_title(self, pubmedid):
        all_xml = self.fetch(pubmedid)
        article = xmltodict.parse(all_xml)
        try:
            title = article["PubmedArticleSet"]["PubmedArticle"]["MedlineCitation"]["Article"]["ArticleTitle"]
        except:
            title = ""
            print("ww")
        return title

    def get_abstract(self, pubmedid):
        all_xml = self.fetch(pubmedid)
        article = xmltodict.parse(all_xml)
        try:
            all_abstract = article["PubmedArticleSet"]["PubmedArticle"]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]
            if type(all_abstract)==str:
                print("ok")
                return all_abstract
            if type(all_abstract)==OrderedDict:
                print("ok")
                return dict(all_abstract)["#text"]
            abstract = " ".join([dict(abstract)["#text"] for abstract in all_abstract])
            print("ok")
        except:
            abstract = ""
            print(pubmedid)
            print("eww")
        return abstract

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
    #results = entrez_search.rule1_query("Electrochemical Biosensor", "cancer")
    #print(len(results))
    #print(entrez_search.fetch(24366930))
    #print(entrez_search.get_abstract(30224278))

    import ast
    pmids = ast.literal_eval(open("../dataset/CancerLiposome").read())
    with open("../outputs/cancer_liposome_results.txt", 'w') as file_writer:
        for pmid in pmids:
            title = entrez_search.get_title(pmid)
            abstract = entrez_search.get_abstract(pmid)
            file_writer.write(pmid+"\t"+str(title)+"\t"+str(abstract)+"\n")