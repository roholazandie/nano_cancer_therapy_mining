from config.nano_cancer_mining_configuration import NanoCancerConfiguration
from pymongo import MongoClient, TEXT
import pymongo


client = MongoClient()


class PubMedSearch(object):

    def __init__(self, config : NanoCancerConfiguration):
        database_name = config.database_config.name
        collection_name = config.database_config.collection_name
        self._db = client[database_name]
        self._collection_name = collection_name



    def rule1_query(self, name1, name2):
        '''
        In case you would like to perform an exact phrase search (logical AND),
         you can do so by specifying double quotes in the search text
        :return:
        '''
        phrase_to_search = "\"" + str(name1) + " " + str(name2) +"\""
        results = self._db[self._collection_name].find({"$text": {"$search": phrase_to_search}})
        for result in results:
            #print(result["article"]["journal"]["title"])
            if "abstract" in result["article"]:
                print(result["article"]["abstract"])


if __name__ == "__main__":
    pubmed_search = PubMedSearch()
    pubmed_search.rule1_query("cancer", "breast")