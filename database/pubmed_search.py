from pymongo import MongoClient, TEXT
import pymongo
client = MongoClient()

#db = client["pubmeddb"]
#collection_name = "articles"


class PubMedSearch(object):

    def __init__(self):
        self._db = client["pubmeddb"]
        self._collection_name = "articles"



    def and_query(self, name1, name2):
        '''
        In case you would like to perform an exact phrase search (logical AND),
         you can do so by specifying double quotes in the search text
        :return:
        '''
        phrase_to_search = "\"" + str(name1) + " " + str(name2) +"\""
        results = self._db[self._collection_name].find({"$text": {"$search": phrase_to_search}})
        for result in results:
            print(result["article"]["journal"]["title"])