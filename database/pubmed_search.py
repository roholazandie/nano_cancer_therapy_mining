from pymongo import MongoClient

from config.nano_cancer_mining_configuration import NanoCancerConfiguration

client = MongoClient()


class PubMedSearch(object):

    def __init__(self, config: NanoCancerConfiguration):
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
        #TODO need to do something to handle multipart words! how to search them. like the first one in cancer.txt file
        phrase_to_search = "\"" + str(name1) + " " + str(name2) + "\""
        results = self._db[self._collection_name].find({"$text": {"$search": phrase_to_search}})
        pubmed_ids = []
        for result in results:
            # print(result["article"]["journal"]["title"])
            if "abstract" in result["article"]:
                # print(result["article"]["abstract"])
                pubmed_ids.append(result["PMID"])

        return pubmed_ids


if __name__ == "__main__":
    import time
    t1 = time.time()
    nano_cancer_configs = NanoCancerConfiguration()
    pubmedsearch = PubMedSearch(nano_cancer_configs)
    pubmed_ids = pubmedsearch.rule1_query("cancer", "breast")
    print(pubmed_ids)
    t2 = time.time()

    print(t2-t1)