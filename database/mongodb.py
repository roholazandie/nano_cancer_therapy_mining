from collections import Counter
import pymongo
from pymongo import MongoClient, IndexModel
from database.xmlread import XMLRead
from config.nano_cancer_mining_configuration import NanoCancerConfiguration

class MongoDB(object):
    def __init__(self, config : NanoCancerConfiguration):
        client = MongoClient()
        self._config = config
        self.db = client[config.database_config.name]

    def populate_db(self):
        """
        get each document in import_dir_name dirctory convert them to json an
        insert them into the collection_name collection of database
        :param import_dir_name: directory name contains all the xml files
        :param collection_name: the name of collection want to store the files
        """
        i = 0
        for article_json in XMLRead(self._config.file_config.raw_xml_dir):
            if 'abstract' in article_json["article"]:
                self.db[self._config.database_config.collection_name].insert_one(article_json)
                i+=1
                if i%10000==0:
                    print(i)


    def read(self, limit=None):
        """
        read limit number of the documents in the collection_name collection
        :param collection_name: the collection that need to be searched
        :param limit: number of documents to be retrieved
        :return:
        """
        documents = []
        for document in self.db[self._config.database_config.collection_name].find().limit(limit):
            documents.append(document)
        return documents

    def query(self, query):
        """
        query the collection_name in database with query string(as json)
        :param collection_name:
        :param query:
        :return:
        """
        documents = []
        for document in self.db[self._config.database_config.collection_name].find(query):
            documents.append(document)

        return documents


    def drop(self, ):
        """
        drop collection_name collection from database
        :param collection_name: the collection that need to be dropped
        """
        self.db[self._config.database_config.collection_name].drop()


    def count(self):
        print(self.db[self._config.database_config.collection_name].find().count())


    def get_all_countries(self):
        """
        search in all the collection and find countries
        finally remove duplicates to get the list of unique countries
        :param collection_name:
        :return:
        """
        countries = []
        for document in self.db[self._config.database_config.collection_name].find():
            if "country" in document:
                countries.append(document["country"])
        return set(countries)

    def top10_journals(self):
        """
        get top10 journals. to do this it collects all the
        journal names and then use Counter to get top 10 journals
        :param collection_name:
        :return:
        """
        journals = []
        for document in self.db[self._config.database_config.collection_name].find():
            if "journal" in document["article"]:
                journals.append(document["article"]["journal"]["title"])
        counter = Counter(journals)
        return counter.most_common(10)


    def create_index(self, collection_name, index_name):
        self.db[collection_name].create_index([(index_name, pymongo.ASCENDING)])


    def create_text_index(self, collection_name, fields):
        #TODO doesn't work with the latest version
        index_model = [(field, pymongo.TEXT) for field in fields]
        index = IndexModel(index_model, name="text_indices")
        self.db[collection_name].create_index(index)



    def drop_index(self, collection_name, index_name):
        self.db[collection_name].drop_index(index_name)



if __name__ == '__main__':
    nano_cancer_configs = NanoCancerConfiguration()
    mongodb = MongoDB(nano_cancer_configs)
    mongodb.populate_db()
