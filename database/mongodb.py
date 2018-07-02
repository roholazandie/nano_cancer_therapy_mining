import time
from collections import Counter
import pymongo
from pymongo import MongoClient, IndexModel
from database.xmlread import XMLRead

class MongoDB(object):
    def __init__(self, database_name):
        client = MongoClient()
        self.db = client[database_name]

    def populate_db(self, import_dir_name, collection_name):
        """
        get each document in import_dir_name dirctory convert them to json an
        insert them into the collection_name collection of database
        :param import_dir_name: directory name contains all the xml files
        :param collection_name: the name of collection want to store the files
        """
        i = 0
        for article_json in XMLRead(import_dir_name):
            self.db[collection_name].insert_one(article_json)
            i+=1
            if i%10000==0:
                print(i)

    def read(self, collection_name, limit=None):
        """
        read limit number of the documents in the collection_name collection
        :param collection_name: the collection that need to be searched
        :param limit: number of documents to be retrieved
        :return:
        """
        documents = []
        for document in self.db[collection_name].find().limit(limit):
            documents.append(document)
        return documents

    def query(self, collection_name, query):
        """
        query the collection_name in database with query string(as json)
        :param collection_name:
        :param query:
        :return:
        """
        documents = []
        for document in self.db[collection_name].find(query):
            documents.append(document)

        return documents

    def drop(self, collection_name):
        """
        drop collection_name collection from database
        :param collection_name: the collection that need to be dropped
        """
        self.db[collection_name].drop()

    def count(self, collection_name, query):
        print(self.db[collection_name].find().count())

    def get_all_countries(self, collection_name):
        """
        search in all the collection and find countries
        finally remove duplicates to get the list of unique countries
        :param collection_name:
        :return:
        """
        countries = []
        for document in self.db[collection_name].find():
            if "country" in document:
                countries.append(document["country"])
        return set(countries)

    def top10_journals(self, collection_name):
        """
        get top10 journals. to do this it collects all the
        journal names and then use Counter to get top 10 journals
        :param collection_name:
        :return:
        """
        journals = []
        for document in self.db[collection_name].find():
            if "journal" in document["article"]:
                journals.append(document["article"]["journal"]["title"])
        counter = Counter(journals)
        return counter.most_common(10)

    def simple_query(self, collection_name, countries, date_from, date_to, top10_journals, top10_authors):
        """
        extracts abstracts and corresponding dates for a query
        that contains countries, date from, date to, top10 journals and top10 authors
        :param collection_name: the name of the collection that need to be searched
        :param countries: list of countries, if empty it searches all the countries
        :param date_from: the date from of articles, if empty it searches from first date
        :param date_to: the date to of articles, if empty it searches to the end date
        :param top10_journals: whether search in top10 journals only or not
        :param top10_authors: whether search in top10 authors or not
        :return: articles abstract and corresponding dates
        """
        abstracts = []
        dates = []

        if top10_journals:
            journal_query = []
            for journal in top10_journals:
                journal_query.append({"article.journal.title": journal})

        if top10_authors:
            author_query = []
            for author in top10_authors:
                lastname, forename, initials = author.split(' ')
                author_query.append({"article.authorlist": {"$elemMatch": {
                                                                    "lastname": lastname,
                                                                    "forename": forename,
                                                                    "initials": initials
                                                          }}})
        query = {}
        if countries:
            query["country"] = {"$in": countries}

        if date_from and date_to:
            query["datecreated"] = {"$gte": date_from, "$lt": date_to}


        if top10_authors:
            query["$and"] = []
            query["$and"].append({"$or": author_query})

        if top10_journals:
            if "$and" not in query:
                query["$and"] = []
            query["$and"].append({"$or": journal_query})

        #print(self.db[collection_name].find(query).count())
        total = 1000000#self.db[collection_name].find(query).count()
        print(total)
        for document in self.db[collection_name].find(query).limit(total):
            if "abstract" in document["article"]:
                abstracts.append(document["article"]["abstract"])
                dates.append(document["datecreated"])
        return dates, abstracts

    def simple_frequency_query(self, collection_name, countries, date_from, date_to, top10_journals, top10_authors):
        abstracts = []
        dates = []

        if top10_journals:
            journal_query = []
            for journal in top10_journals:
                journal_query.append({"article.journal.title": journal})

        if top10_authors:
            author_query = []
            for author in top10_authors:
                lastname, forename, initials = author.split(' ')
                author_query.append({"article.authorlist": {"$elemMatch": {
                    "lastname": lastname,
                    "forename": forename,
                    "initials": initials
                }}})
        query = {}
        if countries:
            query["country"] = {"$in": countries}

        if date_from and date_to:
            query["datecreated"] = {"$gte": date_from, "$lt": date_to}

        if top10_authors:
            query["$and"] = []
            query["$and"].append({"$or": author_query})

        if top10_journals:
            if "$and" not in query:
                query["$and"] = []
            query["$and"].append({"$or": journal_query})

        # query = {"country": {"$in": countries},
        #          "datecreated": {"$gte": date_from, "$lt": date_to},
        #          "$and": [
        #                     {"$or": journal_query},
        #                     {"$or": author_query}
        #             ]
        #          }
        # print(self.db[collection_name].find(query).count())
        total = 1000#self.db[collection_name].find(query).count()
        print(total)
        i = 0
        for document in self.db[collection_name].find(query).limit(total):
            i += 1
            if "abstract" in document["article"]:
                # print(document["article"]["abstract"])
                abstracts.append(document["article"]["abstract"])
                # print(document["datecreated"])
        return abstracts

    def create_index(self, collection_name, index_name):
        self.db[collection_name].create_index([(index_name, pymongo.ASCENDING)])


    def create_text_index(self, collection_name, fields):
        #TODO doesn't work with the latest version
        index_model = [(field, pymongo.TEXT) for field in fields]
        index = IndexModel(index_model, name="text_indices")
        self.db[collection_name].create_index(index)



    def drop_index(self, collection_name, index_name):
        self.db[collection_name].drop_index(index_name)


    def test_query(self, collection_name):
        #a = self.db[collection_name].find({"$text": {"$search": "macrophages"}}, {"score": {"$meta": "textScore"}}).sort([("score", {"$meta": "textScore"})])
        #a = self.db[collection_name].find({"$text": {"$search": "macrophages"}}).count()
        a = self.db[collection_name].find({}).count()

        return a


if __name__ == '__main__':
    t1 = time.time()

    #TODO adding a page for population database
    #dirname = "/home/rohola/Desktop/ftp"
    dirname = "/media/rohola/09183702400/home/Desktop/dataset"
    mongodb = MongoDB("pubmeddb")
    # mongodb.populate_db(dirname, "articles")

    #mongodb.read("articles")
    #mongodb.drop("articles")
    #q = {"article.journal.volume": '23'}
    #mongodb.query("articles", q)
    #mongodb.count("articles", q)
    #countries = mongodb.get_all_countries("articles")
    #print(countries)
    #t2 = time.time()
    #print(t2 - t1)


    # top10_journals = [u'The Journal of biological chemistry', u'Science (New York, N.Y.)', u'PloS one', u'Lancet (London, England)',
    #  u'Proceedings of the National Academy of Sciences of the United States of America', u'Nature',
    #  u'Biochimica et biophysica acta', u'British medical journal',
    #  u'Biochemical and biophysical research communications', u'Physical review letters']
    # top10_journals = None
    # top10_authors = ['Wang Wei W', 'Zhang Wei W', 'Suzuki T T', 'Wang Y Y', 'Suzuki K K', 'Nakamura T T', 'Tanaka K K', 'Takahashi M M', 'Li Wei W', 'Takahashi K K']
    # top10_authors = None
    # countries = None#["United States", "China"]
    # date_from = None#datetime(1970, 1, 1)
    # date_to = None#datetime(2010, 1, 1)
    # t1 = time.time()
    # dates, abstracts = mongodb.simple_query("articles",
    #                                          countries,
    #                                           date_from,
    #                                          date_to,
    #                                          top10_journals,
    #                                           top10_authors)
    #
    # print(abstracts)


    #mongodb.drop_index(collection_name="articles", index_name="article.abstract_text")
    mongodb.create_text_index(collection_name="articles", fields=["article.abstract"])
    #a = mongodb.test_query(collection_name="articles")
    # print(a)
    # for m in a:
    #     print(m)

    t2 = time.time()
    print(t2-t1)
    print(pymongo.__version__)
