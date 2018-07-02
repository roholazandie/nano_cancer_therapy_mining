from pymongo import MongoClient, TEXT
import pymongo
client = MongoClient()

db = client["pubmeddb"]
collection_name = "articles"


def simple_query():
    results = db[collection_name].find({"$text": {"$search": "cat"}})
    for result in results:
        print(result["article"]["journal"]["title"])
        #print(result["article"]["abstract"])



def or_query():
    '''
    By default, the phrase search makes an OR search on all the specified keywords
    :return:
    '''
    words = ["cat", "bird"]
    search_words = " ".join(words)
    results = db[collection_name].find({"$text": {"$search": search_words}})
    for result in results:
        print(result["article"]["journal"]["title"])


def and_query():
    '''
    In case you would like to perform an exact phrase search (logical AND),
     you can do so by specifying double quotes in the search text
    :return:
    '''
    results = db[collection_name].find({"$text": {"$search": "\"cell cancer\""}})
    for result in results:
        print(result["article"]["journal"]["title"])


def negation_query():
    '''
    Prefixing a search keyword with – (minus sign) excludes all the documents that contain the negated term
    For example, try searching for any document which contains the keyword cell but does not contain cancer using the following query
    :return:
    '''
    results = db[collection_name].find({"$text": {"$search": "cell -cancer"}})
    for result in results:
        print(result["article"]["journal"]["title"])



def advanced_search():
    '''
    This search looks for documents which contains cell and cancer but not body
    :return:
    '''
    results = db[collection_name].find({"$text": {"$search": "\"cell cancer\" -body"}})
    for result in results:
        print(result["article"]["journal"]["title"])



def most_relevant_results():
    '''
    Since we are running a text search, we are also interested in getting some statistics about how relevant the resultant documents are.
    For this purpose, we will use the { $meta: "textScore" } expression, which provides information on the processing of the $text operator.
     We will also sort the documents by their textScore using the sort command. A higher textScore indicates a more relevant match.
    :return:
    '''
    results = db[collection_name].find({"$text": {"$search": "dog"}}, {"$score": {"$meta": "textScore"}}).sort([("$score",{"$meta":"textScore"})]).limit(1000)
    for result in results:
        print(result["article"]["abstract"])
        break


if __name__ == "__main__":
    simple_query()
    #or_query()
    #and_query()
    #negation_query()
    #advanced_search()
    #most_relevant_results()