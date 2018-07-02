from pymongo import MongoClient, TEXT
import pymongo
client = MongoClient()
#print(client.address)
print(client.database_names())
#print(client.drop_database("pubmeddb"))

db = client["mytestdatabase"]
collection_name = "c2"
# db[collection_name].insert({"art":{"subject":"Joe owns a dog", "content":"Dogs are man's best friend"}, "likes": 60, "year":2015, "language":"english"})
# db[collection_name].insert({"art":{"subject": "Dogs eat cats and dog eats pigeons too", "content": "Cats are not evil"}, "likes": 30, "year": 2015,"language": "english"})
# db[collection_name].insert({"art":{"subject": "Cats eat rats", "content": "Rats do not cook food"}, "likes": 55, "year": 2014, "language": "english"})
# db[collection_name].insert({"art":{"subject": "Rats eat Joe", "content": "Joe ate a rat"}, "likes": 75, "year": 2014, "language": "english"})

#db[collection_name].create_index([("art.subject", TEXT)])

#a = db[collection_name].find({"$text": {"$search": "dogs"}}, {"score": {"$meta": "textScore"}}).sort([("score",{"$meta":"textScore"})])
a = db[collection_name].find({"$text": {"$search": "cats"}})
for m in a:
    print(m)


print(pymongo.__version__)