## How to populate database?

- prepare your xml files in a folder(the current version can handle both xml and xml.gz)

- write down the path of the xml files in config.yaml (raw_xml_dir)

- create a directory for your database somewhere that you're sure that you have enough space.

- open terminal and write:
```
mongod --dbpath=YOUR_DATABASE_PATH
```
for example:
```
mongod --dbpath="/home/rohola/mongodb-data/"
```

- it should be the same as the database root in config.yaml file.

- specify a name for your database and collection name.

- run the script database.mongodb to start populating the database.

- before start populating make sure that the database doesn't exist. if it exist drop it and then run the script again.

export PYTHONPATH="/home/rohola/codes/nano_cancer_therapy_mining"