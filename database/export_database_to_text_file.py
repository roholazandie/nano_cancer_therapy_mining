from pymongo import MongoClient, TEXT
import spacy
from tqdm import tqdm

client = MongoClient()
db = client["pubmeddb"]
collection_name = "articles"

nlp = spacy.load('en')

def test_db_connection():
    try:
        print(client.server_info())
        print("Connection exist!")
    except:
        raise Exception("ERROR: No connection")



def export(output):
    i = 0
    with open(output, 'w') as file_writer:
        results = db[collection_name].find()
        for result in results:
            i+=1
            if "abstract" in result["article"]:
                file_writer.write(result["article"]["abstract"]+"\n")

            if i%100000==0:
                print(i)



def read_keywords():
    lines = open("dataset/keywords.txt", encoding="utf-8").readlines()
    keywords = [line.rstrip() for line in lines]
    return keywords

def read_all_keywords():
    '''
    this is for keywords from two consecutive files. the keywords are the old one
    and output is the new one. now we add output to all keywords instead of
    using superposition of words
    :return:
    '''
    lines1 = open("dataset/keywords.txt", encoding="utf-8").readlines()
    lines2 = open("dataset/output.txt", encoding="utf-8").readlines()
    lines = lines1+lines2

    keywords = []
    for line in lines:
        words = preprocess_text(line.rstrip())
        keywords.append(" ".join(words))

    keywords = list(set(keywords))
    return keywords


def clean_text(input_file, output_file):
    nlp = spacy.load('en')
    #keywords = read_keywords()
    keywords = read_all_keywords()
    with open(output_file, 'w') as file_writer:
        with open(input_file) as file_reader:
            for i, sen in tqdm(enumerate(file_reader), total=16500000):
                en_abs = nlp(sen)
                words = [token.lemma_ for token in en_abs if not (token.pos_ == "NUM"
                                                                  or token.pos_ == "SYM"
                                                                  or token.pos_ == "PUNCT"
                                                                  or token.pos_ == "CCONJ"
                                                                  or token.like_num
                                                                  or token.is_space
                                                                  or token.is_punct)]

                article = ' '.join(words)
                for keyword in keywords:
                    if keyword in article:
                        print(keyword)
                        article = article.replace(keyword, keyword.replace(" ", "_"))

                file_writer.write(article+"\n")



def preprocess_text(text):
    en_abs = nlp(text)
    words = [token.lemma_ for token in en_abs if not (token.pos_ == "NUM"
                                                      or token.pos_ == "SYM"
                                                      or token.pos_ == "PUNCT"
                                                      or token.pos_ == "CCONJ"
                                                      or token.like_num
                                                      or token.is_space
                                                      or token.is_punct)]
    return words


if __name__ == "__main__":
    #test_db_connection()
    file_name = "dataset/all_abstracts.txt"
    output_file = "dataset/all_abstracts_clean1.txt"
    #print(read_all_keywords())

    #export(output=file_name)
    clean_text(file_name, output_file)
    result = preprocess_text("onions")
    print(result)