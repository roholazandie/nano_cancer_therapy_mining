from gensim.models import KeyedVectors
from gensim.models.word2vec import LineSentence, Word2Vec
from database.export_database_to_text_file import preprocess_text

#TODO needs preprocessing, needs to clean and stem
class WordEmbeddingModel():

    def __init__(self, model_file):
        #self.model = KeyedVectors.load(model_file)
        '''
        '''

    def create_word2vec_model(self, input_file, output_file):
        sentences = LineSentence(input_file)
        model = Word2Vec(sentences, size=100, window=5, min_count=5, workers=4)
        model.save(output_file)
        return model


    def most_similar_words(self, word):
        return self.model.wv.most_similar(word)


    def most_similar_multi_token_words(self, multi_word):
        self.model = KeyedVectors.load(model_file)
        super_imposed_vector = 0
        for word in multi_word.split():
            super_imposed_vector += self.model.wv[word]

        super_imposed_vector = super_imposed_vector/len(multi_word)
        return self.model.wv.most_similar(positive=[super_imposed_vector], topn=10)


if __name__ == "__main__":
    model_file = "embeddings/pubmedword2vec_clean.bin"
    word_embedding_model = WordEmbeddingModel(model_file)

    #create_word2vec_model(input_file="dataset/all_abstracts_clean.txt", output_file="embeddings/pubmedword2vec_clean.bin")


    keywords = ["Cancer",
                "Antitumor",
                "Diagnosis",
                "Detection",
                "Therapy",
                "Imaging",
                "Nanomaterial",
                "Nanoparticle",
                "Nanocomposite",
                "Nanosheet",
                "Nanorod",
                "Nanosphere"]

    drugs = ["Fullerene",
            "Dendrimer",
            "Nanoliposome",
            "Liposome",
            "Nanocapsule",
            "Micelle",
            "Nanocrystal"]

    # with open("keyword_synonyms.txt", 'w') as file_writer:
    #     for keyword in keywords:
    #         file_writer.write(keyword+"\t"+str(most_similar_words(model_file="embeddings/pubmedword2vec_clean.bin", word=keyword.lower()))+"\n\n")


    # with open("drugs_synonyms.txt", 'w') as file_writer:
    #     for keyword in drugs:
    #         file_writer.write(keyword+"\t"+str(most_similar_words(model_file="embeddings/pubmedword2vec_clean.bin", word=keyword.lower()))+"\n\n")

    #print(most_similar_words(model_file="embeddings/pubmedword2vec_clean.bin", word="metalic"))

    with open("dataset/output.txt") as file_reader:
        for line in file_reader:
            cleaned_words = preprocess_text(line)
            cleaned_multi_word = " ".join(cleaned_words)
            try:
                print(word_embedding_model.most_similar_multi_token_words(multi_word=cleaned_multi_word))
            except:
                print("word "+str(cleaned_multi_word)+ "is not in the vocabulary")


