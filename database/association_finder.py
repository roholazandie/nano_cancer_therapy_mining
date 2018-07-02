from database.pubmed_search import PubMedSearch


class AssociationFinder(object):


    def __init__(self):
        self._pubmed_search = PubMedSearch()


    def cancer_nano_particle_association(self, cancer_file, nano_particle_file):
        '''
        In this method we try to find all the articles with the name of
        a specific cancer name AND and a nano particle name and
        finally return the pubmedid of all the articles for any two
        combinations of cancer x nano_particle
        :return:
        '''
        cancer_names = [name.rstrip() for name in open(cancer_file).readlines()]
        nano_particle_names = [name.rstrip() for name in open(nano_particle_file).readlines()]

        cancer_nano_particle_association = dict()
        for cancer_name in cancer_names:
            for nano_particle_name in nano_particle_names:
                cancer_nano_particle_association[cancer_name][nano_particle_name] = self._pubmed_search.and_query(cancer_name, nano_particle_name)



if __name__ == "__main__":
    cancer_file = "/home/rohola/Codes/Python/nano_cancer_therapy_mining/dataset/keywords.txt"
    nano_particle_file = ""
    association_finder = AssociationFinder()
    association_finder.cancer_nano_particle_association(cancer_file= cancer_file, nano_particle_file= nano_particle_file)



