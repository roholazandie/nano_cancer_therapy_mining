from database.pubmed_search import PubMedSearch
from config.nano_cancer_mining_configuration import NanoCancerConfiguration
import numpy as np

class AssociationFinder(object):


    def __init__(self):
        nano_cancer_config = NanoCancerConfiguration()
        self.cancer_file = nano_cancer_config.file_config.cancer_file
        self.nano_particle_file = nano_cancer_config.file_config.nano_particle_file
        self.biosensor_file = nano_cancer_config.file_config.biosensor_file
        self._pubmed_search = PubMedSearch(nano_cancer_config)


    def cancer_nano_particle_association(self):
        '''
        In this method we try to find all the articles with the name of
        a specific cancer name AND and a nano particle name and
        finally return the pubmedid of all the articles for any two
        combinations of cancer x nano_particle
        :return:
        '''
        cancer_names = [name.rstrip() for name in open(self.cancer_file).readlines()]
        nano_particle_names = [name.rstrip() for name in open(self.nano_particle_file).readlines()]


        cancer_nano_particle_association = dict()
        for i, cancer_name in enumerate(cancer_names):
            cancer_nano_particle_association[cancer_name] = dict()
            for j, nano_particle_name in enumerate(nano_particle_names):
                cancer_nano_particle_association[cancer_name][nano_particle_name] = self._pubmed_search.rule1_query(cancer_name.lower(), nano_particle_name.lower())
                if j%10==0:
                    with open("backup"+str(j)+".txt", 'w') as fw:
                        fw.write(str(cancer_nano_particle_association))

        with open("final.txt", 'w') as fw:
            fw.write(str(cancer_nano_particle_association))


    def cancer_biosensor_association(self):
        '''
        In this method we try to find all the pubmedid of all the articles with
        cancer name in cancer_file and biosensor in biosensor_file
        :param cancer_file:
        :param biosensor_file:
        :return:
        '''
        cancer_names = [name.rstrip() for name in open(self.cancer_file).readlines()]
        biosensor_names = [name.rstrip() for name in open(self.biosensor_file).readlines()]

        cancer_biosensor_association = dict()
        for cancer_name in cancer_names:
            for biosensor_name in biosensor_names:
                cancer_biosensor_association[cancer_name][biosensor_name] = self._pubmed_search.rule1_query(cancer_name, biosensor_name)



if __name__ == "__main__":
    association_finder = AssociationFinder()
    association_finder.cancer_nano_particle_association()



