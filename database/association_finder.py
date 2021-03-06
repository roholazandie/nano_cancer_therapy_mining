from database.pubmed_search import PubMedSearch
from config.nano_cancer_mining_configuration import NanoCancerConfiguration
import time

from remote.entrezsearch import EntrezSearch


class AssociationFinder(object):


    def __init__(self):
        nano_cancer_config = NanoCancerConfiguration()
        self.cancer_file = nano_cancer_config.file_config.cancer_file
        self.nano_particle_file = nano_cancer_config.file_config.nano_particle_file
        self.biosensor_file = nano_cancer_config.file_config.biosensor_file
        self._pubmed_search = PubMedSearch(nano_cancer_config)
        self._entrez_search = EntrezSearch()


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



    def cancer_nano_particle_association_bio(self):
        '''
        In this method we try to find all the articles with the name of
        a specific cancer name AND and a nano particle name and
        finally return the pubmedid of all the articles for any two
        combinations of cancer x nano_particle
        :return:
        '''
        cancer_names = sorted(list(
            set([name.rstrip().lower() for name in open(self.cancer_file).readlines() if len(name.rstrip()) != 0])))
        nano_particle_names = sorted(list(set(
            [name.rstrip().lower() for name in open(self.nano_particle_file).readlines() if len(name.rstrip()) != 0])))


        cancer_nano_particle_association = dict()
        for i, cancer_name in enumerate(cancer_names):
            cancer_nano_particle_association[cancer_name] = dict()
            for j, nano_particle_name in enumerate(nano_particle_names):
                try:
                    cancer_nano_particle_association[cancer_name][nano_particle_name] = self._entrez_search.rule1_query(cancer_name, nano_particle_name)
                except:
                    print("sleeping...")
                    cancer_nano_particle_association[cancer_name][nano_particle_name] = []
                    time.sleep(10)
                print("(", i,", ", j, ")")


        with open("cancer_nano_particle_association.txt", 'w') as fw:
            fw.write(str(cancer_nano_particle_association))



    def cancer_biosensor_association(self):
        '''
        In this method we try to find all the pubmedid of all the articles with
        cancer name in cancer_file and biosensor in biosensor_file
        :param cancer_file:
        :param biosensor_file:
        :return:
        '''
        cancer_names = [name.rstrip().lower() for name in open(self.cancer_file).readlines() if len(name.rstrip()) != 0]
        biosensor_names = [name.rstrip().lower() for name in open(self.biosensor_file).readlines() if len(name.rstrip()) != 0]

        cancer_biosensor_association = dict()
        for i, cancer_name in enumerate(cancer_names):
            cancer_biosensor_association[cancer_name] = dict()
            for j, biosensor_name in enumerate(biosensor_names):
                try:
                    cancer_biosensor_association[cancer_name][biosensor_name] = self._entrez_search.rule1_query(cancer_name, biosensor_name)
                except:
                    print("sleeping...")
                    cancer_biosensor_association[cancer_name][biosensor_name] = []
                    time.sleep(10)

                print("(", i, ", ", j, ")")

        with open("cancer_biosensor_association.txt", 'w') as fw:
            fw.write(str(cancer_biosensor_association))

if __name__ == "__main__":
    association_finder = AssociationFinder()
    #association_finder.cancer_nano_particle_association()
    #association_finder.cancer_nano_particle_association_bio()
    association_finder.cancer_biosensor_association()



