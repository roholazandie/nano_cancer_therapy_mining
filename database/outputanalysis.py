from config.nano_cancer_mining_configuration import NanoCancerConfiguration
import numpy as np
from visualization.plotlyvisualize import visualize_associations

class OutputAnalysis():

    def __init__(self):
        nano_cancer_config = NanoCancerConfiguration()
        self.cancer_file = nano_cancer_config.file_config.cancer_file
        self.nano_particle_file = nano_cancer_config.file_config.nano_particle_file
        self.biosensor_file = nano_cancer_config.file_config.biosensor_file


    def simplify_output(self):
        cancer_names = [name.rstrip() for name in open(self.cancer_file).readlines()]
        nano_particle_names = [name.rstrip() for name in open(self.nano_particle_file).readlines()]

        cancer_nanoparticle_dict = eval(open("/home/rohola/Codes/Python/nano_cancer_therapy_mining/dataset/cancer_nanoparticle_association.txt").read())

        with open("/home/rohola/Codes/Python/nano_cancer_therapy_mining/dataset/simple_cancer_nanoparticle_association.txt", 'w') as file_writer:
            for cancer_name in cancer_names:
                for nano_particle_name in nano_particle_names:
                    try:
                        #print(type(cancer_nanoparticle_dict[cancer_name][nano_particle_name]))
                        if len(cancer_nanoparticle_dict[cancer_name][nano_particle_name]) > 0:
                            file_writer.write(cancer_name+"\t"+nano_particle_name+"\t"+str(cancer_nanoparticle_dict[cancer_name][nano_particle_name])+"\n")
                    except:
                        continue



    def visualize_association(self):
        cancer_names = list(set([name.rstrip().lower() for name in open(self.cancer_file).readlines()]))
        nano_particle_names = list(set([name.rstrip().lower() for name in open(self.nano_particle_file).readlines()]))

        cancer_nanoparticle_dict = eval(open(
            "/home/rohola/Codes/Python/nano_cancer_therapy_mining/dataset/cancer_nanoparticle_association.txt").read())

        #cancer_names = cancer_names[0:50]
        #nano_particle_names = nano_particle_names[0:50]

        num_associations = np.zeros((len(cancer_names), len(nano_particle_names)))
        #num_associations = np.zeros((10, 10))
        for i, cancer_name in enumerate(cancer_names):
            for j, nano_particle_name in enumerate(nano_particle_names):
                try:
                    if len(cancer_nanoparticle_dict[cancer_name][nano_particle_name]) > 0:
                        num_associations[i][j] = len(cancer_nanoparticle_dict[cancer_name][nano_particle_name])
                    else:
                        num_associations[i][j] = 0
                except:
                    num_associations[i][j] = 0
                    continue


        #print(np.shape(num_associations))
        #cancer_names = [cancer_name.lower()[0] for cancer_name in cancer_names]
        #nano_particle_names = [nano_particle_name.lower()[0] for nano_particle_name in nano_particle_names]


        visualize_associations(cancer_names, nano_particle_names, num_associations)






if __name__ == "__main__":
    output_analysis = OutputAnalysis()
    #output_analysis.simplify_output()
    output_analysis.visualize_association()