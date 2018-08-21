from config.nano_cancer_mining_configuration import NanoCancerConfiguration
import numpy as np
from visualization.plotlyvisualize import visualize_associations, bar_chart_plot
import operator


class OutputAnalysis():

    def __init__(self):
        nano_cancer_config = NanoCancerConfiguration()
        self.cancer_file = nano_cancer_config.file_config.cancer_file
        self.nano_particle_file = nano_cancer_config.file_config.nano_particle_file
        self.biosensor_file = nano_cancer_config.file_config.biosensor_file
        self.output_dir = nano_cancer_config.file_config.output_dir
        self.dataset_dir = nano_cancer_config.file_config.dataset_dir


    def simplify_output(self):
        cancer_names = sorted(list(
            set([name.rstrip().lower() for name in open(self.cancer_file).readlines() if len(name.rstrip()) != 0])))
        nano_particle_names = sorted(list(set(
            [name.rstrip().lower() for name in open(self.nano_particle_file).readlines() if len(name.rstrip()) != 0])))

        cancer_nanoparticle_dict = eval(open(self.dataset_dir+"final.txt").read())

        frequency_dictionary = dict()
        for cancer_name in cancer_names:
            for nano_particle_name in nano_particle_names:
                try:
                    frequency_dictionary[cancer_name + " " + nano_particle_name] = len(
                        cancer_nanoparticle_dict[cancer_name][nano_particle_name])
                except:
                    frequency_dictionary[cancer_name + " " + nano_particle_name] = 0

        sorted_assiciation_frequency = sorted(frequency_dictionary.items(), key=operator.itemgetter(1), reverse=True)

        with open(self.output_dir+"simple_cancer_nanoparticle_association.txt", 'w') as file_writer:
            for assiciation, value in sorted_assiciation_frequency:
                file_writer.write(assiciation+"\t"+str(value)+"\n")



    def most_frequent_barchart_cancer_nano_particle(self, n_associations):
        cancer_names = sorted(list(
            set([name.rstrip().lower() for name in open(self.cancer_file).readlines() if len(name.rstrip()) != 0])))
        nano_particle_names = sorted(list(set(
            [name.rstrip().lower() for name in open(self.nano_particle_file).readlines() if len(name.rstrip()) != 0])))

        cancer_nanoparticle_dict = eval(open(self.dataset_dir+"final.txt").read())

        frequency_dictionary = dict()
        for cancer_name in cancer_names:
            for nano_particle_name in nano_particle_names:
                try:
                    frequency_dictionary[cancer_name+" "+nano_particle_name] = len(cancer_nanoparticle_dict[cancer_name][nano_particle_name])
                except:
                    frequency_dictionary[cancer_name + " " + nano_particle_name] = 0

        sorted_assiciation_frequency = sorted(frequency_dictionary.items(), key=operator.itemgetter(1), reverse=True)

        association_names = [x[0] for x in sorted_assiciation_frequency]
        values = [x[1] for x in sorted_assiciation_frequency]

        association_names = association_names[:n_associations]
        values = values[:n_associations]

        bar_chart_plot(x=association_names, y=values, output_file=self.output_dir+"barchart_cancer_nano_particle"+str(n_associations))


    def visualize_association(self):
        cancer_names = sorted(list(set([name.rstrip().lower() for name in open(self.cancer_file).readlines() if len(name.rstrip()) != 0])))
        nano_particle_names = sorted(list(set([name.rstrip().lower() for name in open(self.nano_particle_file).readlines() if len(name.rstrip()) != 0])))

        cancer_nanoparticle_dict = eval(open(
            "/home/rohola/Codes/Python/nano_cancer_therapy_mining/dataset/final.txt").read())

        num_associations = np.zeros((len(cancer_names), len(nano_particle_names)))
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

        visualize_associations(nano_particle_names, cancer_names, num_associations, output_file=self.output_dir+"heatmap_cancer_nano_particle")



if __name__ == "__main__":
    output_analysis = OutputAnalysis()
    #output_analysis.simplify_output()
    #output_analysis.visualize_association()
    #output_analysis.most_frequent_barchart_cancer_nano_particle(n_associations=1000)
    #output_analysis.simplify_output()
    output_analysis.visualize_association()