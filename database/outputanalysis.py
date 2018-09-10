from config.nano_cancer_mining_configuration import NanoCancerConfiguration
import numpy as np
from visualization.plotlyvisualize import visualize_associations, bar_chart_plot, scatter_plot, scatter3d_plot
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import scale, normalize
from sklearn.cluster import KMeans
from scipy import linalg
from sklearn import manifold
import operator
import xmltodict
from remote.entrezsearch import EntrezSearch

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

        cancer_nanoparticle_dict = eval(open(self.dataset_dir+"cancer_nano_particle_association_abstract.txt").read())

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

        cancer_nanoparticle_dict = eval(open(self.dataset_dir+"cancer_nano_particle_association_abstract.txt").read())

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
        cancer_names, nano_particle_names, cancer_particle_associations = self.get_cancer_particle_assosiation()

        visualize_associations(nano_particle_names, cancer_names, cancer_particle_associations, output_file=self.output_dir+"heatmap_cancer_nano_particle")



    def svd_decomposition(self):
        cancer_names, nano_particle_names, cancer_particle_associations = self.get_cancer_particle_assosiation()

        # each cancer is a datapoint and paricles act like features, so normalization
        # apply in the first dimension
        # log normalization tansform
        method = "svd"

        cancer_particle_associations = np.log(cancer_particle_associations+1)
        #cancer_particle_associations = normalize(cancer_particle_associations)

        X_projected = self.svd_dimentionality_reduction(cancer_particle_associations, n_component=10)

        labels = self.clustering(X_projected, n_clusters = 5)

        scatter_plot(X_projected[:, 0], X_projected[:, 1], colors = labels, names=cancer_names, output_file= self.output_dir+ "scatter2d_" + method)
        scatter3d_plot(X_projected[:, 0], X_projected[:, 1], X_projected[:, 2], names=cancer_names, colors = labels, output_file=self.output_dir+ "scatter3d_" + method)

    def get_cancer_particle_assosiation(self):
        cancer_names = sorted(list(
            set([name.rstrip().lower() for name in open(self.cancer_file).readlines() if len(name.rstrip()) != 0])))
        nano_particle_names = sorted(list(set(
            [name.rstrip().lower() for name in open(self.nano_particle_file).readlines() if len(name.rstrip()) != 0])))
        cancer_nanoparticle_dict = eval(open(self.dataset_dir + "final.txt").read())
        cancer_particle_associations = np.zeros((len(cancer_names), len(nano_particle_names)))
        for i, cancer_name in enumerate(cancer_names):
            for j, nano_particle_name in enumerate(nano_particle_names):
                try:
                    if len(cancer_nanoparticle_dict[cancer_name][nano_particle_name]) > 0:
                        cancer_particle_associations[i][j] = len(
                            cancer_nanoparticle_dict[cancer_name][nano_particle_name])
                    else:
                        cancer_particle_associations[i][j] = 0
                except:
                    cancer_particle_associations[i][j] = 0
                    continue
        return cancer_names, nano_particle_names, cancer_particle_associations

    def manifold_dim_reduction(self):
        cancer_names, nano_particle_names, cancer_particle_associations = self.get_cancer_particle_assosiation()

        method = "se"

        cancer_particle_associations = np.log(cancer_particle_associations + 1)
        cancer_particle_associations = normalize(cancer_particle_associations)

        X_transformed = self.manifold_learning_transformations(cancer_particle_associations, method=method)


        # X_projected = self.svd_dimentionality_reduction(cancer_particle_associations, n_component=10)
        # labels = self.clustering(X_projected, n_clusters=5)

        labels = self.clustering(X_transformed, n_clusters=5)

        scatter_plot(X_transformed[:,0], X_transformed[:,1], names=cancer_names, colors=labels, output_file= self.output_dir+ "scatter2d_" + method)
        scatter3d_plot(X_transformed[:, 0], X_transformed[:, 1], X_transformed[:, 2],
                       names=cancer_names, colors=labels, output_file= self.output_dir+ "scatter3d_" + method)




    def clustering(self, X, n_clusters=5):
        kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(X)
        return kmeans.labels_



    def manifold_learning_transformations(self, X, method, n_components=3):
        # Perform Locally Linear Embedding Manifold learning
        n_neighbors = 10
        trans_data = {}

        if method == "Modified LLE":
            trans_data = manifold.LocallyLinearEmbedding(n_neighbors, n_components=n_components, method="modified").fit_transform(X)
        elif method == "LLE":
            trans_data = manifold.LocallyLinearEmbedding(n_neighbors, n_components=n_components, method="standard").fit_transform(X)
        elif method == "mds":
            trans_data = manifold.MDS(n_components=n_components).fit_transform(X)
        elif method == "se":
            trans_data = manifold.SpectralEmbedding(n_components=n_components, n_neighbors=n_neighbors).fit_transform(X)
        elif method == "tsne":
            trans_data["tsne"] = manifold.TSNE(n_components=n_components, init='pca', random_state=0).fit_transform(X)

        return trans_data



    def svd_dimentionality_reduction(self, X, n_component):
        svd = TruncatedSVD(n_components=n_component, n_iter=20, random_state=42)
        svd.fit(X)
        #print(svd.singular_values_)
        print(svd.explained_variance_ratio_)
        print(svd.explained_variance_ratio_.cumsum())
        X_projected = svd.fit_transform(X)
        #print(np.shape(X_projected))
        #print(X_projected[1,:])

        return X_projected



class MetaInformation():


    def __init__(self):
        pass


    def get_countries(self, raw_xml_string):
        info_dict = xmltodict.parse(raw_xml_string)
        country = info_dict['PubmedArticleSet']['PubmedArticle']['MedlineCitation']['MedlineJournalInfo']['Country']
        return country


    def get_journal_title(self, raw_xml_string):
        info_dict = xmltodict.parse(raw_xml_string)
        journal_title = info_dict['PubmedArticleSet']['PubmedArticle']['MedlineCitation']['Article']['Journal']['Title']
        return journal_title

    def list_of_authors(self, raw_xml_string):
        info_dict = xmltodict.parse(raw_xml_string)
        authors = info_dict['PubmedArticleSet']['PubmedArticle']['MedlineCitation']['Article']['AuthorList']['Author']
        authors = [(author["LastName"], author["ForeName"], author["Initials"]) for author in authors]
        return authors






if __name__ == "__main__":
    output_analysis = OutputAnalysis()
    #output_analysis.simplify_output()
    #output_analysis.visualize_association()
    #output_analysis.most_frequent_barchart_cancer_nano_particle(n_associations=100)
    #output_analysis.simplify_output()
    #output_analysis.visualize_association()
    #output_analysis.svd_decomposition()
    #output_analysis.manifold_dim_reduction()


    meta_info = MetaInformation()
    entrez_search = EntrezSearch()
    raw_xml_string = entrez_search.fetch(24366930)
    result = meta_info.list_of_authors(raw_xml_string)
    print(result)
