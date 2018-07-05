from config.yaml_config import YamlConfiguration



class DatabaseConfigs():

    def __init__(self, root, name, collection_name):
        self._root = root
        self._name = name
        self._collection_name = collection_name


    @property
    def root(self):
        return self._root

    @property
    def name(self):
        return self._name

    @property
    def collection_name(self):
        return self._collection_name




class FileConfigs():

    def __init__(self, cancer_file_name, nano_particle_file, biosensor_file_name, raw_xml_dir):
        self._cancer_file = cancer_file_name
        self._nano_particle_file = nano_particle_file
        self._biosensor_file = biosensor_file_name
        self._raw_xml_dir = raw_xml_dir


    @property
    def cancer_file(self):
        return self._cancer_file


    @property
    def biosensor_file(self):
        return self._biosensor_file


    @property
    def nano_particle_file(self):
        return self._nano_particle_file

    @property
    def raw_xml_dir(self):
        return self._raw_xml_dir


class NanoCancerConfiguration():

    def __init__(self, ):
        self._yaml_config_file = YamlConfiguration()
        cancer_file = self._yaml_config_file.config["files"]["cancer_names_file"]
        nano_particle_file = self._yaml_config_file.config["files"]["nano_particle_file"]
        biosensor_file = self._yaml_config_file.config["files"]["biosensor_file"]
        raw_xml_dir = self._yaml_config_file.config["files"]["raw_xml_dir"]
        self._file_configs = FileConfigs(cancer_file_name=cancer_file,
                                         nano_particle_file=nano_particle_file,
                                         biosensor_file_name=biosensor_file,
                                         raw_xml_dir=raw_xml_dir)

        root = self._yaml_config_file.config["database"]["root"]
        name = self._yaml_config_file.config["database"]["name"]
        collection_name = self._yaml_config_file.config["database"]["collection_name"]
        self._database_configs = DatabaseConfigs(root=root, name=name, collection_name=collection_name)


    @property
    def file_config(self):
        return self._file_configs

    @property
    def database_config(self):
        return self._database_configs

