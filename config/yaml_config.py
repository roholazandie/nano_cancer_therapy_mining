import yaml


class YamlConfiguration():

    def __init__(self):
        '''
        this class loads yaml configuration
        '''
        self._yaml_dict = None
        self.yaml_file = "../config.yaml"
        self.load_from_file()


    def load_from_file(self):
        with open(self.yaml_file, encoding="utf-8") as file_reader:
            self._yaml_dict = yaml.load(file_reader)


    @property
    def config(self):
        return self._yaml_dict


    def get_value(self, key):
        if key in self._yaml_dict:
            return self._yaml_dict[key]



if __name__ == "__main__":
    yaml_file = "../config.yaml"
    yaml_config_file = YamlConfiguration(yaml_file)
    print(yaml_config_file._yaml_dict)