class Record(object):
    """
    """

    def __init__(self, name):
        """
        @summary: Initialise the Record

        @param name: The name of the molecular species (e.g. GABA etc)
        @type name: strType
        """
        self.name = name
        if '\"' in name:
            name = name.strip('\"')
        elif "\'" in name:
            name = name.strip("\'")
        
        # init attributes (dictionaries)
        # following is used to store actual data
        self.label_isotope_dict = {}
        self.background_dict = {}
        self.position_list = []
        #self.label_isotope_dict[label] = isotope_list


    def add_label_isotope_list(self, label,isotope_list):
        """
        @summary: add a record consisting of a sample name and list of
                  intensity for each isotope
 
        @param label: The sample name (should be unique)
        @type label: strType

        @param isotope_list: a list of intensity of each isotope of the
                             species recorded, starting from m0
        @type isotope_list: listType
        """

        self.label_isotope_dict[label] = isotope_list


    def add_record_position(self, position):
        """ tells the read_mhunter_csv script which
            postions the species was found in"""

        self.position_list.append(position)


    def get_record_positions(self):
        return self.position_list


    def get_label_isotope_dict(self):
        return self.label_isotope_dict


    def add_background_list(self, label, isotope_list):
        """
        @summary: add a record consisting of a sample name and list of
                  intensity for each isotope
 
        @param label: The sample name (should be unique)
        @type label: strType

        @param isotope_list: a list of intensity of each isotope of the
                             species recorded, starting from m0
        @type isotope_list: listType
        """
        self.background_dict[label] = isotope_list    


    def get_background_list(self):
        return self.background_dict


    def set_rt(self, rt):
        """
        @summary: set the retention time
        
        @param rt: the retention time where the molecular species is found
        @type rt: floatType
        """
        self.rt = rt


    def get_rt(self, rt):
        try:
            return self.rt()
        except:
            print("Error: rt not set")
            return 0


    def get_name(self):
        return self.name


    def background_correct(self):
        """
        @summary: does background correction
        """
        full_bg = []
        
        if len(self.background_dict) > 0:
            # first figure out the length of the vector
            max = 0
            for record in self.background_dict.values():
                if len(record) > max:
                    max = len(record)
            # now prepare a zeroed matrix
            for i in range(max):
                full_bg.append(0)
            # do vector addition
            for list_of_bg in self.background_dict.values():
                for i,num in enumerate(list_of_bg):
                    full_bg[i] = full_bg[i] + num
            # finally, find the average
            for value in full_bg:
                value = value/len(self.background_dict)

            # now subtract background from all records
            for is_list in self.label_isotope_dict.values():
                for i,value in enumerate(is_list):
                    value = value - full_bg[i]

            print("background correction done")
        else:
            print("no background results recorded")


class Labelling(object):
    def __init__(self, species, label_dict):
        self.species = species
        self.label_dict = label_dict

    def get_species(self):
        return self.species

    def get_label_dict(self):
        return self.label_dict
