import math
import numpy as np
import pandas as pd
import re

from SauerClass import Record

# From Sauer Book Chapter
# dict_of_mass_abun = { 'C_m0':0.9893, 'C_m1':0.0107, 'H_m0':0.999885, \
#                      'H_m1':0.000115, 'N_m0':0.99632, 'N_m1':0.00368, \
#                      'O_m0':0.99757, 'O_m1':0.00038, 'O_m2':0.00205, \
#                      'Si_m0':0.922297, 'Si_m1':0.046832, 'Si_m2':0.030872 }
# Modified for Joachim
dict_of_mass_abun = { 'C_m0':0.98918, 'C_m1':0.01082, 'H_m0':0.999844, \
                      'H_m1':0.000156, 'N_m0':0.99634, 'N_m1':0.00366, \
                      'O_m0':0.99758, 'O_m1':0.00038, 'O_m2':0.00204, \
                      'Si_m0':0.922297, 'Si_m1':0.046832, 'Si_m2':0.030872 }

dict_of_mass_isotopes = {'C':[12,13], 'H':[1,2], 'N':[14,15], 'O':[16,17,18],\
                             'Si':[28,29,30]}
# From Sauer Book Chapter
# isotopic_mass_vectors = {'C':[0.9893, 0.0107], 'H':[0.999885, 0.000115], \
#                         'N':[0.99632, 0.00368], 'O':[0.99757, 0.00038, \
#                         0.00205], 'Si':[0.922297, 0.046832, 0.030872]}

# Modified for Joachim
isotopic_mass_vectors = {'C':[0.98918, 0.01082], 'H':[0.999844, 0.000116], \
                         'N':[0.99634, 0.00366], 'O':[0.99758, 0.00038, \
                         0.00204], 'Si':[0.922297, 0.046832, 0.030872]}

def full_correction_matrix(species_dict, N, cauchy=True):
    """Implementation of Annika et. al., equation 3. Wrapper function to: 
    1. call correction_matrix_cauchy() or correction_matrix_species() depending on user input.
    Cauchy by default. 
    2. Multiplies the C, O, N, H, and Si matricies (left to right), in that order. 

    PARAMS
    ------
    species_dict: dictionary of 
    N: ??unknown??
    cauchy: bool; specifies whether or not to use correction_matrix_cauchy() (True) or correction_matrix_species() (False)

    RETURNS
    -------
    result: matrix product of all species matricies. 
    """
    #print("Caution: using modified natural abundance figures")
    species_list = ['C', 'O', 'N', 'H', 'Si']
    matrix_list = []

    for species in species_list:
        if species_dict[species] > 0:
            if cauchy:
                matrix = correction_matrix_cauchy(N, species, species_dict[species])
            else:
                matrix = correction_matrix_species(species, species_dict[species], N)
            matrix_list.append(np.matrix(matrix))

    # Repeated matrix multiplication
    result = matrix_list[0]
    for matrix in matrix_list[1:]:
        result = result*matrix

    return result


def correction_matrix_cauchy(N, species, n):
    # It can be shown (Milica Thesis), that the
    # terms along the first column of the correction
    # are equivalent to the cauchy product of the
    # isotopic mass vectors

    initial = isotopic_mass_vectors[species]
    final = initial

    # Essentially the first convolution has
    # been done above by setting initial
    # so only do n-1 more 
    for i in range(n-1): 
        final = np.convolve(final, initial)

    matrix = np.zeros((N,N))

    for i in range(N):
        for j in range(N):
            if i < j:
                matrix[i,j] = 0.0
            elif j == i:
                matrix[i,j] = final[0]
            elif j < i:
                try:
                    #print i-j, final[i-j]
                    matrix[i,j] = final[(i-j)]
                except:
                    matrix[i,j] = 0.0

    #print(species,  matrix)
    return matrix


def correction_matrix_species(species, n,  N):

    matrix = np.zeros((N,N))

    for i in range(N):
        for j in range(N):

            if i == j:
                matrix[i,j] = correction_matrix_element(n,0,0,species)
                
            elif j > i:
                matrix[i,j] = 0.0
            else:
                # only m0 and m1 for C, H, N
                if species in ['C', 'H', 'N']:
                    if (n-(i-j)) >=0:
                        matrix[i,j] = correction_matrix_element((n-(i-j)), \
                                                   (i-j), 0, species)
                    else:
                        matrix[i,j] == 0.0
                # Si and O have m0, m1, m2                        
                elif species in ['Si', 'O']:
                    if i-j==1:
                        matrix[i,j] = correction_matrix_element(n-1, 1, 0, species)
                    elif i-j==2:
                        matrix[i,j] = correction_matrix_element(0, n,0, species) + correction_matrix_element(n-1, 0, 1, species)
                    elif abs(i-j) == 3:
                        matrix[i,j] = correction_matrix_element(0,n-1,1,species)
                    elif abs(i-j) >3:
                        matrix[i,j] = 0.0

    return matrix


def correction_matrix_element(m0, m1, m2, S):
    """
    m0 is the number of non-isotope atoms of this species in the fragment
    m1 is the number of isotope atoms with m1 in the fragment
    m2 is the number of isotope atoms with m2 in the fragment
    S is the species "C", "O", "H", "Si", "S"
    """
    sum_term = np.math.factorial(m0 + m1 + m2)
    
    if S in ['H', 'C', 'N', 'Si', 'O']:  # These have only m0 and m1
        
        if m0 > 0:
            numer_1 = math.pow(dict_of_mass_abun[S+'_m0'], float(m0))
            denom_1 = np.math.factorial(m0)
        else:
            numer_1 = 1
            denom_1 = 1

        if m1 > 0:
            numer_2 = math.pow(dict_of_mass_abun[S+'_m1'], float(m1))
            denom_2 = np.math.factorial(m1)
        else:
            numer_2 = 1
            denom_2 = 1

        if S in ['O', 'Si'] and m2 > 0:
            numer_3 = math.pow(dict_of_mass_abun[S+'_m2'], float(m2))
            denom_3 = np.math.factorial(m2)
        else:
            numer_3 = 1
            denom_3 = 1
                             

        term1 = float(numer_1)/float(denom_1)
        term2 = float(numer_2)/float(denom_2)
        term3 = float(numer_3)/float(denom_3)

    
        return sum_term*term1*term2*term3

    else:
        print("error")
        return 0


def correct_unlabelled(mdva, mdvun, f_unlabelled):
    """Implementation of Nanchen et. al., equation (5). 
    """
    numer = mdva - f_unlabelled*mdvun
    denom = 1 - f_unlabelled

    return numer/denom


def fractional_labelling(mdvaa):
    """
    """
    numer = 0.0
    denom = 0.0

    for i, elem in enumerate(mdvaa):
        numer = numer + i*elem
        denom = denom + (len(mdvaa)-1)*elem

    return numer/denom


def read_mhunter_csv(infile, verbose=False):
    """
    Reads an input file with isotopic labelling information coming
    directly from Agilent Mass Hunter software

    PARAMS
    ------
    infile: str; The name of the input file

    RETURNS
    -------
    record_list: A list of SauerClass.Record containing isotopic labelling information
    """
    record_list = []
    isotope_list = []
    with open(infile, "r") as fp:
        lines = fp.readlines()
    
    # use regex to parse out the species name
    pattern = r"M\+?\d{1,2} Results"
    for i, word in enumerate(lines[0].replace("\n", "").split(',')):
        m_pattern = re.findall(pattern, word)
        if len(m_pattern) > 0:
            species_name = word.replace(" " + m_pattern[0], "")
            
            test = False
            for record in record_list:
                if record.name == species_name:
                    test = True
                    record.add_record_position(i)
                else:
                    pass
            if test == False: # if we haven't seen this before
                record = Record(species_name)
                record.add_record_position(i)
                record_list.append(record)
                
    name_position = -1

    for i, word in enumerate(lines[1].split(',')):
        if word == "Name":
            name_position = i

    for line in lines[2:]:
        words = line.split(',')
        #print words
        if name_position > -1:
            sample_name = words[name_position]
        else:
            sample_name = "NotFound"

        for record in record_list:
            positions = record.get_record_positions()
            isotope_list = []
            for position in positions:
                #if verbose:
                    #print(record.name, " : ", position, " : ", words[position])
                if words[position] != "":
                    try:
                        isotope_list.append(float(words[position]))
                    except(ValueError):
                        print("Could not convert", words[position])
                else:
                    isotope_list.append(0.0)

            record.add_label_isotope_list(sample_name, isotope_list)

    return record_list


def read_hunter_contents(mh_contents):
    """Read a list of lists coming directly from Agilent Mass Hunter software into a
    SauerClass.Record object. Adapted from read_mhunter_csv() for integration into Flask webapp.

    PARAMS
    ------
    mh_contents: list of lists; MH csv input. 

    RETURNS
    -------
    record: a SauerClass.Record
    """
    record_list = []
    isotope_list = []

    # use regex to parse out the species name
    pattern = r"M\+?\d{1,2} Results"
    for i, word in enumerate(mh_contents[0]):
        m_pattern = re.findall(pattern, word)
        if len(m_pattern) > 0:
            species_name = word.replace(" " + m_pattern[0], "")
            
            test = False
            for record in record_list:
                if record.name == species_name:
                    test = True
                    record.add_record_position(i)
                else:
                    pass
            if test == False: # if we haven't seen this before
                record = Record(species_name)
                record.add_record_position(i)
                record_list.append(record)
                
    name_position = -1

    for i, word in enumerate(mh_contents[1]):
        if word == "Name":
            name_position = i

    for line in mh_contents[2:]:
        words = line
        #print words
        if name_position > -1:
            sample_name = words[name_position]
        else:
            sample_name = "NotFound"

        for record in record_list:
            positions = record.get_record_positions()
            isotope_list = []
            for position in positions:
                #if verbose:
                    #print(record.name, " : ", position, " : ", words[position])
                if words[position] != "":
                    try:
                        isotope_list.append(float(words[position]))
                    except(ValueError):
                        print("Could not convert", words[position])
                else:
                    isotope_list.append(0.0)

            record.add_label_isotope_list(sample_name, isotope_list)

    return record_list


def read_atomic_composition(infile):
    """
    @summary: reads the numbers of each atomic species, and the number
    of labelled carbon atoms for each molecular species

    @param infile: the name of the input file
    @type infile: strType

    @return atomic_composition: A dictionary with the atomic composition
                                of each species
    @type atomic_composition: dictType

    @return N_dict: A dictionary of the number of labelled carbons for each
                    species
    @type N_dict: dictType
    """
    with open(infile, "r") as fp:
        lines = fp.readlines()

    atomic_composition = {}
    N_dict = {}

    for line in lines:
        atomic_dict = {}
        words = line.split(',')

        N_dict[words[0]] = int(words[2])

        for word in words[1].split(' '):
            parts = re.split('(\d+)', word)
            atomic_dict[parts[0].lstrip().strip('\"')] = int(parts[1])
        atomic_composition[words[0].strip('\"')] = atomic_dict

    # get rid of all "s and 's
    for key in N_dict.keys():
        name = re.sub(r'\W+', '', key)
        
    for key in atomic_composition.keys():
        name = re.sub(r'\W+', '', key)
        
    return atomic_composition, N_dict



def read_atomic_contents(contents):
    """
    Reads the numbers of each atomic species, and the number
    of labelled carbon atoms for each molecular species. Adapted from read_atomic_contents() for 
    integration into Flask Webapp.

    @param infile: the name of the input file
    @type infile: strType

    @return atomic_composition: A dictionary with the atomic composition
                                of each species
    @type atomic_composition: dictType

    @return N_dict: A dictionary of the number of labelled carbons for each
                    species
    @type N_dict: dictType
    """
    atomic_composition = {}
    N_dict = {}

    for line in contents:
        atomic_dict = {}
        words = line

        N_dict[words[0]] = int(words[2])

        for word in words[1].split(' '):
            parts = re.split('(\d+)', word)
            atomic_dict[parts[0].lstrip().strip('\"')] = int(parts[1])
        atomic_composition[words[0].strip('\"')] = atomic_dict

    # get rid of all "s and 's
    for key in N_dict.keys():
        name = re.sub(r'\W+', '', key)
        
    for key in atomic_composition.keys():
        name = re.sub(r'\W+', '', key)
        
    return atomic_composition, N_dict
            

def calculate_labelling(record, N_dict, atomic_dict, verbose=False):
    """
    Calculates percentage labelling based on the method of Sauer et al.

    PARAMS
    ------
    record: a class containing a species name, and m0 to mn values of isotopic contribution. 
    Class.Record
    N_dict:A dictionary of the number of labelled carbons for each species (which is the rank of the mdva matrix), where each key is a sample. 
    atomic_dict: A dictionary with the atomic composition of each species
    verbose: bool; verbosity parameter

    RETURNS
    -------
    results_dict: a dictionary of sample name: labelling % pairs
    """
    
    name = record.get_name()
    species_dict = atomic_dict[name]
    #print species_dict

    N = N_dict[name] + 1
    if verbose:
        print(name, N)
    
    convolved_matrix = full_correction_matrix(species_dict, N, cauchy=True)
    if verbose:
        print(convolved_matrix)
    inverse = convolved_matrix.I
    #print inverse
    results_dict = {}

    #now find ms_matrix
    correction_dict = record.get_background_list()
    for key, value in correction_dict.items():
        #print value
        #for val in value:
            #print float(val)/float(sum(value))
        try:
            ms_matrix = np.matrix(([[float(num)/float(sum(value[0:N]))] \
                                        for num in value[0:N]]))
        except(ZeroDivisionError):
            ms_matrix = np.matrix(([[0.0] for i in range(0,N)]))
        #print 'ms', ms_matrix
        #print 'inverse', inverse
        try:
            mdva = inverse*ms_matrix/(sum(inverse*ms_matrix))[0,0]
            labelling = fractional_labelling(mdva)
        except(ValueError):
            print("Error in record", name)
            print("<br />")
            print("size inverse: " ,inverse.shape[0], inverse.shape[1])
            print("size ms_matrix: ",ms_matrix.shape[0], ms_matrix.shape[1])
            print("<br /> ")
        #print "background\n"
        #print  mdva
        #print "\n"

        # the labeling_mdva_list has the fractional labeling as
        # it's first element, then the labeling into each m channel
        # m0, m1, m2 etc, from the mdva matrix
        
        label_mdva_list = []
        label_mdva_list.append(labelling[0,0])
        for number in mdva:
            label_mdva_list.append(float(number[0,0]))
        results_dict[key + "*"] = label_mdva_list
    
    working_dict = record.get_label_isotope_dict()
    
    #temporary counter
    i = 0

    for key, value in working_dict.items():
        #print value
        #for val in value:
            #print float(val)/float(sum(value))
        try:
            ms_matrix = np.matrix(([[float(num)/float(sum(value[0:N]))] for num in value[0:N]]))
        except(ZeroDivisionError):
            ms_matrix = np.matrix(([[0.0] for i in range(0,N)]))
        #except(TypeError):
        #    print i, num, key, value
        #print ms_matrix
        mdva = inverse*ms_matrix/(sum(inverse*ms_matrix))[0,0]
        #print "main\n"
        #print mdva
        #print "\n"
        labelling = fractional_labelling(mdva)
        #print labelling[0,0]

        label_mdva_list = []
        label_mdva_list.append(labelling[0,0])
        for number in mdva:
            label_mdva_list.append(float(number[0,0]))
        results_dict[key] = label_mdva_list

        i = i+1
        
    return results_dict


if __name__ == "__main__":
    species_dict = {'C':8, 'O':2, 'N':1, 'H':26, 'Si':2}
    N = 4
    
    my_matrix = full_correction_matrix(species_dict, N, cauchy=True)
    
    print(my_matrix)

    # this is the m0, m1 etc values of the fragment as measured
    ms_matrix = np.matrix(([0.6228],[0.1517], [0.0749], [0.1507]))

    inverse = my_matrix.I
    mdva =  inverse*ms_matrix/sum(inverse*ms_matrix)[0,0]
    
    # Accounts for unlabelled portion
    mdvun = np.matrix(([0.9682],[0.0314],[0.0003],[0.0]))
    f_unlabelled = 0.01

    mdvaa = correct_unlabelled(mdva, mdvun, f_unlabelled)

    labelling = fractional_labelling(mdvaa)
    labelling_unc = fractional_labelling(mdva)

    print('mdva:')
    print(mdva)
    print('')
    print('mdvaa')
    print(mdvaa)
    print('')
    print('Labelling: %1.2f' %labelling[0,0])
    print('Labelling uncorrected: %1.2f' %labelling_unc[0,0])
