import cPickle as pickle

def load_pickle(fname):
    with open(fname, 'rb') as input:
        return pickle.load(input)

def save_pickle(obj, fname='pickle_file.dat'):
    with open(fname, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
