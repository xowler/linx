import numpy

def Load(fn):
    f = numpy.load(fn)
    d =  f['arr_0']
    f.close()
    return d


def PCA(data,copy=False):
    import sklearn.decomposition
    if copy: assert 0
    m = numpy.mean(data,axis=0)
    data = data - m
    s = numpy.std(data,axis=0)
    data = data/s
    pca = sklearn.decomposition.PCA()
    data = pca.fit_transform(data)

    d = {
        "c":pca.components_,
        "ev":pca.explained_variance_,
        "evr":pca.explained_variance_ratio_,
        "m_":pca.mean_,
        "nc":pca.n_components,
        "m":m,
        "s":s,
    }
    return data, d


def PCAFromNPZ(fn):
    pca = sklearn.decomposition.PCA()
    d = numpy.load(fn)
    pca.components_=d["c"]
    pca.explained_variance_=d["ev"]
    pca.explained_variance_ratio_=d["evr"]
    pca.mean_=d["m_"]
    pca.n_components=d["nc"]
    pca.m=d["m"]
    pca.s=d["s"]
    pca.fit = None
    pca.fit_transform = None
    pca.tranform = pca.transform
    return pca


def PCATransformFromNPZ(fn,data):
    d = numpy.load(fn)
    data = data - d['m']
    data = data/d['s']
    return numpy.dot(data, d['c'].T)




