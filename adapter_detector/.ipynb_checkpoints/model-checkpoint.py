import os

os.environ['KERAS_BACKEND'] = 'tensorflow'
from keras.models import load_model


def load_adapter_model(model_fn=None):
    if model_fn is None:
        path = os.path.split(__file__)[0]
        model_fn = os.path.join(path,
                                'data',
                                'adapter_model.h5')
    if not os.path.exists(model_fn):
        raise OSError('Model file {} does not exist'.format(model_fn))
    model = models.load_model(model_fn)
    return model
