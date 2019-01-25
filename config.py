import os


class Config:
    CWD = os.path.dirname(os.path.abspath(__file__))
    DB_PATH = ':memory:'

    SQLALCHEMY_DATABASE_URI = 'sqlite:///%s' % DB_PATH
    SQLALCHEMY_BINDS = {
        'nist'   : 'sqlite:///' + os.path.join(CWD, 'database', 'nist.sqlite'),
        'pubchem': 'sqlite:///' + os.path.join(CWD, 'database', 'pubchem.sqlite'),
    }
    SQLALCHEMY_COMMIT_ON_TEARDOWN = True
    SQLALCHEMY_TRACK_MODIFICATIONS = True
