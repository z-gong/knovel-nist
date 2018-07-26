import os


class Config:
    CWD = os.path.dirname(os.path.abspath(__file__))
    DB_NAME = 'nist.sqlite'
    DB_PATH = os.path.join(CWD, 'database', DB_NAME)

    SQLALCHEMY_DATABASE_URI = 'sqlite:///%s' % DB_PATH
    SQLALCHEMY_BINDS = {
        'nist': 'sqlite:///' + os.path.join(CWD, 'database', 'nist.sqlite')
    }
    SQLALCHEMY_COMMIT_ON_TEARDOWN = True
    SQLALCHEMY_TRACK_MODIFICATIONS = True