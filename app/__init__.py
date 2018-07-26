from flask import Flask
from flask_sqlalchemy import SQLAlchemy

from config import Config

app = Flask(__name__)
app.config.from_object(Config)
app.jinja_env.auto_reload = True

db = SQLAlchemy(app)

from .main import main as main_blueprint

app.register_blueprint(main_blueprint)


@app.template_filter()
def r(f, n):
    ''' round float to n decimal'''
    if f == None:
        return None
    else:
        return round(f, n)

@app.template_filter()
def rint(f):
    ''' round float to int'''
    if f == None:
        return None
    else:
        return int(round(f))