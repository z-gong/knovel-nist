#!/usr/bin/env python3
# coding=utf-8

import sys

sys.path.append('../ms-tools')

from app import app

app.run(host='0.0.0.0', port=4000, debug=True, threaded=True)