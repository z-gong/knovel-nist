#!/usr/bin/env python3
# coding=utf-8

import os
import sys
import time
import json
import requests
import random


def gen_headers(cid):
    return {'Pragma'         : 'no-cache',
            'Accept-Encoding': 'gzip, deflate, br',
            'Accept-Language': 'en-US,en;q=0.9,zh-CN;q=0.8,zh;q=0.7',
            'User-Agent'     : 'Mozilla/5.0 (Windows NT 6.1; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.139 Safari/537.36',
            'Accept'         : '*/*',
            'Referer'        : 'https://app.knovel.com/web/poc/ms/profile.v?cid=%(cid)s' % ({'cid': cid}),
            'Cookie'         : 'AMCV_4D6368F454EC41940A4C98A6@AdobeOrg=-1330315163|MCIDTS|17672; SilentAuthAttemptCounter-prod={"count":1}',
            'Connection'     : 'keep-alive',
            'Cache-Control'  : 'no-cache',
            }


def get_session_id(headers):
    r = requests.post('https://app.knovel.com/api/',
                      json={"APIKEY"          : "EE109DAC-CDA6-4299-8095-F6D6DBF96380",
                            "APPLICATION_NAME": "web",
                            "CLIENT_TRACK_ID" : "knovel",
                            "METHOD"          : {
                                "NAME"  : "login/authenticate_silent",
                                "PARAMS": {
                                    "REMEMBER_ME_DISABLE": False
                                },
                                "TYPE"  : "POST"
                            },
                            "OPTIONAL_INPUT"  : {
                                "ANALYTICS_ID": "096efccc-bcac-4fce-a5b4-f015b3ca9624",
                                "DOMAIN_NAME" : "moc.levonk.ppa",
                                "CUID"        : "CTcc4160ef-f1d9-dad6-7163-7eab720e21ce-2075"
                            },
                            "SESSION"         : {}
                            },
                      headers=headers
                      )

    jout = r.json()
    session_id = jout.get('BODY').get('SESSION').get('ID')
    return session_id


def get_ekey(cid, session_id, headers):
    r = requests.post('https://app.knovel.com/api/',
                      json={"APIKEY"          : "EE109DAC-CDA6-4299-8095-F6D6DBF96380",
                            "APPLICATION_NAME": "web",
                            "CLIENT_TRACK_ID" : "knovel",
                            "METHOD"          : {
                                "NAME"  : "entitlements/get-license",
                                "PARAMS": {
                                    "CONTENT_ID"    : cid,
                                    "CONTENT_FORMAT": "dbkms",
                                    "BOOK"          : False,
                                    "META"          : False,
                                    "CHILD"         : False,
                                    "BROWSE_FROM"   : "",
                                    "ROW_CONTENT_ID": "",
                                    "ACTION"        : "read"
                                },
                                "TYPE"  : "GET"
                            },
                            "OPTIONAL_INPUT"  : {
                                "DOMAIN_NAME": "moc.levonk.ppa"
                            },
                            "SESSION"         : {
                                "ID": session_id,
                            }
                            },
                      headers=headers
                      )

    jout = r.json()
    # print(jout)
    ekey = jout.get('BODY').get('EKEY')
    return ekey


def get_constant(cid, ekey, session_id, headers):
    r = requests.post('https://app.knovel.com/api/',
                      json={"APIKEY"          : "EE109DAC-CDA6-4299-8095-F6D6DBF96380",
                            "APPLICATION_NAME": "web",
                            "CLIENT_TRACK_ID" : "knovel",
                            "METHOD"          : {
                                "NAME"  : "kms-purechem/compound-constant",
                                "PARAMS": {
                                    "CONTENT_ID"      : cid,
                                    "SOURCE_ID"       : "all",
                                    "DATA_QUALITY"    : "dqRECM",
                                    "EKEY"            : ekey,
                                    "UNIT_PERSIST_IDS": "{}",
                                    "OUTPUT_FORMAT"   : "JSON",
                                    "TEMPLATE"        : "PROFILE"
                                },
                                "TYPE"  : "GET"
                            },
                            "OPTIONAL_INPUT"  : {
                                "DOMAIN_NAME": "moc.levonk.ppa"
                            },
                            "SESSION"         : {
                                "ID": session_id,
                            }
                            },
                      headers=headers
                      )

    jout = r.json()
    return jout


def get_prop_data(cid, prop_id, phase_id, ekey, session_id, headers):
    r = requests.post('https://app.knovel.com/api/',
                      json={"APIKEY"          : "EE109DAC-CDA6-4299-8095-F6D6DBF96380",
                            "APPLICATION_NAME": "web",
                            "CLIENT_TRACK_ID" : "knovel",
                            "METHOD"          : {
                                "NAME"     : "kms-purechem/compound-thermodata-property-data",
                                "PARAMS"   : {
                                    "CONTENT_ID"       : cid,
                                    "SOURCE_ID"        : "all",
                                    "DATA_QUALITY"     : "dqRECM",
                                    "PROPERTY_ID"      : prop_id,
                                    "PROPERTY_PHASE_ID": phase_id,
                                    "UNIT_PERSIST_IDS" : "",
                                    "EKEY"             : ekey,
                                    "OUTPUT_FORMAT"    : "JSON",
                                    "APPLY_DOWNSAMPLE" : "TRUE",
                                    "DOWNSAMPLE_ON"    : "isobar",
                                    "PARAMETERS"       : "",
                                    "TEMPLATE"         : "PROFILE"
                                },
                                "TYPE"     : "GET",
                                "CLIENT_IP": "202.120.51.74"
                            },
                            "OPTIONAL_INPUT"  : {
                                "USER_AGENT"  : "Mozilla/5.0 (Windows NT 6.1; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.139 Safari/537.36",
                                "WEB_SERVER"  : "10.0.20.245",
                                "REFERRER"    : "https://app.knovel.com/web/poc/ms/profile.v?cid=%(cid)s" % (
                                    {'cid': cid}),
                                "ROUTE"       : "/",
                                "ANALYTICS_ID": "",
                                "DOMAIN_TOKEN": "",
                                "DOMAIN_NAME" : "moc.levonk.ppa"
                            },
                            "SESSION"         : {
                                "ID": session_id,
                            }
                            },
                      headers=headers
                      )

    jout = r.json()
    return jout


def print_prop_data(json_dict):
    data = json_dict.get('BODY').get('COMPOUND').get('DATA')
    n_points = data.get('DATA_POINTS_COUNT')
    if n_points > 0:
        first = data.get('DATA_POINTS')[0]
        print(n_points, first)
    else:
        print(n_points)


prop_list = [
    ['prVP', 'prVAPRLG'],  # vapor pressure
    # # ['prEN', 'prENTHLG'],  # enthalpy of saturated liquid
    # # ['prEN', 'prENTHIG'],  # enthalpy of ideal gas
    ['prDENS', 'prDENLG'],  # density of saturated liquid
    ['prDENS', 'prDENGL'],  # density of saturated gas
    # # ['prDENS', 'prDENL'],  # density of liquid
    ['prESU', 'prENSUCG'],  # enthalpy of sublimation and vaporization

    # # ['prHCACP', 'prHECACOIG'],  # heat capacity at constant pressure of ideal gas
    # # ['prHCACP', 'prHECACOL'],  # heat capacity at constant pressure of liquid
    ['prHCASP', 'prHCASPLG'],  # heat capacity at saturated pressure of liquid
    # # ['prCOIE', 'prCOISEXL'],  # coefficient of isobaric expansion of liquid
    # # ['prIC', 'prISCOL'],  # isothermal compressibility of liquid
    ['prSOS', 'prSPSOLG'],  # speed of sound of saturated liquid
    ['prDV', 'prDYVILG'],  # dynamic viscosity of saturated liquid
    # # ['prDV', 'prDYVIL'], # dynamic viscosity of liquid
    ['prSUTE', 'prSURTENLG'],  # surface tension
    ['prTC', 'prTHCOLG'],  # thermal conductivity of saturated liquid
    # # ['prTC', 'prTHCOL'], # thermal conductivity of liquid
]

DATA_DIR = '/home/gongzheng/knovel-nist-crawler-json/'

with open(sys.argv[1]) as f:
    lines = f.read().splitlines()

for line in lines:
    line = line.strip()
    if line.startswith('#') or line == '':
        continue
    words = line.split()

    no = int(words[0])
    cid = words[7]

    # if no <= 1861:
    #     continue

    filename = DATA_DIR + 'constant/%i-%s' % (no, cid)

    headers = gen_headers(cid)

    GET_KEY = False
    while not GET_KEY:
        session_id = get_session_id(headers)
        ekey = get_ekey(cid, session_id, headers)
        if ekey == '':
            print('Blocked! Retry after 600s...')
            time.sleep(600)
        else:
            GET_KEY = True

    j_constant = get_constant(cid, ekey, session_id, headers)
    with open(filename, 'w') as f:
        f.write(json.dumps(j_constant))

    for prop_id, phase_id in prop_list:
        filename = DATA_DIR + 'prop-data/%i-%s-%s-%s' % (no, cid, prop_id, phase_id)
        if os.path.exists(filename):
            continue

        data_dict = get_prop_data(cid, prop_id, phase_id, ekey, session_id, headers)
        print('%5i %10s %10s %10s ' % (no, cid, prop_id, phase_id), end='')
        time.sleep(2 + random.random() * 2)

        try:
            print_prop_data(data_dict)
        except Exception as e:
            print(repr(e))
        else:
            with open(filename, 'w') as f:
                f.write(json.dumps(data_dict))

    time.sleep(2 + random.random() * 5)
