#!/usr/bin/env python3

import fire
import os
import sys
import shutil
import requests
import tarfile

sys.path.append('..')
from app.models import *


def discovery_results(txt):
    '''Get molecules using discovery-results API'''
    r = requests.post('https://app.knovel.com/api/',
                      json={
                          "APIKEY"          : "EE109DAC-CDA6-4299-8095-F6D6DBF96380",
                          "APPLICATION_NAME": "web",
                          "CLIENT_TRACK_ID" : "knovel",
                          "METHOD"          : {
                              "NAME"     : "kms-purechem/discovery-results",
                              "PARAMS"   : {
                                  "OPTIONS": "{\"pagination\":{\"limit\":5,\"offset\":0},\"queries\":{\"q\":[]},\"complex_property_filters\":{\"property\":[]},\"constant_property_filters\":{\"property\":[]},\"taxonomy_filters\":{\"taxonomy\":[]},\"unit_family\":{\"family\":\"imperial\"},\"columns\":{\"sort\":{\"property\":\"weight\",\"dir\":\"desc\"},\"column\":[\"prMATII\",\"prCT\",\"prCP\",\"prCD\"]}}"
                              },
                              "TYPE"     : "GET",
                              "CLIENT_IP": "202.120.51.74"
                          },
                          "OPTIONAL_INPUT"  : {
                              "USER_AGENT"  : "Mozilla/5.0 (Windows NT 6.1; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/64.0.3282.186 Safari/537.36",
                              "WEB_SERVER"  : "10.0.20.245",
                              "REFERRER"    : "http://devint-tester.knovel.net/",
                              "ROUTE"       : "/",
                              "ANALYTICS_ID": "",
                              "DOMAIN_TOKEN": "",
                              "DOMAIN_NAME" : "moc.levonk.ppa"
                          },
                          "SESSION"         : {
                              "ID": "5a9789c0-755e-7b25-f0a7-cac490e9eab9"
                          }
                      }
                      )
    j = r.json()

    data_list = j['BODY']['RESULTS']['data']
    print(len(data_list), data_list[0])

    text = ''
    for i, data in enumerate(data_list):
        cas = data['casrn']
        formula = data['hill_form']
        cid = data['content_id']
        name = data['name']
        kid = data['knovel_id']
        weight = data['prMATII']
        tc = data['prCT'] or ''
        pc = data['prCP'] or ''
        dc = data['prCD'] or ''
        line = '%i|%s|%s|%s|%s|%f|%s|%s|%s|%s\n' % (i + 1, cid, kid, cas, formula, weight, name, tc, pc, dc)
        text += line

    with open(txt, 'w') as f:
        f.write(text)


def _get_metadata(cid):
    '''Get SMILES using compound-meta API'''
    r = requests.post('https://app.knovel.com/api/',
                      json={
                          "APIKEY"          : "EE109DAC-CDA6-4299-8095-F6D6DBF96380",
                          "APPLICATION_NAME": "web",
                          "CLIENT_TRACK_ID" : "knovel",
                          "METHOD"          : {
                              "NAME"     : "kms-purechem/compound-metadata",
                              "PARAMS"   : {
                                  "CONTENT_ID": cid,
                                  "SOURCE_ID" : "srKN"
                              },
                              "TYPE"     : "GET",
                              "CLIENT_IP": "202.120.51.74"
                          },
                          "OPTIONAL_INPUT"  : {
                              "USER_AGENT"  : "Mozilla/5.0 (Windows NT 6.1; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/64.0.3282.186 Safari/537.36",
                              "WEB_SERVER"  : "10.0.20.245",
                              "REFERRER"    : "http://devint-tester.knovel.net/",
                              "ROUTE"       : "/",
                              "ANALYTICS_ID": "",
                              "DOMAIN_TOKEN": "",
                              "DOMAIN_NAME" : "moc.levonk.ppa"
                          },
                          "SESSION"         : {
                              "ID": "5a9789c0-755e-7b25-f0a7-cac490e9eab9"
                          }
                      }
                      )
    j = r.json()
    return j


def get_metadata(txt, txt_out):
    with open(txt) as f:
        lines = f.read().splitlines()

    f = open(txt_out, 'w')

    for line in lines:
        if line == '' or line.startswith('#'):
            continue
        words = line.split('|')
        no = words[0]
        cid = words[1]
        try:
            j = _get_metadata(cid)
            smiles = j['BODY']['COMPOUND']['METADATA']['SMILES']
            inchi = j['BODY']['COMPOUND']['METADATA']['InChi']
        except:
            print('\nERROR:', line)
        else:
            sys.stdout.write('\r\t%s %s %s %s' % (no, cid, smiles, inchi))
            sys.stdout.flush()

            f.write('%s|%s|%s\n' % (line, smiles, inchi))
            f.flush()

    f.close()


def create_db():
    db.create_all()


def insert_props():
    '''
    Insert properties
    Only saturated properties are inserted because they have the largest reliability
    '''

    prop_list = [
        ['prVP', 'prVAPRLG', 'pvap-lg'],  # vapor pressure
        ['prDENS', 'prDENLG', 'density-lg'],  # density of saturated liquid
        ['prDENS', 'prDENGL', 'density-gl'],  # density of saturated gas
        ##['prDENS', 'prDENL', 'density-l'], # density of liquid
        ['prESU', 'prENSUCG', 'hvap-lg'],  # enthalpy of sublimation and vaporization
        ##['prHCACP', 'prHECACOIG', 'cp-ig'], # heat capacity at constant pressure of ideal gas
        ##['prHCACP', 'prHECACOL', 'cp-l'], # heat capacity at constatn pressure of liquid
        ['prHCASP', 'prHCASPLG', 'cp-lg'],  # heat capacity at saturated pressure of liquid
        ##['prCOIE', 'prCOISEXL', 'expansion-l'], # coefficient of isobaric expansion of liquid
        ##['prIC', 'prISCOL', 'compressibility-l'], # isothermal compressibility of liquid
        ['prSOS', 'prSPSOLG', 'sound-lg'],  # speed of sound of saturated liquid
        ['prDV', 'prDYVILG', 'viscosity-lg'],  # dynamic viscosity of saturated liquid
        ##['prDV', 'prDYVIL', 'viscosity-l'], # dynamic viscosity of liquid
        ['prSUTE', 'prSURTENLG', 'st-lg'],  # surface tension
        ['prTC', 'prTHCOLG', 'tc-lg'],  # thermal conductivity of saturated liquid
        ##['prTC', 'prTHCOL', 'tc-l'], # thermal conductivity of liquid
    ]

    for prop in prop_list:
        p = NistProperty()
        p.property_id = prop[0]
        p.phase_id = prop[1]
        p.name = prop[2]

        db.session.add(p)

    db.session.commit()


def insert_mol(txt):
    '''
    Insert molecules to DB
    '''
    with open(txt) as f:
        lines = f.read().splitlines()

    for i, line in enumerate(lines):
        words = line.split('|')
        content_id, knovel_id, cas, formula, weight, name, tc, pc, dc, smiles, inchi = words[1:]
        cas = cas or None  # CAS missing
        tc = float(tc) if tc != '' else None
        pc = float(pc) if pc != '' else None
        dc = float(dc) if dc != '' else None

        try:
            m = pybel.readstring('smi', smiles)
        except:
            try:
                m = pybel.readstring('inchi', inchi)
            except:
                print('SMILES, InChI Error')
                print(line)
                continue

        if m.formula != formula:
            print('Formula Error')
            print(line)
            continue

        mol = NistMolecule.query.filter_by(content_id=content_id).first()
        if mol is None:
            mol = NistMolecule()
        mol.content_id = content_id
        mol.knovel_id = knovel_id
        mol.cas = cas  # could be None
        mol.formula = formula
        mol.name = name
        mol.tc = tc
        mol.pc = pc
        mol.dc = dc
        mol.inchi = inchi
        mol.smiles = m.write('can').strip()
        mol.weight = m.molwt
        mol.n_heavy = m.OBMol.NumHvyAtoms()

        db.session.add(mol)

    db.session.commit()


def insert_constant():
    '''Insert Constant'''

    tar = tarfile.open('/home/gongzheng/workspace/MGI/NIST/constant.tar.gz')
    names = tar.getnames()

    cids = [name.split('/')[-1].split('-')[-1] for name in names]
    cids = list(filter(lambda x: x.startswith('kc'), cids))
    print(len(cids), cids[:100])

    mols = db.session.query(NistMolecule).filter(NistMolecule.constant_inserted == False) \
        .filter(NistMolecule.content_id.in_(cids))

    print(mols.count())
    mol: NistMolecule
    for i, mol in enumerate(mols):
        if i % 10 == 0:
            sys.stdout.write('\r\t%i' % i)
            sys.stdout.flush()
        filename = mol.content_id
        for name in names:
            if name.endswith(filename):
                filename = name
                break
        else:
            print('File not exist %s' % filename)
            continue

        with tar.extractfile(tar.getmember(filename)) as f:
            j = json.load(f)
        j_constant = j['BODY']['CONSTANT_DATA']['CONSTANT']
        j_tc = j_constant.get('prCT')
        j_dc = j_constant.get('prCD')
        j_pc = j_constant.get('prCP')
        j_tb = j_constant.get('prBP')
        j_tt = j_constant.get('prTPT')
        j_hfus = j_constant.get('prEFU')
        if j_tc is not None:
            mol.tc = j_tc['srNT']['VALUES'][0]['v']
            mol.tc_u = j_tc['srNT']['VALUES'][0].get('u')
            mol.tc_has_exp = j_tc['srNT']['HAS_EXPERIMENTAL']
        if j_dc is not None:
            mol.dc = j_dc['srNT']['VALUES'][0]['v']
            mol.dc_u = j_dc['srNT']['VALUES'][0].get('u')
            mol.dc_has_exp = j_dc['srNT']['HAS_EXPERIMENTAL']
        if j_pc is not None:
            mol.pc = j_pc['srNT']['VALUES'][0]['v']
            mol.pc_u = j_pc['srNT']['VALUES'][0].get('u')
            mol.pc_has_exp = j_pc['srNT']['HAS_EXPERIMENTAL']
        if j_tb is not None:
            mol.tb = j_tb['srNT']['VALUES'][0]['v']
            mol.tb_u = j_tb['srNT']['VALUES'][0].get('u')
            mol.tb_has_exp = j_tb['srNT']['HAS_EXPERIMENTAL']
        if j_tt is not None:
            mol.tt = j_tt['srNT']['VALUES'][0]['v']
            mol.tt_u = j_tt['srNT']['VALUES'][0].get('u')
            mol.tt_has_exp = j_tt['srNT']['HAS_EXPERIMENTAL']
        if j_hfus is not None:
            mol.hfus = j_hfus['srNT']['VALUES'][0]['v']
            mol.hfus_u = j_hfus['srNT']['VALUES'][0].get('u')
            mol.hfus_has_exp = j_hfus['srNT']['HAS_EXPERIMENTAL']
        mol.constant_inserted = True

    db.session.commit()


def _get_available_props(cid):
    r = requests.post('https://app.knovel.com/api/',
                      json={
                          "APIKEY"          : "EE109DAC-CDA6-4299-8095-F6D6DBF96380",
                          "APPLICATION_NAME": "web",
                          "CLIENT_TRACK_ID" : "knovel",

                          "METHOD"          : {
                              "NAME"  : "kms-purechem/compound-thermodata-properties",
                              "PARAMS": {
                                  "CONTENT_ID": cid,
                                  "SOURCE_ID" : "all"
                              },
                              "TYPE"  : "GET"
                          },
                          "OPTIONAL_INPUT"  : {
                              "DOMAIN_NAME": "moc.levonk.ppa"
                          },
                          "SESSION"         : {
                              "ID": "5b8291cc-fbba-4dd5-c6fb-360f362d7cde"
                          }
                      }
                      )
    j = r.json()
    return j


def get_available_props():
    for nist in NistMolecule.query:
        filename = '/home/gongzheng/workspace/MGI/NIST/properties/' + nist.content_id
        if os.path.exists(filename):
            continue

        print(nist.id, nist.content_id)
        j = _get_available_props(nist.content_id)

        with open(filename, 'w') as f:
            f.write(json.dumps(j))


def insert_available_props():
    '''
    Insert available properties to DB
    The tar file contains availability of all properties
    But we only insert property information we care
    '''
    tar = tarfile.open('/home/gongzheng/workspace/MGI/NIST/properties.tar.gz')
    names = tar.getnames()

    mols = db.session.query(NistMolecule)
    properties = db.session.query(NistProperty).all()

    for i, mol in enumerate(mols):
        if i % 10 == 0:
            sys.stdout.write('\r\t%i' % i)
            sys.stdout.flush()
        filename = mol.content_id
        for name in names:
            if name.endswith(filename):
                filename = name
                break
        else:
            print('File not exist %s' % filename)
            continue

        with tar.extractfile(tar.getmember(filename)) as f:
            j = json.load(f)

        for j_prop_category in j['BODY']['COMPOUND']['PROPERTY_METADATA'].values():
            for j_prop in j_prop_category:
                for p in properties:
                    if p.phase_id not in j_prop['PHASES'].keys():
                        continue

                    has_data = NistHasData()
                    has_data.molecule = mol
                    has_data.property = p
                    has_data.property_name = p.name
                    if 'dqEXP' in j_prop['DATA_QUALITIES'].keys():
                        has_data.has_exp = True
                    if 'dqRECM' in j_prop['DATA_QUALITIES'].keys():
                        has_data.has_rec = True
                    db.session.add(has_data)
    db.session.commit()


def insert_prop_data():
    '''Insert Thermodata to DB'''

    tar = tarfile.open('/home/gongzheng/workspace/MGI/NIST/prop-data.tar.gz')
    names = tar.getnames()

    mols = db.session.query(NistMolecule).filter(NistMolecule.constant_inserted == True) \
        .filter(NistMolecule.data_inserted == False)

    print(mols.count())
    for i, mol in enumerate(mols):
        if i % 10 == 0:
            sys.stdout.write('\r\t%i' % i)
            sys.stdout.flush()

        for p in db.session.query(NistProperty):
            filename = '%s-%s-%s' % (mol.content_id, p.property_id, p.phase_id)
            for name in names:
                if name.endswith(filename):
                    filename = name
                    break
            else:
                print('File not exit %s' % filename)
                continue

            with tar.extractfile(tar.getmember(filename)) as f:
                j = json.load(f)
            data_points = j['BODY']['COMPOUND']['DATA']['DATA_POINTS']
            data_list = []
            for point in data_points:
                data_list.append([point['T'], point['P'], point['V']])
            data_list.sort(key=lambda x: (x[0], x[1]))
            for data in data_list:
                d = NistData()
                d.molecule = mol
                d.property = p
                d.t = data[0]
                d.p = data[1]
                d.value = data[2]['v']
                d.uncertainty = data[2]['u']
                db.session.add(d)

        mol.data_inserted = True

        if i % 100 == 0:
            db.session.commit()

    db.session.commit()


def fit():
    pass


def interpolate():
    '''Spline Interpolate'''
    from scipy import interpolate

    pname_list = ['pvap-lg', 'density-lg', 'density-gl', 'hvap-lg',
                  'cp-lg', 'st-lg', 'sound-lg', 'viscosity-lg', 'tc-lg']
    props = db.session.query(NistProperty).filter(NistProperty.name.in_(pname_list))

    mols = db.session.query(NistMolecule).filter(NistMolecule.data_inserted == True) \
        .filter(NistMolecule.spline_inserted == False)

    print(mols.count())
    for i, mol in enumerate(mols):
        sys.stdout.write('\r\t%i' % i)
        sys.stdout.flush()

        for p in props:
            t_list = []
            v_list = []
            u_list = []

            datas = mol.datas.filter(NistData.property == p)
            if mol.tt is not None:
                datas = datas.filter(NistData.t > mol.tt + 1)  # ignore temperatures lower than triple point + 1

            for d in datas:
                if d.value != 0:  # surface tension can be 0
                    dev = (d.uncertainty - d.value) / d.value * 100
                    if dev > 50:  # ignore points with dev large than 50 %
                        continue

                t_list.append(d.t)
                v_list.append(d.value)
                u_list.append(d.uncertainty)
            if len(t_list) <= 3:
                continue

            t_list, v_list, u_list = zip(*sorted(zip(t_list, v_list, u_list)))
            spl_v = interpolate.splrep(t_list, v_list)
            spl_u = interpolate.splrep(t_list, u_list)
            spline = NistSpline(molecule=mol, property=p)
            spline.t_min = min(t_list)
            spline.t_max = max(t_list)
            spline.coef_v = json.dumps([spl_v[0].tolist(), spl_v[1].tolist(), spl_v[2]])
            spline.coef_u = json.dumps([spl_u[0].tolist(), spl_u[1].tolist(), spl_u[2]])
            db.session.add(spline)

        mol.spline_inserted = True

        if i % 10 == 0:
            db.session.commit()

    db.session.commit()


def _organize_constant():
    olddir = '/home/gongzheng/workspace/MGI/NIST/constant_old'
    names = os.listdir(olddir)

    mols = db.session.query(NistMolecule)
    mol: NistMolecule
    i = 0
    for mol in mols:
        f = Formula(mol.formula)
        if not set(f.atomdict.keys()) <= {'C', 'H', 'O', 'N', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'Si', 'B'}:
            continue

        i += 1
        if i % 10 == 0:
            sys.stdout.write('\r\t%i' % i)
            sys.stdout.flush()

        filename = mol.content_id
        for name in names:
            if name.endswith(filename):
                filename = name
                break
        else:
            print('File not exist %s' % filename)
            continue

        newdir = '/home/gongzheng/workspace/MGI/NIST/constant'
        newname = '%i-%s' % (mol.id, mol.content_id)
        shutil.move(os.path.join(olddir, filename), os.path.join(newdir, newname))


def _organize_prop_data():
    olddir = '/home/gongzheng/workspace/MGI/NIST/prop-data_old'
    names = os.listdir(olddir)

    p_dict = dict()
    for p in db.session.query(NistProperty).all():
        p_dict[p.name] = p

    mols = db.session.query(NistMolecule)
    mol: NistMolecule
    i = 0
    for mol in mols:
        f = Formula(mol.formula)
        if not set(f.atomdict.keys()) <= {'C', 'H', 'O', 'N', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'Si', 'B'}:
            continue

        i += 1
        if i % 10 == 0:
            sys.stdout.write('\r\t%i' % i)
            sys.stdout.flush()

        for has_data in db.session.query(NistHasData).filter(NistHasData.molecule == mol).filter(
                NistHasData.has_rec == True):
            p = p_dict[has_data.property_name]

            filename = '%s-%s-%s' % (mol.content_id, p.property_id, p.phase_id)
            for name in names:
                if name.endswith(filename):
                    filename = name
                    break
            else:
                print('File not exit %s' % filename)
                continue

            newdir = '/home/gongzheng/workspace/MGI/NIST/prop-data'
            newname = '%i-%s-%s-%s' % (mol.id, mol.content_id, p.property_id, p.phase_id)
            shutil.move(os.path.join(olddir, filename), os.path.join(newdir, newname))


if __name__ == '__main__':
    fire.Fire()
