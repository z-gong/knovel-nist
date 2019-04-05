#!/usr/bin/env python3

import os
import sys
import pybel
import fire

sys.path.append('..')
sys.path.append('../../ms-tools')

from mstools.formula import Formula


def select_molecules(sdf_dir, out_dir, begin=0, end=1000_000):
    '''
    Filter molecules in sdf.gz files by charge, n_heavy, element, components
    Then write the selected molecules in sdf files

    :param sdf_dir:
    :param out_dir:
    :param begin:
    :param end:
    :return:
    '''
    sdf_list = list(filter(lambda x: x.endswith('.sdf.gz'), os.listdir(sdf_dir)))
    print(len(sdf_list))

    for i, sdf in enumerate(sdf_list[begin: end]):
        sys.stdout.write('\r\t%i / %i' % (i + begin, len(sdf_list)))
        sys.stdout.flush()
        sdf_out = pybel.Outputfile('sdf', os.path.join(out_dir, 'CHONFClBr-%04i.sdf' % (i + 1)))

        for m in pybel.readfile('sdf', os.path.join(sdf_dir, sdf), opt={'P': None}):
            try:
                cid = int(m.data['PUBCHEM_COMPOUND_CID'])
                formula = m.data['PUBCHEM_MOLECULAR_FORMULA']
                name = m.data['PUBCHEM_IUPAC_NAME']
                smiles = m.data['PUBCHEM_OPENEYE_ISO_SMILES']
                inchi = m.data['PUBCHEM_IUPAC_INCHI']
                cactvs = m.data['PUBCHEM_CACTVS_SUBSKEYS']  # base64 encoded
                weight = float(m.data['PUBCHEM_MOLECULAR_WEIGHT'])
                charge = int(m.data['PUBCHEM_TOTAL_CHARGE'])
                n_heavy = int(m.data['PUBCHEM_HEAVY_ATOM_COUNT'])
            except:
                continue

            # Ignore ion
            if charge != 0:
                continue

            f = Formula(formula)
            # Ignore large molecule
            if f.n_heavy > 19:
                continue
            # Limit element
            atom_set = set(f.atomdict.keys())
            if 'C' not in atom_set or atom_set & {'H', 'F', 'Cl', 'Br'} == set() or \
                    not atom_set <= {'C', 'H', 'O', 'N', 'F', 'Cl', 'Br'}:
                continue
            # Kick out mixture
            if smiles.find('.') > -1:
                continue

            mol = pybel.readstring('smi', smiles)
            if mol.formula != formula or mol.charge != charge:
                print('SMILES formula Error:', cid)
                continue

            sdf_out.write(m)
            continue

        sdf_out.close()


def save_smiles(sdf_dir, smi_out, begin=0, end=1000_000):
    '''
    Save the SMILES of molecules in sdf files into single smi file

    :param sdf_dir:
    :param smi_out:
    :param begin:
    :param end:
    :return:
    '''
    sdf_list = list(filter(lambda x: x.endswith('.sdf'), os.listdir(sdf_dir)))
    print(len(sdf_list))

    f_out = open(smi_out, 'w')
    for i, sdf in enumerate(sdf_list[begin: end]):
        sys.stdout.write('\r\t%i / %i' % (i + begin, len(sdf_list)))
        sys.stdout.flush()

        for m in pybel.readfile('sdf', os.path.join(sdf_dir, sdf), opt={'P': None}):
            smiles = m.data['PUBCHEM_OPENEYE_ISO_SMILES']
            f_out.write('%s\t%s\n' % (smiles, m.title))

    f_out.close()


def select_unique(sim, threshold):
    '''
    Select unique molecules using the similarity file produced by chemfp

    :param sim:
    :param threshold:
    :return:
    '''
    bad_set = set()

    f_out = open('unique.txt', 'w')
    N = 0
    for line in open(sim):
        if line.strip() == '' or line.startswith('#'):
            continue

        N += 1
        if N % 1000 == 0:
            sys.stdout.write('\r\t%i' % N)
        words = line.strip().split()
        good = int(words[1])
        if good in bad_set:
            continue

        for i in range((len(words) - 2) // 2):
            bad = int(words[2 + 2 * i])
            score = float(words[2 + 2 * i + 1])
            if bad > good and score >= threshold:
                bad_set.add(bad)

        f_out.write('%i\n' % good)


def select_fps(fps, ids):
    '''
    Save a subset of fps from id listed in ids file

    :param fps:
    :param ids:
    :return:
    '''
    fps_list = [''] * 150_000_000

    f_out = open('select.fps', 'w')

    N = 0
    for line in open(fps):
        if line.startswith('#'):
            f_out.write(line)
            continue

        N += 1
        if N % 1000 == 0:
            sys.stdout.write('\r\t%i' % N)
        words = line.strip().split()
        fps_list[int(words[1])] = words[0]

    with open(ids) as f:
        for cid in f.read().splitlines():
            f_out.write('%s\t%s\n' % (fps_list[int(cid)], cid))


def insert_db(sdf_dir, unique_cids, begin=0, end=1000_000):
    '''
    Write selected molecules from sdf files into database

    :param sdf_dir:
    :param unique_cids:
    :param begin:
    :param end:
    :return:
    '''
    from app.models_pubchem import PubchemMolecule, db

    sdf_list = list(filter(lambda x: x.endswith('.sdf'), os.listdir(sdf_dir)))
    print(len(sdf_list))

    with open(unique_cids) as f:
        cids = list(map(int, f.read().splitlines()))

    for i, sdf in enumerate(sdf_list[begin: end]):
        sys.stdout.write('\r\t%i / %i' % (i + begin, len(sdf_list)))
        sys.stdout.flush()

        m_list = [None] * 150_000_000
        for m in pybel.readfile('sdf', os.path.join(sdf_dir, sdf), opt={'P': None}):
            cid = int(m.title)
            m_list[cid] = m

        for cid in cids:
            m = m_list[cid]
            if m is None:
                continue

            formula = m.data['PUBCHEM_MOLECULAR_FORMULA']
            name = m.data['PUBCHEM_IUPAC_NAME']
            smiles = m.data['PUBCHEM_OPENEYE_ISO_SMILES']
            # inchi = m.data['PUBCHEM_IUPAC_INCHI']
            n_heavy = int(m.data['PUBCHEM_HEAVY_ATOM_COUNT'])

            mol = pybel.readstring('smi', smiles)

            pubchem = PubchemMolecule(cid=cid,
                                      formula=formula,
                                      name=name,
                                      smiles=mol.write('can').strip(),
                                      # inchi=inchi,
                                      weight=mol.molwt,
                                      n_heavy=n_heavy,
                                      )

            db.session.add(pubchem)

        db.session.commit()


def print_mols_for_msdweb():
    '''
    Write molecule information into a txt file, so that can be imported into msdweb
    Information needed are SMILES and IUPAC
    '''
    from app.models_pubchem import PubchemMolecule, db

    print('#No\tCID\tFormula\tsmiles\tiupac')
    ### Only select molecules with n_heavy >= 4 and n_C >= 2
    ### Ignore chiral isomers
    mols = db.session.query(PubchemMolecule) \
        .filter(PubchemMolecule.remark == None) \
        .filter(PubchemMolecule.smiles.notilike('%@%')) \
        .filter(PubchemMolecule.n_heavy >= 4)

    i = 0
    for mol in mols:
        f = Formula.read(mol.formula)
        if f.atomdict.get('C') < 2:
            continue

        i += 1
        print(i, mol.cid, mol.formula, mol.smiles, mol.name, sep='\t')


if __name__ == '__main__':
    fire.Fire()
