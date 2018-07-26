#!/usr/bin/env python3
# coding=utf-8

import sys, pybel
from pybel import Smarts

sys.path.append('..')
sys.path.append('../../ms-tools')

from mstools.formula import Formula

from app.models import NistMolecule, NistGroup, NistMoleculeGroup, db

smarts_bad = {
    '[#6;v0,v1,v2,v3]'                 : 'radicalC',
    '*=*=*'                            : '*=*=*',
    '[F,Cl,Br,I]~[!#6]'                : '[F,Cl,Br,I]~[!#6]',
    '*#*~[!#6]'                        : '*#*~[!#6]',
    '[NX2,NX4]'                        : '[NX2,NX4]',
    'O~[NX3](~[!$([OX1])])~[!$([OX1])]': 'O~[NX3](~[!$([OX1])])~[!$([OX1])]',
    'O~O'                              : 'peroxide',
    'N~N'                              : 'N~N',
    '[O,N]*[O,N;H1,H2]'                : '[O,N]*[O,N;H1,H2]',
    'C=C~[O,N;H1,H2]'                  : 'C=C~[O,N;H1,H2]',
    'O=C~*~C=O'                        : 'beta-dicarbonyl',
    'C1=CC=C*1'                        : 'cyclopentadiene',
    'a=*'                              : 'a=*',
    '[nX3;r6]'                         : 'pyridine-N-oxide',
    '[R3]'                             : '[R3]',
    '[r3,r4;R2]'                       : '[r3,r4;R2]',
    '[r3,r4;#6X3]'                     : '[r3,r4;#6X3]',
    '[r3,r4,r5]~[F,Cl,Br,I]'           : '[r3,r4,r5]~[F,Cl,Br,I]',
    '[r3,r4,r5]~[#7,#8]'               : '[r3,r4,r5]~[#7,#8]',
}
smarts_good = {
    'n'                          : 'N-aromatic',
    'C=C~[O,N;H0]'               : 'C=C~[O,N;H0]',
    '[C,c]O[NX3](~[OX1])~[OX1]'  : 'nitrate',

    '[r3]'                       : 'ring3',
    '[r4]'                       : 'ring4',
    '[r5]'                       : 'ring5',
    '[CX3]=[CX3]'                : 'alkene',
    '[CX2]#[CX2]'                : 'alkyne',
    'c1ccccc1'                   : 'benzene',
    'c1ccc2c(c1)cccc2'           : 'polybenzene',
    # 'COC'                       : 'ether',
    # 'cO[C,c]'                   : 'Ph-ether',
    '[C,c]O[C,c]'                : 'ether',
    '[O;H0]C[O;H0]'              : 'acetal',
    '[C,c]C(=O)[C,c]'            : 'ketone',
    '[C,c]C(=O)[H]'              : 'aldehyde',
    # 'C[OH]'                     : 'alcohol',
    # 'c[OH]'                     : 'Ph-ol',
    '[C,c][OH]'                  : 'alcohol',
    'C(=O)O[C,c;!$(C=O)]'        : 'ester',
    'C(=O)[OH]'                  : 'acid',
    # 'C[NX3;!$(NC=O);!$(N~O);!$(Nc)]': 'amine',
    # 'c[NX3;!$(NC=O);!$(N~O)]'       : 'Ph-amine',
    '[C,c][NX3;!$(NC=O);!$(N~O)]': 'amine',
    'C(=O)[NX3;!$(NO)]'          : 'amide',
    '[NX1]#[CX2][C,c]'           : 'nitrile',
    '[C,c][NX3](~[OX1])~[OX1]'   : 'nitro',
    '[F,Cl,Br,I]'                : 'halogen',
    '[CX4][F,Cl,Br,I]'           : 'CX4-X',
    'c[F,Cl,Br,I]'               : 'c-X',
    'C=C[F,Cl,Br,I]'             : 'C=C-X',
    'O=C[F,Cl,Br,I]'             : 'O=C-X',
    'F'                          : 'F',
    'Cl'                         : 'Cl',
    'Br'                         : 'Br',
    'I'                          : 'I',
}
smarts_dict = dict(**smarts_bad, **smarts_good)

group_dict = {}


def init_group():
    db.create_all()
    db.session.query(NistMoleculeGroup).delete()
    db.session.query(NistGroup).delete()

    for smarts, name in smarts_dict.items():
        group = NistGroup(name=name, smarts=smarts)
        if smarts in smarts_bad.keys():
            group.bad = True
        group_dict[name] = group
        db.session.add(group)

    diol = NistGroup(name='diol', smarts='diol')
    group_dict['diol'] = diol
    db.session.add(diol)

    CH = NistGroup(name='hydrocarbon', smarts='hydrocarbon')
    group_dict['hydrocarbon'] = CH
    db.session.add(CH)

    alkane = NistGroup(name='alkane', smarts='alkane')
    group_dict['alkane'] = alkane
    db.session.add(alkane)

    db.session.commit()


def match_group():
    print('Match group from smiles...')
    for molecule in db.session.query(NistMolecule).filter(NistMolecule.smiles != None):
        if molecule.id % 100 == 0:
            sys.stdout.write('\r\t%i' % molecule.id)
            sys.stdout.flush()
        if molecule.smiles.find('.') > -1:  # ignore complex
            continue
        py_mol = pybel.readstring('smi', molecule.smiles)
        formula = py_mol.formula
        atom_set = set(Formula.read(formula).atomdict.keys())
        if not ({'C', 'H'} <= atom_set and atom_set <= {'C', 'H', 'O', 'N', 'F', 'Cl', 'Br', 'I'}):
            continue

        if atom_set == {'C', 'H'}:
            CH = group_dict['hydrocarbon']
            mol_group = NistMoleculeGroup(molecule_id=molecule.id, group_id=CH.id)
            db.session.add(mol_group)

            for s in ['[CX2]', '[CX3]', 'c', '[#6;v0,v1,v2,v3]']:
                if pybel.Smarts(s).findall(py_mol) != []:
                    break
            else:
                alkane = group_dict['alkane']
                mol_group = NistMoleculeGroup(molecule_id=molecule.id, group_id=alkane.id)
            db.session.add(mol_group)

        if Smarts('[OH]').findall(py_mol).__len__() > 1:
            diol = group_dict['diol']
            mol_group = NistMoleculeGroup(molecule_id=molecule.id, group_id=diol.id)
            db.session.add(mol_group)

        for smarts, name in smarts_dict.items():
            if Smarts(smarts).findall(py_mol):
                group = group_dict[name]
                mol_group = NistMoleculeGroup(molecule_id=molecule.id, group_id=group.id)
                db.session.add(mol_group)

    db.session.commit()


def add_remark():
    NistMolecule.query.update({'remark': None})
    db.session.commit()
    return None

    molecules = NistMolecule.query.filter(NistMolecule.n_heavy <= 19)
    #####
    print('Select CHON')
    mol_list = molecules.all()
    for name in smarts_bad.values():
        exclude_molecules = molecules.filter(NistMolecule.groups.any(NistGroup.name == name)).all()
        mol_list = [mol for mol in mol_list if mol not in exclude_molecules]

    mol_list = filter(lambda x: len(x.groups) >= 1, mol_list)
    for mol in mol_list:
        mol.remark = 'CHON'

    #####
    print('Select alkane')
    mol_list = molecules.filter(NistMolecule.groups.any(NistGroup.name == 'alkane')).all()
    for mol in mol_list:
        mol.remark = 'alkane'

    #####
    print('Select amide')
    mol_list = molecules.filter(NistMolecule.groups.any(NistGroup.name == 'amide')).all()
    for mol in mol_list:
        mol.remark = 'amide'

    db.session.commit()


if __name__ == '__main__':
    init_group()
    match_group()
    add_remark()
