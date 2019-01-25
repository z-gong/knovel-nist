#!/usr/bin/env python3
# coding=utf-8

import sys, pybel
from pybel import Smarts

sys.path.append('..')
sys.path.append('../../ms-tools')

from mstools.formula import Formula

from app.models import db

if sys.argv[1] == 'nist':
    from app.models import NistMolecule as Molecule, NistGroup as Group, NistMoleculeGroup as MoleculeGroup
elif sys.argv[1] == 'pubchem':
    from app.models_pubchem import PubchemMolecule as Molecule, PubchemGroup as Group, \
        PubchemMoleculeGroup as MoleculeGroup

smarts_bad = {
    '[#6;v0,v1,v2,v3]'             : 'radicalC',
    '*=*=*'                        : '*=*=*',
    '*#*~*#*'                      : '*#*~*#*',
    '[F,Cl,Br]~[!#6]'              : '[F,Cl,Br]~[!#6]',
    '*#*~[!#6]'                    : '*#*~[!#6]',
    '[NX2,NX4]'                    : '[NX2,NX4]',
    'O~N(~[!$([OX1])])~[!$([OX1])]': 'O~N(~[!$([OX1])])~[!$([OX1])]',
    'O~O'                          : 'peroxide',
    'N~N'                          : 'N~N',
    '[O,N]*[O,N;H1,H2]'            : '[O,N]*[O,N;H1,H2]',
    'C=C~[O,N;H1,H2]'              : 'C=C~[O,N;H1,H2]',
    'O=C~*~C=O'                    : 'beta-dicarbonyl',

    'a=*'                          : 'a=*',
    'o'                            : 'o',
    '[n;r5]'                       : '[n;r5]',
    '[nX3;r6]'                     : 'pyridine-N-oxide',
    '[$(nnn),$(nnan),$(nanan)]'    : 'triazine(zole)',

    '[R3]'                         : '[R3]',
    '[r3,r4;R2]'                   : '[r3,r4;R2]',
    '[r3,r4;#6X3]'                 : '[r3,r4;#6X3]',
    '[r3,r4]~[!#6]'                : '[r3,r4]~[!#6]',

    'O[NX3](~[OX1])~[OX1]'         : 'nitrate',
    'O=C[NX3]'                     : 'amide',
    'O=C[F,Cl,Br]'                 : 'acyl-halide',
    'c1ccc2c(c1)cccc2'             : 'polybenzene',
}
smarts_good = {
    'C=C~[O,N;H0]'               : 'C=C~[O,N;H0]',
    '[r5]~[!#6]'                 : '[r5]~[!#6]',
    '[r3]'                       : 'ring3',
    '[r4]'                       : 'ring4',
    '[r5]'                       : 'ring5',
    '[CX3]=[CX3]'                : 'alkene',
    '[CX2]#[CX2]'                : 'alkyne',
    'c1ccccc1'                   : 'benzene',
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
    '[NX1]#[CX2][C,c]'           : 'nitrile',
    '[C,c][NX3](~[OX1])~[OX1]'   : 'nitro',
    'n'                          : 'N-aromatic',
    '[F,Cl,Br]'                  : 'halogen',
    # 'F'                          : 'F',
    # 'Cl'                         : 'Cl',
    # 'Br'                         : 'Br',
    '[CX4][F,Cl,Br]'             : 'CX4-X',
    'c[F,Cl,Br]'                 : 'c-X',
    'C=C[F,Cl,Br]'               : 'C=C-X',
}
smarts_dict = dict(**smarts_bad, **smarts_good)

group_dict = {}


def init_group():
    db.create_all()
    db.session.query(MoleculeGroup).delete()
    db.session.query(Group).delete()

    for smarts, name in smarts_dict.items():
        group = Group(name=name, smarts=smarts)
        if smarts in smarts_bad.keys():
            group.bad = True
        group_dict[name] = group
        db.session.add(group)

    diol = Group(name='diol', smarts='diol')
    group_dict['diol'] = diol
    db.session.add(diol)

    CH = Group(name='hydrocarbon', smarts='hydrocarbon')
    group_dict['hydrocarbon'] = CH
    db.session.add(CH)

    alkane = Group(name='alkane', smarts='alkane')
    group_dict['alkane'] = alkane
    db.session.add(alkane)

    db.session.commit()


def match_group():
    print('Match group from smiles...')
    for molecule in db.session.query(Molecule).filter(Molecule.smiles != None):
        if molecule.id % 100 == 0:
            sys.stdout.write('\r\t%i' % molecule.id)
            sys.stdout.flush()
        if molecule.smiles.find('.') > -1:  # ignore complex
            continue
        py_mol = pybel.readstring('smi', molecule.smiles)
        py_mol.removeh()
        formula = py_mol.formula
        atom_set = set(Formula.read(formula).atomdict.keys())
        if not atom_set <= {'C', 'H', 'O', 'N', 'F', 'Cl', 'Br'}:
            continue
        if not ('C' in atom_set and {'H', 'F', 'Cl', 'Br'} & atom_set != set()):
            continue

        if atom_set == {'C', 'H'}:
            CH = group_dict['hydrocarbon']
            mol_group = MoleculeGroup(molecule_id=molecule.id, group_id=CH.id)
            db.session.add(mol_group)

            for s in ['[CX2]', '[CX3]', 'c', '[#6;v0,v1,v2,v3]']:
                if pybel.Smarts(s).findall(py_mol) != []:
                    break
            else:
                alkane = group_dict['alkane']
                mol_group = MoleculeGroup(molecule_id=molecule.id, group_id=alkane.id)
            db.session.add(mol_group)

        if Smarts('[OH]').findall(py_mol).__len__() > 1:
            diol = group_dict['diol']
            mol_group = MoleculeGroup(molecule_id=molecule.id, group_id=diol.id)
            db.session.add(mol_group)

        for smarts, name in smarts_dict.items():
            if Smarts(smarts).findall(py_mol):
                group = group_dict[name]
                mol_group = MoleculeGroup(molecule_id=molecule.id, group_id=group.id)
                db.session.add(mol_group)

    db.session.commit()


if __name__ == '__main__':
    init_group()
    match_group()
